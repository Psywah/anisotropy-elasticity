#include "VariationalFormsIsoElas.h"
#include "HyperElasticityBIso.h"
using namespace dolfin;
#define Length 2.0
#define R 0.5

// Sub domain for clamp at left side
class LeftSide : public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        if(std::abs(x[0]) < DOLFIN_EPS)
            return true;
        return false;
    }
};

// Sub domain for enforcing at right side
class RightSide : public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        if (std::abs(x[0]-Length) < DOLFIN_EPS)
            return true;
        return false;
    }
};


// Pressure boundary condition
class Traction : public Expression
{
    public:

        Traction(const Mesh& mesh, const double* T, const char direction) 
            : Expression(3), mesh(mesh), t(T), d(direction){}

        void eval(Array<double>& values, const Array<double>& x,
                const ufc::cell& ufc_cell) const
        {
            //dolfin_assert(ufc_cell.local_facet >= 0);

            //Cell cell(mesh, ufc_cell.index);
            //Point n = cell.normal(ufc_cell.local_facet);
            values[0] = 0.;
            values[1] = 0.;
            values[2] = 0.;

            if( d =='x')
                values[0] = -*t;
            else if(d == 'y')
                values[1] = -*t;
            else if(d == 'z')
                values[2] = -*t;
            else {
                double r = std::sqrt(x[1]*x[1]+x[2]*x[2])+DOLFIN_EPS;
                values[0]=0.0;
                values[1] = -x[2]/R**t;
                values[2] = x[1]/R**t;
            }
        }

    private:

        const Mesh& mesh;
        const double* t;
        const char d;
};

//class MaterialCoef : public Expression
//{
//    private:
//        double adv, med, pla;
//    public:
//
//    MaterialCoef(double adv, double med, double pla) : Expression(), adv(adv), med(med), pla(pla){}
//
//        void eval(Array<double>& values, const Array<double>& x) const
//        {
//            double r = sqrt(x[0]*x[0]+x[1]*x[1]);
//            if ( r > R+Thk_med)
//                values[0] = adv;
//            else if (r > R)
//                values[0] = med;
//            else
//                values[0] = pla;
//        }
//
//};
//




class Clamp : public Expression
{
    public:

        Clamp() : Expression(3) {}

        void eval(Array<double>& values, const Array<double>& x) const
        {
            values[0] = 0.0;
            values[1] = 0.0;
            values[2] = 0.0;
        }

};


VariationalFormsIsoElas::VariationalFormsIsoElas(
        Mesh& mesh, MeshFunction<std::size_t>& sub_domains_mark,
        MeshFunction<std::size_t>& boundary_mark,
        Parameters& _para): para(_para)
{
    
    _mesh = reference_to_no_delete_pointer(mesh);
    std::string elas_model = para["elas_model"];
    switch (elas_model[0]) {
        case 'B':
            info("inital forms of model B_iso");
            _V = std::make_shared<HyperElasticityBIso::Form_Res::TestSpace>(mesh);
            _u = std::make_shared<Function>(_V);
            _F = std::make_shared<HyperElasticityBIso::Form_Res>(_V);
            dolfin_assert(_F);
            _J = std::make_shared<HyperElasticityBIso::Form_Jac>(_V, _V);
            _VMS = std::make_shared<HyperElasticityBIso::Form_L_vms::TestSpace>(mesh);
            _p = std::make_shared<Function>(_VMS);
            _a = std::make_shared<HyperElasticityBIso::Form_Mass_vms>(_VMS,_VMS);
            _L = std::make_shared<HyperElasticityBIso::Form_L_vms>(_VMS);
            break;
        default :
            dolfin_error("VariationalFormsIsoElas.cpp",
                         "Switch Elasticity Model",
                         "Error option for Elasticity Model");
    }
    info("# of Vertices %d, # of Cells %d, # of Unknowns %d", mesh.num_vertices(),
         mesh.num_cells(), _V->dim());
    // Set material parameters
    double C1 = (double)(para["C1"]),
           Epsilon1  = (double)(para["Epsilon1"]),
           Epsilon2= (double)(para["Epsilon2"]);

    std::shared_ptr<Constant> coef_c1 = std::make_shared<Constant>(C1);
    std::shared_ptr<Constant> coef_epsilon1 = std::make_shared<Constant>(Epsilon1);
    std::shared_ptr<Constant> coef_epsilon2 = std::make_shared<Constant>(Epsilon2);

    // Define source and boundary traction functions
    pressure = para["pressure_boundary_condition"]; t = 1.0;
    final_pressure = pressure;
    std::string traction_direction = para["traction_direction"];
    std::shared_ptr<Traction> traction = 
                     std::make_shared<Traction>(mesh,&pressure, traction_direction[0]);
    std::shared_ptr<Constant> B = std::make_shared<Constant>(0.0, 0.0, 0.0);
     
    std::map<std::string, std::shared_ptr<const GenericFunction>> coef_list
        = { {"u",       _u},
            {"c1",      coef_c1},
            {"epsilon1",coef_epsilon1},
            {"epsilon2",coef_epsilon2},
            {"B",       B},
            {"T",       traction}};
    info("number of coefficients %d", coef_list.size());


    _F->set_some_coefficients(coef_list);
    _F->ds = boundary_mark;
    _F->dx = sub_domains_mark;
    _J->set_some_coefficients(coef_list);
    _J->dx = sub_domains_mark;
    _L->set_some_coefficients(coef_list);
    _L->dx = sub_domains_mark;

    // Define Dirichlet boundaries
    std::shared_ptr<LeftSide> left = std::make_shared<LeftSide>();
    std::shared_ptr<RightSide> right = std::make_shared<RightSide>();

    // Define Dirichlet boundary functions
    std::shared_ptr<Constant> zero = std::make_shared<Constant>(0.0,0.0,0.0);
    std::string  method("pointwise");

    // Create Dirichlet boundary conditions
    std::shared_ptr<DirichletBC> bcl = std::make_shared<DirichletBC>((_V), zero, left, method);
    bcs.push_back(bcl);

}

    
void VariationalFormsIsoElas::save_solution()
{
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << double(para["pressure_boundary_condition"]);
    std::string idx = std::string("Pres") + ss.str();
    // Save solution in VTK format
    File file(std::string("displacement")+idx+std::string(".pvd"));
    file << *_u;
    File file_backup(std::string("backup_solution")+(".xml"));
    file_backup << *_u;
}

void VariationalFormsIsoElas::load_solution(std::string str)
{
    File file_backup(str);
    file_backup >> *_u;
}

void VariationalFormsIsoElas::save_von_misec_stress()
{
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << double(para["pressure_boundary_condition"]);
    std::string idx = std::string("Pres") + ss.str();

    solve(*_a==*_L,*_p);
    _mesh->move(*_u);
    File filep(std::string("stress")+idx+std::string(".pvd"));
    filep << *_p;
}

void VariationalFormsIsoElas::plot_solution()
{
    plot(*_u);
    plot(*_p);
    interactive();

}


