#include "VariationalFormsCarotidIso.h"
#include "HyperElasticityCarotidIso.h"
#include "P1VectorElement.h"
using namespace dolfin;

// Sub domain for clamp at left end
class InletXaxisPoint : public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        if( (std::abs(x[2] - 0.0085466e3 ) < 1.e-5) &&
                (std::abs(x[1] + 0.00122233e3) < 1.e-5) 
          )
        {   
            //std::cout<<"found X point"<<std::endl;
            //std::cout<<x[0]<<"  "<<x[1]<<"  "<<x[2]<<std::endl;
            return true;
        }
        return false;
    }
};
class InletYaxisPoint : public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        if( (std::abs(x[2] - 0.0085466e3 ) < 1.e-5) &&
                (std::abs(x[0] - 0.0272064e3) < 1.e-5) 
          )
        {   
            //std::cout<<"found Y point"<<std::endl;
            //std::cout<<x[0]<<"  "<<x[1]<<"  "<<x[2]<<std::endl;
            return true;
        }
        return false;
    }
};

// Pressure boundary condition
class PressureNormal : public Expression
{
    public:

        PressureNormal(const Mesh& mesh, double* pres) : Expression(3), mesh(mesh), pressure(pres){}

        void eval(Array<double>& values, const Array<double>& x,
                const ufc::cell& ufc_cell) const
        {
            dolfin_assert(ufc_cell.local_facet >= 0);

            Cell cell(mesh, ufc_cell.index);
            Point n = cell.normal(ufc_cell.local_facet);
            values[0] = -*pressure*n[0];
            values[1] = -*pressure*n[1];
            values[2] = -*pressure*n[2];
        }

    private:

        const Mesh& mesh;
        double* pressure;


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


VariationalFormsCarotidIso::VariationalFormsCarotidIso(
        Mesh& mesh, MeshFunction<std::size_t>& sub_domains_mark,
        MeshFunction<std::size_t>& boundary_mark,
        Parameters& _para): para(_para)
{

    _mesh = reference_to_no_delete_pointer(mesh);
    std::string elas_model = para["elas_model"];
    switch (elas_model[0]) {
        case 'C':
            info("inital forms of model CarotidIso");
            _V = std::make_shared<HyperElasticityCarotidIso::Form_Res::TestSpace>(_mesh);
            _u = std::make_shared<Function>(_V);
            _F = std::make_shared<HyperElasticityCarotidIso::Form_Res>(_V);
            _obj = std::make_shared<HyperElasticityCarotidIso::Form_Pi>(_mesh);
            dolfin_assert(_F);
            _J = std::make_shared<HyperElasticityCarotidIso::Form_Jac>(_V, _V);
            _VMS = std::make_shared<HyperElasticityCarotidIso::Form_L_vms::TestSpace>(_mesh);
            _p = std::make_shared<Function>(_VMS);
            _a = std::make_shared<HyperElasticityCarotidIso::Form_Mass_vms>(_VMS,_VMS);
            _L = std::make_shared<HyperElasticityCarotidIso::Form_L_vms>(_VMS);
            _V_p1 = std::make_shared<P1VectorElement::FunctionSpace>(_mesh);
            _u_p1 = std::make_shared<Function>(_V_p1);
            break;
        default :
            dolfin_error("VariationalFormsCarotidIso.cpp",
                    "Switch Elasticity Model",
                    "Error option for Elasticity Model");
    }
    info("# of Vertices %d, # of Cells %d, # of Unknowns %d", mesh.num_vertices(),
            mesh.num_cells(), _V->dim());
    // Set material parameters
    double C1 = (double)(para["C1"]),
           Epsilon1  = (double)(para["Epsilon1"]),
           Epsilon2= (double)(para["Epsilon2"]),
           Beta1 = (double)(para["Beta1"]),
           Eta1 = (double)(para["Eta1"]),
           Delta1 = (double)(para["Delta1"]);
    double Delta2 = Beta1+2*Eta1+Delta1;

    std::shared_ptr<Constant> coef_c1 = std::make_shared<Constant>(C1);
    std::shared_ptr<Constant> coef_epsilon1 = std::make_shared<Constant>(Epsilon1);
    std::shared_ptr<Constant> coef_epsilon2 = std::make_shared<Constant>(Epsilon2);
    std::shared_ptr<Constant> coef_beta1 = std::make_shared<Constant>(Beta1);
    std::shared_ptr<Constant> coef_eta1 = std::make_shared<Constant>(Eta1);
    std::shared_ptr<Constant> coef_delta1 = std::make_shared<Constant>(Delta1);
    std::shared_ptr<Constant> coef_delta2 = std::make_shared<Constant>(Delta2);

    // Define source and boundary traction functions
    pressure = para["pressure_boundary_condition"]; t = 1.0;
    final_pressure = pressure;
    std::shared_ptr<PressureNormal> pressure_normal = 
        std::make_shared<PressureNormal>(mesh,&pressure);
    std::shared_ptr<Constant> B = std::make_shared<Constant>(0.0, 0.0, 0.0);

    std::map<std::string, std::shared_ptr<const GenericFunction>> coef_list
        = { {"u",       _u},
            {"c1",      coef_c1},
            {"epsilon1",coef_epsilon1},
            {"epsilon2",coef_epsilon2},
            {"beta1",      coef_beta1},
            {"eta1",      coef_eta1},
            {"delta1",      coef_delta1},
            {"delta2",      coef_delta2},
            {"B",       B},
            {"T",       pressure_normal}};
    info("number of coefficients %d", coef_list.size());


    _F->set_some_coefficients(coef_list);
    _F->ds = reference_to_no_delete_pointer(boundary_mark);
    _F->dx = reference_to_no_delete_pointer(sub_domains_mark);
    _obj->set_some_coefficients(coef_list);
    _obj->ds = reference_to_no_delete_pointer(boundary_mark);
    _obj->dx = reference_to_no_delete_pointer(sub_domains_mark);
    _J->set_some_coefficients(coef_list);
    _J->dx = reference_to_no_delete_pointer(sub_domains_mark);
    _L->set_some_coefficients(coef_list);
    _L->dx = reference_to_no_delete_pointer(sub_domains_mark);

    // Define Dirichlet boundaries
    //std::shared_ptr<LeftSide> left = std::make_shared<LeftSide>();
    //std::shared_ptr<RightSide> right = std::make_shared<RightSide>();

    // Define Dirichlet boundary functions
    std::shared_ptr<Constant> zero = std::make_shared<Constant>(0.0);
    std::string  method("pointwise");


    std::shared_ptr<InletXaxisPoint> inletXPoint =std::make_shared<InletXaxisPoint>();
    std::shared_ptr<InletYaxisPoint> inletYPoint =std::make_shared<InletYaxisPoint>();



    // Create Dirichlet boundary conditions
    std::shared_ptr<DirichletBC> bcInlet = std::make_shared<DirichletBC>((*_V)[2], zero,
            reference_to_no_delete_pointer(boundary_mark), 1);
    std::shared_ptr<DirichletBC> bcOutlet = std::make_shared<DirichletBC>((*_V)[2], zero,
            reference_to_no_delete_pointer(boundary_mark), 2);
    std::shared_ptr<DirichletBC> bcx = std::make_shared<DirichletBC>((*_V)[0], zero, inletXPoint, method);
    std::shared_ptr<DirichletBC> bcy = std::make_shared<DirichletBC>((*_V)[1], zero, inletYPoint, method);
    bcs.push_back(bcInlet);
    bcs.push_back(bcOutlet);
    bcs.push_back(bcx);
    bcs.push_back(bcy);

}


void VariationalFormsCarotidIso::save_solution()
{
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << double(para["pressure_boundary_condition"]);
    std::string idx = std::string("Pres") + ss.str();
    // Save solution in VTK format
    File file(std::string("./result/displacement")+idx+std::string(".pvd"));
    file << *_u;

    std::string str(std::string("backup_solution")+(".xml"));
    backup_solution(str);
}

void VariationalFormsCarotidIso::backup_solution(std::string str)
{
    File file_backup(str);
    file_backup <<*_u;
}

void VariationalFormsCarotidIso::load_solution(std::string str)
{
    if(! File::exists(str)) return;
    File file_backup(str);
    file_backup >> *_u;
}

void VariationalFormsCarotidIso::save_von_misec_stress()
{
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << double(para["pressure_boundary_condition"]);
    std::string idx = std::string("Pres") + ss.str();

    solve(*_a==*_L,*_p);
    if (_V->element()->ufc_element()->degree() != 1)
    {
        _u_p1->interpolate(*_u);
        ALE::move(*_mesh, *_u_p1);
    }
    else ALE::move(*_mesh, *_u);
    File filep(std::string("./result/stress")+idx+std::string(".pvd"));
    filep << *_p;
}

void VariationalFormsCarotidIso::plot_solution()
{
    plot(*_u);
    plot(*_p);
    interactive();

}


