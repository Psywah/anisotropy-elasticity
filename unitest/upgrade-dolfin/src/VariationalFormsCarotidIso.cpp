#include "VariationalFormsCarotidIso.h"
#include "HyperElasticityBIso.h"
using namespace dolfin;

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
        case 'B':
            info("inital forms of model B_iso");
            _V = std::make_shared<HyperElasticityBIso::Form_Res::TestSpace>(_mesh);
            _u = std::make_shared<Function>(_V);
            _F = std::make_shared<HyperElasticityBIso::Form_Res>(_V);
            dolfin_assert(_F);
            _J = std::make_shared<HyperElasticityBIso::Form_Jac>(_V, _V);
            _VMS = std::make_shared<HyperElasticityBIso::Form_L_vms::TestSpace>(_mesh);
            _p = std::make_shared<Function>(_VMS);
            _a = std::make_shared<HyperElasticityBIso::Form_Mass_vms>(_VMS,_VMS);
            _L = std::make_shared<HyperElasticityBIso::Form_L_vms>(_VMS);
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
           Epsilon2= (double)(para["Epsilon2"]);

    std::shared_ptr<Constant> coef_c1 = std::make_shared<Constant>(C1);
    std::shared_ptr<Constant> coef_epsilon1 = std::make_shared<Constant>(Epsilon1);
    std::shared_ptr<Constant> coef_epsilon2 = std::make_shared<Constant>(Epsilon2);

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
            {"B",       B},
            {"T",       pressure_normal}};
    info("number of coefficients %d", coef_list.size());


    _F->set_some_coefficients(coef_list);
    _F->ds = reference_to_no_delete_pointer(boundary_mark);
    _F->dx = reference_to_no_delete_pointer(sub_domains_mark);
    _J->set_some_coefficients(coef_list);
    _J->dx = reference_to_no_delete_pointer(sub_domains_mark);
    _L->set_some_coefficients(coef_list);
    _L->dx = reference_to_no_delete_pointer(sub_domains_mark);

    // Define Dirichlet boundaries
    //std::shared_ptr<LeftSide> left = std::make_shared<LeftSide>();
    //std::shared_ptr<RightSide> right = std::make_shared<RightSide>();

    // Define Dirichlet boundary functions
    std::shared_ptr<Constant> zero = std::make_shared<Constant>(0.0,0.0,0.0);
    std::string  method("pointwise");

    // Create Dirichlet boundary conditions
    std::shared_ptr<DirichletBC> bcInlet = std::make_shared<DirichletBC>((_V), zero,
            reference_to_no_delete_pointer(boundary_mark), 1);
    std::shared_ptr<DirichletBC> bcOutlet = std::make_shared<DirichletBC>((_V), zero,
            reference_to_no_delete_pointer(boundary_mark), 2);
    bcs.push_back(bcInlet);
    bcs.push_back(bcOutlet);

}

    
void VariationalFormsCarotidIso::save_solution()
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

void VariationalFormsCarotidIso::load_solution(std::string str)
{
    File file_backup(str);
    file_backup >> *_u;
}

void VariationalFormsCarotidIso::save_von_misec_stress()
{
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << double(para["pressure_boundary_condition"]);
    std::string idx = std::string("Pres") + ss.str();

    solve(*_a==*_L,*_p);
    ALE::move(*_mesh, *_u);
    File filep(std::string("stress")+idx+std::string(".pvd"));
    filep << *_p;
}

void VariationalFormsCarotidIso::plot_solution()
{
    plot(*_u);
    plot(*_p);
    interactive();

}


