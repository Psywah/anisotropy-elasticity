#include "VariationalFormsAnisoElas.h"
#include "HyperElasticityC.h"
#include "HyperElasticityB.h"
#include "HyperElasticityA.h"
#include "P1VectorElement.h"
using namespace dolfin;




// Sub domain for clamp at left end
class LeftPoint : public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        if( //(std::abs(x[2]) < DOLFIN_EPS) &&
                (std::abs(x[0]+R+Thk_med+Thk_adv) < DOLFIN_EPS) 
          )
        {   
            //std::cout<<"found left point"<<std::endl;
            //std::cout<<x[0]<<"  "<<x[1]<<"  "<<x[2]<<std::endl;
            return true;
        }
        return false;
    }
};

// Sub domain for clamp at right end
class RightPoint : public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        if (//(std::abs(x[2]) < DOLFIN_EPS) &&
                (std::abs(x[0]-R-Thk_med-Thk_adv) < DOLFIN_EPS) 
           )
        {   
            //std::cout<<"found right point"<<std::endl;
            //std::cout<<x[0]<<"  "<<x[1]<<"  "<<x[2]<<std::endl;
            return true;
        }
        return false;
    }
};

// Sub domain for clamp at botton end
class BottomPoint : public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        if(//(std::abs(x[2]) < DOLFIN_EPS) &&
                (std::abs(x[1]+R+Thk_med+Thk_adv) < DOLFIN_EPS) 
          )
        {   
            //std::cout<<"found bottom point"<<std::endl;
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

class MaterialCoef : public Expression
{
    private:
        double adv, med, pla;
    public:

    MaterialCoef(double adv, double med, double pla) : Expression(), adv(adv), med(med), pla(pla){}

        void eval(Array<double>& values, const Array<double>& x) const
        {
            double r = sqrt(x[0]*x[0]+x[1]*x[1]);
            if ( r > R+Thk_med)
                values[0] = adv;
            else if (r > R)
                values[0] = med;
            else
                values[0] = pla;
        }

};


class Orient : public Expression
{
    private:
        double phi_adv, phi_med;
    public:

        Orient(double sign, double Phi_adv, double Phi_med) : Expression(3)
            {phi_adv = sign*DOLFIN_PI*Phi_adv/180.0; phi_med = sign*DOLFIN_PI*Phi_med/180.0;}

        void eval(Array<double>& values, const Array<double>& x) const
        {
            double r = sqrt(x[0]*x[0]+x[1]*x[1]);
            if ( r > R+Thk_med)
            {values[0] = cos(phi_adv)*(-x[1])/r; values[1] = cos(phi_adv)*x[0]/r; values[2] = sin(phi_adv);}
            else if (r > R)
            {values[0] = cos(phi_med)*(-x[1])/r; values[1] = cos(phi_med)*x[0]/r; values[2] = sin(phi_med);}
            else
            {values[0] = 0.0; values[1] = 0.0; values[2] = 0.0;}

        }
};




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


VariationalForms::VariationalForms(
        Mesh& mesh, MeshFunction<std::size_t>& sub_domains_mark,
        MeshFunction<std::size_t>& boundary_mark,
        Parameters& _para): para(_para)
{
    
    _mesh = reference_to_no_delete_pointer(mesh);
    _V_p1 = std::make_shared<P1VectorElement::FunctionSpace>(_mesh);
    _u_p1 = std::make_shared<Function>(_V_p1);
    std::string elas_model = para["elas_model"];
    switch (elas_model[0]) {
        case 'C':
            // FunctionSpace
            info("inital forms of model C");
            _V = std::make_shared<HyperElasticityC::Form_Res::TestSpace>(_mesh);
            _u = std::make_shared<Function>(_V);
            _F = std::make_shared<HyperElasticityC::Form_Res>(_V);
            _obj = std::make_shared<HyperElasticityC::Form_Pi>(_mesh);
            dolfin_assert(_F);
            _J = std::make_shared<HyperElasticityC::Form_Jac>(_V, _V);
            _VMS = std::make_shared<HyperElasticityC::Form_L_vms::TestSpace>(_mesh);
            _p = std::make_shared<Function>(_VMS);
            _a = std::make_shared<HyperElasticityC::Form_Mass_vms>(_VMS,_VMS);
            _L = std::make_shared<HyperElasticityC::Form_L_vms>(_VMS);
            break;
        case 'A':
            info("inital forms of model A");
            _V = std::make_shared<HyperElasticityA::Form_Res::TestSpace>(_mesh);
            _u = std::make_shared<Function>(_V);
            _F = std::make_shared<HyperElasticityA::Form_Res>(_V);
            _obj = std::make_shared<HyperElasticityA::Form_Pi>(_mesh);
            dolfin_assert(_F);
            _J = std::make_shared<HyperElasticityA::Form_Jac>(_V, _V);
            _VMS = std::make_shared<HyperElasticityA::Form_L_vms::TestSpace>(_mesh);
            _p = std::make_shared<Function>(_VMS);
            _a = std::make_shared<HyperElasticityA::Form_Mass_vms>(_VMS,_VMS);
            _L = std::make_shared<HyperElasticityA::Form_L_vms>(_VMS);
            break;
        case 'B':
            info("inital forms of model B");
            _V = std::make_shared<HyperElasticityB::Form_Res::TestSpace>(_mesh);
            _u = std::make_shared<Function>(_V);
            _F = std::make_shared<HyperElasticityB::Form_Res>(_V);
            _obj = std::make_shared<HyperElasticityB::Form_Pi>(_mesh);
            dolfin_assert(_F);
            _J = std::make_shared<HyperElasticityB::Form_Jac>(_V, _V);
            _VMS = std::make_shared<HyperElasticityB::Form_L_vms::TestSpace>(_mesh);
            _p = std::make_shared<Function>(_VMS);
            _a = std::make_shared<HyperElasticityB::Form_Mass_vms>(_VMS,_VMS);
            _L = std::make_shared<HyperElasticityB::Form_L_vms>(_VMS);
            break;
        default :
            dolfin_error("VariationalForms.cpp",
                         "Switch Elasticity Model",
                         "Error option for Elasticity Model");
    }
    info("# of Vertices %d, # of Cells %d, # of Unknowns %d", dolfin::MPI::sum(MPI_COMM_WORLD,num_v), dolfin::MPI::sum(MPI_COMM_WORLD,num_c), _V->dim());
    // Set material parameters
    double C1_adv = (double)(para["C1_adv"]),
           C1_med  = (double)(para["C1_med"]),
           C1_pla  = (double)(para["C1_pla"]),
           Epsilon1_adv  = (double)(para["Epsilon1_adv"]),
           Epsilon1_med  = (double)(para["Epsilon1_med"]),
           Epsilon1_pla  = (double)(para["Epsilon1_pla"]),
           Epsilon2_adv  = (double)(para["Epsilon2_adv"]),
           Epsilon2_med  = (double)(para["Epsilon2_med"]),
           Epsilon2_pla  = (double)(para["Epsilon2_pla"]),
           Alpha3_adv  = (double)(para["Alpha3_adv"]),
           Alpha3_med = (double)(para["Alpha3_med"]),
           Alpha4_adv  = (double)(para["Alpha4_adv"]),
           Alpha4_med  = (double)(para["Alpha4_med"]),
           Alpha5_adv  = (double)(para["Alpha5_adv"]),
           Alpha5_med  = (double)(para["Alpha5_med"]),
           Alpha5_pla  = (double)(para["Alpha5_pla"]),
           K1_adv  = (double)(para["K1_adv"]),
           K2_adv  = (double)(para["K2_adv"]),
           K1_med  = (double)(para["K1_med"]),
           K2_med  = (double)(para["K2_med"]),
           Beta1  = (double)(para["Beta1"]),
           Eta1  = (double)(para["Eta1"]),
           Delta1  = (double)(para["Delta1"]),
           Phi_adv  = (double)(para["Phi_adv"]),
           Phi_med  = (double)(para["Phi_med"]),
           Delta2;
    Delta2 = (Beta1+Eta1*2+Delta1);
    std::shared_ptr<MaterialCoef> coef_c1 = std::make_shared<MaterialCoef>(C1_adv, C1_med, C1_pla);
    std::shared_ptr<MaterialCoef> coef_epsilon1 = std::make_shared<MaterialCoef>(Epsilon1_adv, Epsilon1_med, Epsilon1_pla);
    std::shared_ptr<MaterialCoef> coef_epsilon2 = std::make_shared<MaterialCoef>(Epsilon2_adv, Epsilon2_med, Epsilon2_pla);
    std::shared_ptr<MaterialCoef> coef_alpha3 = std::make_shared<MaterialCoef>(Alpha3_adv, Alpha3_med, 0.0);
    std::shared_ptr<MaterialCoef> coef_alpha4 = std::make_shared<MaterialCoef>(Alpha4_adv, Alpha4_med, 0.0);
    std::shared_ptr<MaterialCoef> coef_alpha5 = std::make_shared<MaterialCoef>(Alpha5_adv, Alpha5_med, Alpha5_pla);
    std::shared_ptr<MaterialCoef> coef_k1 = std::make_shared<MaterialCoef>(K1_adv, K1_med, 0.0);
    std::shared_ptr<MaterialCoef> coef_k2 = std::make_shared<MaterialCoef>(K2_adv, K2_med, 0.0);
    std::shared_ptr<Constant> coef_beta1 = std::make_shared<Constant>(Beta1);
    std::shared_ptr<Constant> coef_eta1 = std::make_shared<Constant>(Eta1);
    std::shared_ptr<Constant> coef_delta1 = std::make_shared<Constant>(Delta1);
    std::shared_ptr<Constant> coef_delta2 = std::make_shared<Constant>(Delta2);
    std::shared_ptr<Orient> vec1 = std::make_shared<Orient>(1.0,Phi_adv,Phi_med), 
                            vec2 = std::make_shared<Orient>(-1.0,Phi_adv,Phi_med);
    // Define source and boundary traction functions
    pressure = para["pressure_boundary_condition"]; t = 1.0;
    final_pressure = pressure;
    std::shared_ptr<PressureNormal> pressure_normal = 
                     std::make_shared<PressureNormal>(mesh,&pressure);
    std::shared_ptr<Constant> B = std::make_shared<Constant>(0.0, 0.0, 0.0);
     
    std::map<std::string, std::shared_ptr<const GenericFunction>> coef_list
        = { {"u",       _u},
            {"vec1",    vec1},
            {"vec2",    vec2},
            {"alpha3",  coef_alpha3},
            {"alpha4",  coef_alpha4},
            {"beta1",   coef_beta1},
            {"eta1",    coef_eta1},
            {"delta1",  coef_delta1},
            {"delta2",  coef_delta2},
            {"c1",      coef_c1},
            {"epsilon1",coef_epsilon1},
            {"epsilon2",coef_epsilon2},
            {"B",       B},
            {"alpha5",  coef_alpha5},
            {"k1",      coef_k1},
            {"k2",      coef_k2},
            {"T",       pressure_normal}};
    // add coef for model A`
    double Alpha1_adv  = (double)(para["Alpha1_adv"]),
           Alpha1_med = (double)(para["Alpha1_med"]),
           Alpha2_adv  = (double)(para["Alpha2_adv"]),
           Alpha2_med  = (double)(para["Alpha2_med"]);
    std::shared_ptr<MaterialCoef> coef_alpha1 = std::make_shared<MaterialCoef>(Alpha1_adv, Alpha1_med, 0.0);
    std::shared_ptr<MaterialCoef> coef_alpha2 = std::make_shared<MaterialCoef>(Alpha2_adv, Alpha2_med, 0.0);
    coef_list["alpha1"] = coef_alpha1;
    coef_list["alpha2"] = coef_alpha2;
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
    std::shared_ptr<LeftPoint> left = std::make_shared<LeftPoint>();
    std::shared_ptr<RightPoint> right = std::make_shared<RightPoint>();
    std::shared_ptr<BottomPoint> bottom = std::make_shared<BottomPoint>();

    // Define Dirichlet boundary functions
    std::shared_ptr<Constant> zero = std::make_shared<Constant>(0.0);
    std::string  method("pointwise");

    // Create Dirichlet boundary conditions
    std::shared_ptr<DirichletBC> bcl = std::make_shared<DirichletBC>((*_V)[1], zero, left, method);
    std::shared_ptr<DirichletBC> bcr = std::make_shared<DirichletBC>((*_V)[1], zero, right, method);
    std::shared_ptr<DirichletBC> bct = std::make_shared<DirichletBC>((*_V)[0], zero, bottom, method);
    std::shared_ptr<DirichletBC> bc_cross = std::make_shared<DirichletBC>((*_V)[2], zero, 
            reference_to_no_delete_pointer(boundary_mark), 2);
    bcs.push_back(bcl);
    bcs.push_back(bcr);
    bcs.push_back(bct);
    bcs.push_back(bc_cross);

}

    
void VariationalForms::save_solution()
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

void VariationalForms::load_solution(std::string str)
{
    File file_backup(str);
    file_backup >> *_u;
}

void VariationalForms::save_von_misec_stress()
{
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << double(para["pressure_boundary_condition"]);
    std::string idx = std::string("Pres") + ss.str();

    solve(*_a==*_L,*_p);
    ALE::move(*_mesh, *_u);
    File filep(std::string("stress")+idx+std::string(".pvd"));
    filep << *_p;
}

void VariationalForms::plot_solution()
{
    plot(*_u);
    plot(*_p);
    interactive();

}


/*
void save_parameters()
{
    Parameters para_newton("newton_solver");
    para_newton = NewtonSolver::default_parameters();

    Parameters para_snes("snes_solver");
    para_snes = PETScSNESSolver::default_parameters();

    list_linear_solver_methods();
    list_krylov_solver_preconditioners();

    Parameters para;

    para.add("nonlinear_solver", "newton"); // "newton" or "snes"
    para.add("pressure_boundary_condition", 24.0); // "newton" or "snes"
    para.add(para_newton);
    para.add(para_snes);

    File file_para("aniso_elas_parameters.xml");
    file_para << para;
}*/
