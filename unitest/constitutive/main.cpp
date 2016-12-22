#include <dolfin.h>
#include "HyperElasticity.h"
#include "Von_Misec_Stress.h"

using namespace dolfin;

#define R 0.01
#define Thk_med 0.00132
#define Thk_adv 0.00096
#define Depth 0.005
#define C1_adv 6.8e3
#define C1_med 11.8e3
#define Epsilon1_adv 52.5e3
#define Epsilon1_med 84.2e3
#define Epsilon2_adv 10.
#define Epsilon2_med 5.94
#define Alpha3_adv 1005.2e3
#define Alpha3_med 49999.1e3
#define Alpha4_adv 6.3
#define Alpha4_med 4.1
#define Beta1 80.e3
#define Eta1 250.e3
#define Delta1 2000.e3
#define Delta2 (Beta1+Eta1*2+Delta1)
#define Phi_adv 57.1
#define Phi_med 28.8
#define Pressure 1500

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

        PressureNormal(const Mesh& mesh) : Expression(3), mesh(mesh) {}

        void eval(Array<double>& values, const Array<double>& x,
                const ufc::cell& ufc_cell) const
        {
            dolfin_assert(ufc_cell.local_facet >= 0);

            Cell cell(mesh, ufc_cell.index);
            Point n = cell.normal(ufc_cell.local_facet);

            const double pressure = Pressure;
            values[0] = -pressure*n[0];
            values[1] = -pressure*n[1];
            values[2] = -pressure*n[2];
        }

    private:

        const Mesh& mesh;

};

class MCoef_c1 : public Expression
{
    public:

        MCoef_c1() : Expression() {}

        void eval(Array<double>& values, const Array<double>& x) const
        {
            double r = sqrt(x[0]*x[0]+x[1]*x[1]);
            if ( r > R+Thk_med)
                values[0] = C1_adv;
            else if (r > R)
                values[0] = C1_med;
            else
                values[0] = C1_med;
        }

};

class MCoef_alpha3 : public Expression
{
    public:

        MCoef_alpha3() : Expression() {}

        void eval(Array<double>& values, const Array<double>& x) const
        {
            double r = sqrt(x[0]*x[0]+x[1]*x[1]);
            if ( r > R+Thk_med)
                values[0] = Alpha3_adv;
            else if (r > R)
                values[0] = Alpha3_med;
            else
                values[0] = 0.0;
        }

};

class MCoef_alpha4 : public Expression
{
    public:

        MCoef_alpha4() : Expression() {}

        void eval(Array<double>& values, const Array<double>& x) const
        {
            double r = sqrt(x[0]*x[0]+x[1]*x[1]);
            if ( r > R+Thk_med)
                values[0] = Alpha4_adv;
            else if (r > R)
                values[0] = Alpha4_med;
            else
                values[0] = 0.0;
        }

};

class MCoef_epsilon1 : public Expression
{
    public:

        MCoef_epsilon1() : Expression() {}

        void eval(Array<double>& values, const Array<double>& x) const
        {
            double r = sqrt(x[0]*x[0]+x[1]*x[1]);
            if ( r > R+Thk_med)
                values[0] = Epsilon1_adv;
            else {
                if (r > R)
                    values[0] = Epsilon1_med;
                else
                    values[0] = Epsilon1_med;
            }
        }

};

class MCoef_epsilon2 : public Expression
{
    public:

        MCoef_epsilon2() : Expression() {}

        void eval(Array<double>& values, const Array<double>& x) const
        {
            double r = sqrt(x[0]*x[0]+x[1]*x[1]);
            if ( r > R+Thk_med)
                values[0] = Epsilon2_adv;
            else if (r > R)
                values[0] = Epsilon2_med;
            else
                values[0] = Epsilon2_med;
        }

};


class Orient : public Expression
{
    private:
        double phi_adv, phi_med;
    public:

        Orient(double sign) : Expression(3) {phi_adv = sign*DOLFIN_PI*Phi_adv/180.0; phi_med = sign*DOLFIN_PI*Phi_med/180.0;}

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


int main()
{
    // Create mesh and define function space
    Mesh mesh("../../mesh/tube-2layer.xml");
    MeshFunction<std::size_t> sub_domains_mark(mesh, 
            "../../mesh/tube-2layer-domains-marker.xml" );
    MeshFunction<std::size_t> boundary_mark(mesh, 
            "../../mesh/tube-2layer-boundary-marker.xml");

    HyperElasticity::FunctionSpace V(mesh);


    // Define Dirichlet boundaries
    LeftPoint left;
    RightPoint right;
    BottomPoint bottom;

    // Define Dirichlet boundary functions
    Constant zero(0.0);
    std::string  method("pointwise");

    // Create Dirichlet boundary conditions
    DirichletBC bcl(*(V[1]), zero, left, method);
    DirichletBC bcr(*(V[1]), zero, right, method);
    DirichletBC bct(*(V[0]), zero, bottom, method );
    DirichletBC bc_cross(*(V[2]), zero, boundary_mark, 2);
    std::vector<const DirichletBC*> bcs = {{&bcl, &bcr, &bct, &bc_cross}};

    // Define source and boundary traction functions
    PressureNormal pressure_normal(mesh);
    Constant B(0.0, 0.0, 0.0);
    Constant T(0.1,  0.0, 0.0);


    // Set material parameters
    MCoef_c1 coef_c1;
    MCoef_epsilon1 coef_epsilon1;
    MCoef_epsilon2 coef_epsilon2;
    MCoef_alpha3 coef_alpha3;
    MCoef_alpha3 coef_alpha4;
    Constant coef_beta1(Beta1);
    Constant coef_eta1(Eta1);
    Constant coef_delta1(Delta1);
    Constant coef_delta2(Delta2);
    Orient vec1(1.0),vec2(-1.0);

    /*
    HyperElasticity::CoefficientSpace_vec1 V_vec(mesh);
    Function fun_vec1(V_vec);
    Function fun_vec2(V_vec);
    fun_vec1.interpolate(vec1);
    fun_vec2.interpolate(vec2);
    // Save solution in VTK format
    File fileOrient1("1orient.pvd");
    File fileOrient2("2orient.pvd");
    fileOrient1 << fun_vec1;
    fileOrient2 << fun_vec2;

    HyperElasticity::CoefficientSpace_c1 DG(mesh);
    Function fun_c1(DG);
    fun_c1 = coef_c1;
    // Save solution in VTK format
    File fileC1("C1.pvd");
    fileC1 << fun_c1;
    //plot(fun_vec1);
    //plot(fun_vec2);
    //interactive();
    */

    // Define solution function
    Function u(V);
    u = Constant(0.0,0.0,0.0);


    // Create (linear) form defining (nonlinear) variational problem
    HyperElasticity::ResidualForm F(V);
    F.vec1 = vec1; F.vec2 = vec2;
    F.alpha3 = coef_alpha3; F.alpha4 = coef_alpha4;
    F.beta1 = coef_beta1; F.eta1 = coef_eta1; F.delta1=coef_delta1; F.delta2=coef_delta2;
    F.c1 = coef_c1; F.epsilon1 = coef_epsilon1; F.epsilon2 = coef_epsilon2; F.u = u;
    F.B = B; F.T = pressure_normal;
    F.ds = boundary_mark;
    F.dx = sub_domains_mark;

    Vector b;
    assemble(b, F);
    Function w( reference_to_no_delete_pointer(V),reference_to_no_delete_pointer(b));

    //for(std::size_t i = 0; i  < b.size(); ++i)
    //    std:cout<< b[i]<<std:endl;
    //std::cout << b.max()<<std::endl;
    //std::cout << b.min()<<std::endl;
    std::cout << b.size()<<std::endl;
    //info(b,true);

    
    // Create jacobian dF = F' (for use in nonlinear solver).
    HyperElasticity::JacobianForm J(V, V);
    J.vec1 = vec1; J.vec2 = vec2;
    J.alpha3 = coef_alpha3; J.alpha4 = coef_alpha4;
    J.beta1 = coef_beta1; J.eta1 = coef_eta1; J.delta1=coef_delta1; J.delta2=coef_delta2;
    J.c1 = coef_c1; J.epsilon1 = coef_epsilon1; J.epsilon2 = coef_epsilon2; J.u = u;
    J.dx = sub_domains_mark;
    //Matrix A;
    //assemble(A, J);
    //info(A,true);
    

    // Solve nonlinear variational problem F(u; v) = 0
    solve(F == 0, u, bcs, J);
    

    std::string idx("P1e3");
    // Save solution in VTK format
    File file(std::string("displacement")+idx+std::string(".pvd"));
    file << u;
    Von_Misec_Stress::FunctionSpace VMS(mesh);
    Von_Misec_Stress::BilinearForm a(VMS,VMS);
    Von_Misec_Stress::LinearForm L(VMS);
    L.vec1 = vec1; L.vec2 = vec2;
    L.alpha3 = coef_alpha3; L.alpha4 = coef_alpha4;
    L.beta1 = coef_beta1; L.eta1 = coef_eta1; L.delta1=coef_delta1; L.delta2=coef_delta2;
    L.c1 = coef_c1; L.epsilon1 = coef_epsilon1; L.epsilon2 = coef_epsilon2; L.u = u;
    L.dx = sub_domains_mark;
    Function p(VMS);
    solve(a==L,p);
    File filep(std::string("stress")+idx+std::string(".pvd"));
    filep << p;

    // Plot solution
    plot(p);

    interactive();
    
    


    return 0;
}
