#include <dolfin.h>
#include "VariationalForms.h"
#include "NonlinearVariationalProblem.h"
#include "NonlinearVariationalSolver.h"

using namespace dolfin;


void list_petsc_snes_methods()
{
  // Pretty-print list of methods
  Table t("PETSc snes method", false);
  for (auto solver : PETScSNESSolver::methods())
    t(solver.first, "Description") = solver.second;
  cout << t.str(true) << endl;
}

void list_petsc_pre_methods()
{
  // Pretty-print list of methods
  Table t("PETSc preconsitioner method", false);
  for (auto pre : PETScPreconditioner::preconditioners())
    t(pre.first, "Description") = pre.second;
  cout << t.str(true) << endl;
}

void list_petsc_ksp_methods()
{
  // Pretty-print list of methods
  Table t("PETSc KrylovSolver method", false);
  for (auto pre : PETScKrylovSolver::methods())
    t(pre.first, "Description") = pre.second;
  cout << t.str(true) << endl;
}


void pseudo_time_steping(double dt, VariationalForms& forms, 
                         NonlinearVariationalSolver& solver);

int main()
{
    // list_linear_solver_methods();
    // list_krylov_solver_preconditioners();
    #ifdef HAS_PETSC
    info("has PETSc");
    list_petsc_snes_methods();
    list_petsc_ksp_methods();
    list_petsc_pre_methods();
    parameters["linear_algebra_backend"] = "PETSc";
    #endif


    Parameters para_material("user_defined_parameters");
    File para_file_material("../aniso_elas_parameters.xml");
    para_file_material >> para_material;
    info(para_material, true);
    Parameters para_nls("nls_parameters");
    File para_file_nls("../nonlinear_solver_parameters.xml");
    para_file_nls >> para_nls;
    
    // Create mesh and define function space
    std::string prefix("../../../mesh/tube-4components-short");
    Mesh mesh(prefix + std::string(".xml"));
    MeshFunction<std::size_t> sub_domains_mark(mesh, 
            prefix+ std::string("-domains-marker.xml") );
    MeshFunction<std::size_t> boundary_mark(mesh, 
            prefix+ std::string("-boundary-marker.xml"));
    //plot(mesh);plot(sub_domains_mark);plot(boundary_mark);
    //interactive();

    /*std::vector<double>& coord = mesh.coordinates();
    for(std::size_t i = 0; i < coord.size(); i++)
        coord[i]*= 1000.;
        */

    Timer t1("Inital Forms"); info("Initial Forms");
    t1.start();
    VariationalForms forms(mesh, sub_domains_mark, boundary_mark, para_material);
    t1.stop();

    //solve(F == 0, u, bcs, J, para);
    Timer t2("Inital Nonlinear Problem"); info("Initial Nonlinear Problem");
    t2.start();
    NonlinearVariationalProblem problem(forms._F, forms._u, forms.bcs, forms._J);
    t2.stop();

    Timer t3("Initial Nonlinear Solver"); info("Initial Nonlinear Solver");
    t3.start();
    NonlinearVariationalSolver solver(problem);
    solver.parameters.update(para_nls);
    t3.stop();

    Timer t4("Solve Nonlinear Problem"); info("Solve Nonlinear Problem");
    double dt = para_nls["dt"];
    t4.start();
    //pseudo_time_steping(dt, forms, solver);
    solver.solve();
    t4.stop();

    forms.save_solution();
    forms.save_von_misec_stress();

    std::set<TimingType> type = {TimingType::wall,TimingType::user, TimingType::system};
    //list_timings(TimingClear::clear,type); 

    forms.plot_solution();


    return 0;

}

void pseudo_time_steping(double dt, VariationalForms& forms, 
                         NonlinearVariationalSolver& solver)
{
    double t = dt;
    while(t < 1.0 + dt -DOLFIN_EPS)
    {
        info("Solve for time %f, current time step %f, End time 1.0", t, dt);
        forms.set_pseudo_time((t>1.0)? 1.0:t);
        solver.solve();
        t += dt;
        //forms.plot_solution();
    }

}


