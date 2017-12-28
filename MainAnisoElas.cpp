//#define HAS_PETSC
#include <dolfin/log/LogManager.h>
#include <dolfin/log/log.h>
#include "NonlinearVariationalSolver.h"
#include "VariationalFormsAnisoElas.h"
#include "NonlinearVariationalProblem.h"
#include <dolfin.h>

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

int main(int argc, char** argv)
{
  #ifdef HAS_PETSC
    parameters["linear_algebra_backend"] = "PETSc";
    PetscInitialize(&argc,&argv,NULL,NULL);
    PetscOptionsInsertFile(PETSC_COMM_WORLD,NULL,"./petsc_options",PETSC_TRUE);
    if(dolfin::MPI::rank(MPI_COMM_WORLD) == 0){
        info("has PETSc");
        list_petsc_snes_methods();
        list_petsc_ksp_methods();
        list_petsc_pre_methods();
    }
    #endif
    LogManager::logger().set_log_active(dolfin::MPI::rank(MPI_COMM_WORLD) == 0);


    Parameters para_material("user_defined_parameters");
    File para_file_material("../parameters/elas_parameters_aniso.xml");
    para_file_material >> para_material;
    info(para_material, true);
    Parameters para_nls("nls_parameters");
    File para_file_nls("./solver_parameters_aniso.xml");
    para_file_nls >> para_nls;
    
    //HDF5File filer(MPI_COMM_WORLD,"/home/gongshihua/work/mesh/tube-2layerL2.h5","r");
    HDF5File filer(MPI_COMM_WORLD,"/home/gongshihua/work/mesh/tube-4components.h5","r");
    
   /* 
    // produce a initial guess
    Parameters para_nls_c("nls_parameters_c");
    File para_file_nls_c("./solver_parameters_aniso_c.xml");
    para_file_nls_c >> para_nls_c;
    Mesh mesh_c;
    int meshID_c  = 1;
    filer.read(mesh_c,std::string("mesh") + std::to_string(meshID_c) ,false);

    // Create mesh functions over the cells and acets
    MeshFunction<std::size_t> sub_domains_mark_c(reference_to_no_delete_pointer(mesh_c), mesh_c.topology().dim());
    MeshFunction<std::size_t> boundary_mark_c(reference_to_no_delete_pointer(mesh_c), mesh_c.topology().dim() - 1);

    filer.read(sub_domains_mark_c,std::string("subdomains_mark") + std::to_string(meshID_c) );
    filer.read(boundary_mark_c,std::string("facet_mark") + std::to_string(meshID_c) );

    Timer t1_c("Inital Forms"); info("Initial Forms");
    t1_c.start();
    VariationalForms forms_c(mesh_c, sub_domains_mark_c, boundary_mark_c, para_material);
    t1_c.stop();

    //solve(F == 0, u, bcs, J, para);
    Timer t2_c("Inital Nonlinear Problem"); info("Initial Nonlinear Problem");
    t2_c.start();
    //NonlinearVariationalProblem problem(forms._F, forms._u, forms.bcs, forms._J, forms._obj);
    NonlinearVariationalProblem problem_c(forms_c._F, forms_c._u, forms_c.bcs, forms_c._J);
    t2_c.stop();

    Timer t3_c("Initial Nonlinear Solver"); info("Initial Nonlinear Solver");
    t3_c.start();
    NonlinearVariationalSolver solver_c(problem_c);
    solver_c.parameters.update(para_nls_c);
    t3_c.stop();

    Timer t4_c("Solve Nonlinear Problem"); info("Solve Nonlinear Problem");
    //solver_c.solve();
    info("Coarse Solver done. \n\n\n\n\n\n\n\n");
    */
    
    
    


    Mesh mesh;
    int meshID  = (int)para_nls["meshID"];
    filer.read(mesh,std::string("mesh") + std::to_string(meshID) ,false);

    // Create mesh functions over the cells and acets
    MeshFunction<std::size_t> sub_domains_mark(reference_to_no_delete_pointer(mesh), mesh.topology().dim());
    MeshFunction<std::size_t> boundary_mark(reference_to_no_delete_pointer(mesh), mesh.topology().dim() - 1);

    filer.read(sub_domains_mark,std::string("subdomains_mark") + std::to_string(meshID) );
    filer.read(boundary_mark,std::string("facet_mark") + std::to_string(meshID) );
    filer.close();


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
    //NonlinearVariationalProblem problem(forms._F, forms._u, forms.bcs, forms._J, forms._obj);
    NonlinearVariationalProblem problem(forms._F, forms._u, forms.bcs, forms._J);
    t2.stop();
    
    /*
    forms_c.load_solution("u1.xml");
    forms_c._u->set_allow_extrapolation(true);
    forms._u->interpolate(*(forms_c._u));
    forms_c.save_solution("u1");
    forms.save_solution(std::string("u") + std::to_string(meshID) );
    //forms_c.backup_solution("u1.xml");
    forms.backup_solution(std::string("u") + std::to_string(meshID) + std::string(".xml"));
    return 0;
    */
    
    
    
    
    //forms.load_solution("u3.xml");
    forms.load_solution(std::string("u") + std::to_string(meshID) + std::string(".xml"));
    //forms.save_solution();
    //return 0;

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
    list_timings(TimingClear::clear,type); 

    //forms.plot_solution();


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


