#include <dolfin.h>
#include "VariationalForms.h"
#include "NonlinearVariationalProblem.h"
#include "NonlinearVariationalSolver.h"

using namespace dolfin;


void list_petsc_snes_methods()
{
  // Get methods
  //std::map<std::string, std::string> methods = PETScSNESSolver::_methods;

  // Pretty-print list of methods
  Table t("PETSc snes method", false);
  for (auto solver : PETScSNESSolver::methods())
    t(solver.first, "Description") = solver.second;
  cout << t.str(true) << endl;
}

int main()
{
    // list_linear_solver_methods();
    // list_krylov_solver_preconditioners();
    #ifdef HAS_PETSC
    info("has PETSc");
    // list_petsc_snes_methods();
    parameters["linear_algebra_backend"] = "PETSc";
    #endif


    Parameters para("user_defined_parameters");
    File para_file("../aniso_elas_parameters.xml");
    para_file >> para;
    info(para("problem_parameters"));
    
    // Create mesh and define function space
    std::string prefix("../../../mesh/tube-4components");
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

    
    VariationalForms forms(mesh, sub_domains_mark, boundary_mark, para("problem_parameters"));

    //solve(F == 0, u, bcs, J, para);
    NonlinearVariationalProblem problem(forms._F, forms._u, forms.bcs, forms._J);
    NonlinearVariationalSolver solver(problem);
    solver.parameters.update(para("nonlinear_solver_parameters"));
    solver.solve();

    forms.save_solution();
    forms.save_von_misec_stress();
    forms.plot_solution();


    return 0;

}


