#ifndef __NONLINEAR_VARIATIONAL_FORMS_ISO_ELAS_H
#define __NONLINEAR_VARIATIONAL_FORMS_ISO_ELAS_H
#include <dolfin.h>


namespace dolfin
{

  /// This is a base class for variational problems which can return the
  /// BilinearForm LinearForm and boundary condition

  class VariationalFormsIsoElas
  {
  public:

    /// Constructor
    VariationalFormsIsoElas(Mesh& mesh, MeshFunction<std::size_t>& sub_domains_mark,
                        MeshFunction<std::size_t>& boundary_mark,
                        Parameters& para);

    /// Destructor
    virtual ~VariationalFormsIsoElas() {}

    void save_solution();

    void load_solution(std::string str="backup_solution.xml");
    void save_von_misec_stress();
    void plot_solution();

    std::shared_ptr<Mesh> _mesh;

    // function space and forms
    std::shared_ptr<FunctionSpace> _V;
    std::shared_ptr<Form> _J;
    std::shared_ptr<Form> _F;

    // boundary conditions
    std::vector<std::shared_ptr<const DirichletBC>> bcs;

    // solution
    std::shared_ptr<Function> _u;

    // von Misec Stress
    std::shared_ptr<FunctionSpace> _VMS;
    std::shared_ptr<Form> _a;
    std::shared_ptr<Form> _L;
    std::shared_ptr<Function> _p;

    // pseudo time stepping
    double t, pressure,final_pressure;
    void set_pseudo_time(double _t){t=_t;pressure=t*final_pressure;}

    Parameters para;

  };

}

#endif


