// Copyright (C) 2008-2011 Anders Logg and Garth N. Wells
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Marie E. Rognes 2011
// Modified by Corrado Maurini 2013
//
// First added:  2011-01-14 (2008-12-26 as VariationalProblem.cpp)
// Last changed: 2013-11-21

//#define HAS_PETSC
#include <petscsys.h>
#include <dolfin/common/NoDeleter.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/function/Function.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/fem/assemble.h>
#include <dolfin/fem/Form.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/PETScVector.h>
#include "NonlinearVariationalProblem.h"
#include "NonlinearVariationalSolver.h"
#include <dolfin/common/Timer.h>

//#define SAVE_RESIDULE_VEC

#define NONLINEAR_ELIMINATION
#if defined(SAVE_RESIDULE_VEC) || defined(NONLINEAR_ELIMINATION)
#include <dolfin/io/dolfin_io.h>
#endif

using namespace dolfin;

//-----------------------------------------------------------------------------
NonlinearVariationalSolver::
NonlinearVariationalSolver(NonlinearVariationalProblem& problem)
  : _problem(reference_to_no_delete_pointer(problem))
{
  // Set parameters
  parameters = default_parameters();
  File para("nls_parameters.xml");
  para <<parameters;
}
//-----------------------------------------------------------------------------
NonlinearVariationalSolver::
NonlinearVariationalSolver(std::shared_ptr<NonlinearVariationalProblem> problem)
  : _problem(problem)
{
  // Set parameters
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
std::pair<std::size_t, bool>
NonlinearVariationalSolver::solve(const Function& lb,
                                  const Function& ub)
{
  return solve(lb.vector(), ub.vector());
}
//-----------------------------------------------------------------------------
std::pair<std::size_t, bool>
NonlinearVariationalSolver::solve(std::shared_ptr<const Function> lb,
                                  std::shared_ptr<const Function> ub)
{
  return solve(*lb,*ub);
}
//-----------------------------------------------------------------------------
std::pair<std::size_t, bool>
NonlinearVariationalSolver::solve(const GenericVector& lb,
                                  const GenericVector& ub)
{
  return solve(reference_to_no_delete_pointer(lb),
               reference_to_no_delete_pointer(ub));
}
//-----------------------------------------------------------------------------
std::pair<std::size_t, bool>
NonlinearVariationalSolver::solve(std::shared_ptr<const GenericVector> lb,
                                  std::shared_ptr<const GenericVector> ub)
{
  // Set bounds and solve
  this->_problem->set_bounds(lb,ub);
  return solve();
}
//-----------------------------------------------------------------------------
std::pair<std::size_t, bool> NonlinearVariationalSolver::solve()
{
  begin("Solving nonlinear variational problem.");

  // Check that the Jacobian has been defined
  dolfin_assert(_problem);
  if (!_problem->has_jacobian())
  {
    dolfin_error("NonlinearVariationalSolver.cpp",
                 "solve nonlinear variational problem",
                 "The Jacobian form has not been defined");
  }
  // Check that the dolfin is configured with petsc is bounds are set
#ifndef HAS_PETSC
  if (_problem->has_lower_bound() || _problem->has_upper_bound())
  {
    dolfin_error("NonlinearVariationalSolver.cpp",
                 "solve nonlinear variational problem",
                 "Needs PETSc to solve bound constrained problems");
  }
#endif
  // Get problem data
  dolfin_assert(_problem);
  std::shared_ptr<Function> u(_problem->solution());

  // Create nonlinear problem
  if (!nonlinear_problem)
  {
    nonlinear_problem = std::make_shared<NonlinearDiscreteProblem>(_problem,
                                                    reference_to_no_delete_pointer(*this));
  }
  // Create nonlinear coarse problem
  if (!nonlinear_coarse_problem)
  {
    nonlinear_coarse_problem = std::make_shared<NonlinearCoarseDiscreteProblem>(_problem,
                                                    reference_to_no_delete_pointer(*this));
  }

  std::pair<std::size_t, bool> ret;
  if (std::string(parameters["nonlinear_solver"]) == "newton")
  {
    if (_problem->has_lower_bound() && _problem->has_upper_bound())
    {
      dolfin_error("NonlinearVariationalSolver.cpp",
                   "solve nonlinear variational problem",
                   "Set the \"nonlinear_solver\" parameter to \"snes\" or remove bounds");
    }
    // Create Newton solver and set parameters
    if (!newton_solver)
      newton_solver = std::shared_ptr<NewtonSolver>(new NewtonSolver());

    // Pass parameters to Newton solver
    newton_solver->parameters.update(parameters("newton_solver"));

    // Solve nonlinear problem using Newton's method
    dolfin_assert(u->vector());
    dolfin_assert(nonlinear_problem);
    ret = newton_solver->solve(*nonlinear_problem, *u->vector());
  }
  #ifdef HAS_PETSC
  else if (std::string(parameters["nonlinear_solver"]) == "snes")
  {
    // Create SNES solver and set parameters
    if (!snes_solver)
    {
      // Create Newton solver and set parameters
      snes_solver = std::make_shared<PETScSNESSolver>();
    }
    snes_solver->parameters.update(parameters("snes_solver"));

    // Solve nonlinear problem using PETSc's SNES
    dolfin_assert(u->vector());
    dolfin_assert(nonlinear_problem);
    if (_problem->has_lower_bound() && _problem->has_upper_bound())
    {
      ret = snes_solver->solve(*nonlinear_problem, *u->vector(),
                               *_problem->lower_bound(),
                               *_problem->upper_bound());
    }
    else
    {
      #ifdef SAVE_RESIDULE_VEC
       info("begin snes solve, save each residual");
       File u_hist("u_hist.pvd");
       File du_hist("du_hist.pvd");
       File r_hist("r_hist.pvd");
       Function u0(u->function_space()), du(u->function_space()), res(u->function_space());
       std::size_t iter = 0 ;
       bool converged = false;
       while(iter< 30 && !converged)
       {
           nonlinear_problem->F(*res.vector(), *u->vector());
           *(du.vector()) = *(u->vector());
           ret = snes_solver->solve(*nonlinear_problem, *u->vector());
           *(du.vector()) *= -1.0;
           *(du.vector()) += *(u->vector());
           ++iter;
           converged = std::get<1>(ret);
           info("save residule of iterations (%d).\n", iter);
           u_hist << *u ;
           du_hist << du ;
           r_hist << res ;
       }
      #elif defined NONLINEAR_ELIMINATION
       if(double(parameters["NL_res_r"]) >1)
       {
           info("begin snes solve");
           ret = snes_solver->solve(*nonlinear_problem, *u->vector());
       }else
       {
           std::string prefixdir("./results/");
           std::shared_ptr<File> u_hist = std::make_shared<File>(prefixdir+ "u/u_hist.pvd");
           std::shared_ptr<File> r_hist = std::make_shared<File>(prefixdir+ "r/r_hist.pvd");
           std::shared_ptr<File> rc_hist = std::make_shared<File>(prefixdir+ "rc/rc_hist.pvd");
           std::shared_ptr<File> bad_dof = std::make_shared<File>(prefixdir+ "dof/bad_dof_hist.pvd");
           std::shared_ptr<File> bad_cell = std::make_shared<File>(prefixdir+ "celldof/bad_cell_hist.pvd");
           std::shared_ptr<File> bad_cell_ov = std::make_shared<File>(prefixdir+ "ovdof/bad_ov_hist.pvd");
           info("begin snes solve with nonlinear elimination");
           Parameters para_coarse = parameters("snes_solver");
           para_coarse["absolute_tolerance"] = (double)(parameters["NL_atol"]);
           para_coarse["relative_tolerance"] = (double)(parameters["NL_rtol"]);
           para_coarse["line_search"] = std::string(parameters["NL_line_search"]);
           para_coarse["report"] = (bool)(parameters["NL_report"]);
           para_coarse["maximum_iterations"] = (int)(parameters["NL_maximum_iterations"]);
           double rho0 = parameters["NL_res_r"];
           double rho1 = parameters["NL_infty_r"];
           double rho2 = parameters["NL_size_r"];
           int overlap = parameters["NL_overlap"];

           Function res(u->function_space());
           double norm_res0,norm_res;
           double obj0,obj;
           if(_problem->has_object())
               nonlinear_problem->Object(*res.vector(),obj0, *u->vector());
           nonlinear_problem->F(*res.vector(), *u->vector());
           norm_res0 =res.vector()->norm("l2");

           bool converged = false;
           bool accept = false;
           std::size_t iter = 0, max_iter = parameters("snes_solver")["maximum_iterations"] ;
           parameters("snes_solver")["maximum_iterations"] = 1;
           while(!converged && iter< max_iter)
           {
               if(iter!=0 && iter%30==0)
               {
                   std::string str(std::string("backup_solution")+(".xml"));
                   File backup_file(str);
                   backup_file << *u;
               }
               snes_solver->parameters.update(parameters("snes_solver"));
               ret = snes_solver->solve(*nonlinear_problem, *u->vector());

               
               *u_hist << *u;
               ++iter;
               converged = std::get<1>(ret);
               if(_problem->has_object())
                   nonlinear_problem->Object(*res.vector(),obj, *u->vector());
               nonlinear_problem->F(*res.vector(), *u->vector());
               *r_hist << res ;
               norm_res =res.vector()->norm("l2");
               info(ANSI_COLOR_MAGENTA "%d iterations " ANSI_COLOR_RESET  "of global nonlinear iteration:", iter);
               info("       last residual: %5.2e, vs current " ANSI_COLOR_MAGENTA "residual: %5.3e " ANSI_COLOR_RESET,
                       norm_res0, norm_res);
               if(_problem->has_object())
                   info("       last object: %5.2e, vs current " ANSI_COLOR_MAGENTA "object: %5.3e " ANSI_COLOR_RESET,
                           obj0, obj);
               if(converged)
                   break;
               begin("Check Nonlinearity.");

               // find all good dof
               double infty = res.vector()->norm("linf");
               info("finding bad dofs [ > %f*%.2f (infty norm * rho1)]", infty, rho1);
               PETScMatrix J;
               nonlinear_problem->J(J,*u->vector());
               nonlinear_coarse_problem->construct_coarse_IS(*res.vector(), infty*rho1, J, overlap);

               if(norm_res < norm_res0*rho0)
               {
                   info("Reduction of residual " ANSI_COLOR_GREEN "%.2f " ANSI_COLOR_RESET "< rho0(%.2f), No Need for nonlinear elimination",
                           norm_res/norm_res0, rho0);
                   norm_res0 = norm_res;
                   if(_problem->has_object())
                       obj0 = obj;
                   nonlinear_coarse_problem->clear_dofs_values();
                   nonlinear_coarse_problem->save_coarse(bad_dof,bad_cell,bad_cell_ov);
                   *rc_hist << res ;
                   end();
                   continue;
               }
               info("Reduction of residual " ANSI_COLOR_RED "%.2f " ANSI_COLOR_RESET "> rho0(%.2f), Need for nonlinear elimination",
                       norm_res/norm_res0, rho0);

               nonlinear_coarse_problem->save_coarse(bad_dof,bad_cell,bad_cell_ov);

               std::size_t size_total = res.vector()->size();
               std::size_t size_good = nonlinear_coarse_problem->size_dirichlet_dofs_global();

               if( size_total - size_good > rho2*size_total || size_total -size_good ==0)
               {
                   info("Size of coarse problem %d > %d*%.2f (total_size * rho2), No Need for nonlinear elimination",
                           size_total-size_good, size_total, rho2);
                   norm_res0 = norm_res;
                   if(_problem->has_object())
                       obj0 = obj;
                   *rc_hist << res ;
                   end();
                   continue;
               }
               info("Size of coarse problem %d < %d * %.2f (total_size * rho_2), Need for nonlinear elimination",
                       size_total-size_good, size_total, rho2);
               PETScVector x_copy(u->vector()->down_cast<PETScVector>());

               snes_solver->parameters.update(para_coarse);
               ret = snes_solver->solve(*nonlinear_coarse_problem, *u->vector());
               if(1)// restricted
               {
                   nonlinear_coarse_problem->restricted_update(*u->vector(), x_copy);
               }

               double obj_c;
               if(_problem->has_object())
                   nonlinear_problem->Object(*res.vector(),obj_c, *u->vector());
               nonlinear_problem->F(*res.vector(), *u->vector());
               *rc_hist << res ;
               double norm_res_c = res.vector()->norm("l2");

               if(_problem->has_object())
               {
                   norm_res0 = norm_res_c;
                   obj0 = obj_c;
                   info("elimination result: residual %f --> %f , object %f --> %f", norm_res, norm_res_c, obj, obj_c);
               }
               else{
                   if(iter <10 || !accept || norm_res_c<norm_res)
                   {
                       norm_res0 = norm_res_c;
                       accept=true;
                       if(iter < 10)
                           info("elimination at the first 10 steps %f < %f, " ANSI_COLOR_GREEN "accept " ANSI_COLOR_RESET "nonlinear elimination",
                                   norm_res_c,norm_res);
                       else if(!accept)
                           info("elimination since last not accept %f < %f, " ANSI_COLOR_GREEN "accept " ANSI_COLOR_RESET "nonlinear elimination",
                                   norm_res_c,norm_res);
                       else
                           info("elimination Residual norm with correction %f < %f, " ANSI_COLOR_GREEN "accept " ANSI_COLOR_RESET "nonlinear elimination",
                                   norm_res_c,norm_res);
                   }
                   else
                   {
                       norm_res0 =norm_res;
                       *u->vector() = x_copy;
                       accept = false;
                       info("Residual norm with correction %f > %f, " ANSI_COLOR_RED "not accept " ANSI_COLOR_RESET "nonlinear elimination",
                               norm_res_c,norm_res);
                   }
               }
               info("Solved Coarse Problem");
               end();

           }
       }
#else
       info("begin snes solve");
       ret = snes_solver->solve(*nonlinear_problem, *u->vector());
#endif
    }
  }
#endif
  else
  {
      dolfin_error("NonlinearVariationalSolver.cpp",
              "solve nonlinear variational problem",
              "Unknown nonlinear solver type");
  }

  end();
  return ret;
}

//-----------------------------------------------------------------------------
// Implementation of NonlinearDiscreteProblem
//-----------------------------------------------------------------------------
NonlinearVariationalSolver::NonlinearDiscreteProblem::
NonlinearDiscreteProblem(std::shared_ptr<NonlinearVariationalProblem> problem,
        std::shared_ptr<NonlinearVariationalSolver> solver)
: _problem(problem), _solver(solver)
{
    // Do nothing
}
//-----------------------------------------------------------------------------
NonlinearVariationalSolver::NonlinearDiscreteProblem::~NonlinearDiscreteProblem()
{
    // Do nothing
}
//-----------------------------------------------------------------------------
void NonlinearVariationalSolver::
NonlinearDiscreteProblem::F(GenericVector& b, const GenericVector& x)
{
    Timer t("Assemble Residual");
    // Get problem data
    dolfin_assert(_problem);
    std::shared_ptr<const Form> F(_problem->residual_form());
    std::vector<std::shared_ptr<const DirichletBC>> bcs(_problem->bcs());

    // Assemble right-hand side
    dolfin_assert(F);
    Assembler assembler;
    assembler.assemble(b, *F);

    // Apply boundary conditions
    for (std::size_t i = 0; i < bcs.size(); i++)
    {
        dolfin_assert(bcs[i]);
        bcs[i]->apply(b, x);
    }

    // Print vector
    dolfin_assert(_solver);
    const bool print_rhs = _solver->parameters["print_rhs"];
    if (print_rhs)
        info(b, true);
}
//-----------------------------------------------------------------------------
void NonlinearVariationalSolver::
NonlinearDiscreteProblem::Object(GenericVector& b,double & val,  const GenericVector& x)
{
    if(_problem->has_object())
    {
        Timer t("Assemble Object");
        // Get problem data
        dolfin_assert(_problem);
        std::shared_ptr<const Form> obj(_problem->object_form());

        // Assemble right-hand side
        dolfin_assert(obj);
        val = dolfin::assemble(*obj);
    }
    else{
        Timer t("Assemble Residual");
        // Get problem data
        dolfin_assert(_problem);
        std::shared_ptr<const Form> F(_problem->residual_form());
        std::vector<std::shared_ptr<const DirichletBC>> bcs(_problem->bcs());

        // Assemble right-hand side
        dolfin_assert(F);
        Assembler assembler;
        assembler.assemble(b, *F);

        // Apply boundary conditions
        for (std::size_t i = 0; i < bcs.size(); i++)
        {
            dolfin_assert(bcs[i]);
            bcs[i]->apply(b, x);
        }
        b.apply("insert");

        // Print vector
        dolfin_assert(_solver);
        const bool print_rhs = _solver->parameters["print_rhs"];
        if (print_rhs)
            info(b, true);
        val = b.norm("l2");
    }
}
//-----------------------------------------------------------------------------
    void
NonlinearVariationalSolver::NonlinearDiscreteProblem::J(GenericMatrix& A,
        const GenericVector& x)
{
    Timer t("Assemble Jacobian");
    // Get problem data
    dolfin_assert(_problem);
    std::shared_ptr<const Form> J(_problem->jacobian_form());
    std::vector<std::shared_ptr<const DirichletBC>> bcs(_problem->bcs());

    // Assemble left-hand side
    dolfin_assert(J);
    Assembler assembler;
    assembler.assemble(A, *J);

    // Apply boundary conditions
    for (std::size_t i = 0; i < bcs.size(); i++)
    {
        dolfin_assert(bcs[i]);
        bcs[i]->apply(A);
    }

    // Print matrix
    dolfin_assert(_solver);
    const bool print_matrix = _solver->parameters["print_matrix"];
    if (print_matrix)
        info(A, true);
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Implementation of NonlinearCoarseDiscreteProblem
//-----------------------------------------------------------------------------
NonlinearVariationalSolver::NonlinearCoarseDiscreteProblem::
NonlinearCoarseDiscreteProblem(std::shared_ptr<NonlinearVariationalProblem> problem,
        std::shared_ptr<NonlinearVariationalSolver> solver)
: _problem(problem), _solver(solver)
{
    // Do nothing
}
//-----------------------------------------------------------------------------
NonlinearVariationalSolver::NonlinearCoarseDiscreteProblem::~NonlinearCoarseDiscreteProblem()
{
    // Do nothing
}
//-----------------------------------------------------------------------------
void NonlinearVariationalSolver::
NonlinearCoarseDiscreteProblem::F(GenericVector& b, const GenericVector& x)
{
    Timer t("Assemble Residual");
    // Get problem data
    dolfin_assert(_problem);
    std::shared_ptr<const Form> F(_problem->residual_form());
    std::vector<std::shared_ptr<const DirichletBC>> bcs(_problem->bcs());

    // Assemble right-hand side
    dolfin_assert(F);
    Assembler assembler;
    assembler.assemble(b, *F);

    // Apply boundary conditions
    for (std::size_t i = 0; i < bcs.size(); i++)
    {
        dolfin_assert(bcs[i]);
        bcs[i]->apply(b, x);
    }
    b.set(values.data(), values.size(), dirichlet_dofs.data());
    b.apply("insert");

    // Print vector
    dolfin_assert(_solver);
    const bool print_rhs = _solver->parameters["print_rhs"];
    if (print_rhs)
        info(b, true);
}
//-----------------------------------------------------------------------------
void NonlinearVariationalSolver::
NonlinearCoarseDiscreteProblem::Object(GenericVector& b,double & val,  const GenericVector& x)
{
    if(_problem->has_object())
    {
        Timer t("Assemble Object");
        // Get problem data
        dolfin_assert(_problem);
        std::shared_ptr<const Form> obj(_problem->object_form());

        // Assemble right-hand side
        dolfin_assert(obj);
        val = dolfin::assemble(*obj);
    }
    else{
        Timer t("Assemble Residual");
        // Get problem data
        dolfin_assert(_problem);
        std::shared_ptr<const Form> F(_problem->residual_form());
        std::vector<std::shared_ptr<const DirichletBC>> bcs(_problem->bcs());

        // Assemble right-hand side
        dolfin_assert(F);
        Assembler assembler;
        assembler.assemble(b, *F);

        // Apply boundary conditions
        for (std::size_t i = 0; i < bcs.size(); i++)
        {
            dolfin_assert(bcs[i]);
            bcs[i]->apply(b, x);
        }
        b.set(values.data(), values.size(), dirichlet_dofs.data());
        b.apply("insert");

        // Print vector
        dolfin_assert(_solver);
        const bool print_rhs = _solver->parameters["print_rhs"];
        if (print_rhs)
            info(b, true);
        val = b.norm("l2");
    }
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
    void
NonlinearVariationalSolver::NonlinearCoarseDiscreteProblem::J(GenericMatrix& A,
        const GenericVector& x)
{
    Timer t("Assemble Jacobian");
    // Get problem data
    dolfin_assert(_problem);
    std::shared_ptr<const Form> J(_problem->jacobian_form());
    std::vector<std::shared_ptr<const DirichletBC>> bcs(_problem->bcs());

    // Assemble left-hand side
    dolfin_assert(J);
    Assembler assembler;
    assembler.assemble(A, *J);

    // Apply boundary conditions
    for (std::size_t i = 0; i < bcs.size(); i++)
    {
        dolfin_assert(bcs[i]);
        bcs[i]->apply(A);
    }
    A.ident(dirichlet_dofs.size(),dirichlet_dofs.data());
    A.apply("insert");

    // Print matrix
    dolfin_assert(_solver);
    const bool print_matrix = _solver->parameters["print_matrix"];
    if (print_matrix)
        info(A, true);
}
//-----------------------------------------------------------------------------

void NonlinearVariationalSolver::NonlinearCoarseDiscreteProblem::add_dirichlet_dof(dolfin::la_index dof, double value)
{
    dirichlet_dofs.push_back(dof); values.push_back(value);
}


void NonlinearVariationalSolver::NonlinearCoarseDiscreteProblem::clear_dofs_values()
{
    dirichlet_dofs.clear(); 
    values.clear(); 
    bad_dofs.clear(); 
    dirichlet_dofs0.clear();
    bad_cell_dofs.clear();
    bad_cell_ov_dofs.clear();
}


void NonlinearVariationalSolver::NonlinearCoarseDiscreteProblem::set_all_dofs()
{
    clear_dofs_values();
    std::pair<std::size_t, std::size_t> range = _problem->solution()->vector()->local_range();
    for(dolfin::la_index dof = range.first; dof<range.second; dof++)
    {
        //double value;
        //u.get(&value, 1, &dof);
        //add_good_dof(dof, value);
        add_dirichlet_dof(dof, 0.0);
    }
    dirichlet_dofs0 = dirichlet_dofs;
}

std::size_t NonlinearVariationalSolver::NonlinearCoarseDiscreteProblem::size_dirichlet_dofs_global()
{
    std::size_t local_size=dirichlet_dofs.size();
    return dolfin::MPI::sum<std::size_t>(PETSC_COMM_WORLD, local_size);
}

std::size_t NonlinearVariationalSolver::NonlinearCoarseDiscreteProblem::size_bad_dofs()
{
    std::size_t local_size=bad_dofs.size();
    return dolfin::MPI::sum<std::size_t>(PETSC_COMM_WORLD,local_size);
}

void NonlinearVariationalSolver::NonlinearCoarseDiscreteProblem::restricted_update(GenericVector& u, const GenericVector& u_pre)
{
    //double* value= new double[bad_dofs.size()];
    //u.get(value, bad_dofs.size(), bad_dofs.data());
    //u = u_pre;
    //u.set(value, bad_dofs.size(), bad_dofs.data());
    //info("last updated value%5.2e, #bad_dofs: %d", value[bad_dofs.size()-1], bad_dofs.size());
    //u.apply("insert");
    //delete[] value;

    //double* value= new double[bad_cell_dofs.size()];
    //u.get(value, bad_cell_dofs.size(), bad_cell_dofs.data());
    //u = u_pre;
    //u.set(value, bad_cell_dofs.size(), bad_cell_dofs.data());
    //info("last updated value%5.2e, #bad_dofs: %d", value[bad_dofs.size()-1], bad_dofs.size());
    //u.apply("insert");
    //delete[] value;
    //
    //u += u_pre;
    //u *=.5;
}

void NonlinearVariationalSolver::NonlinearCoarseDiscreteProblem::save_coarse(std::shared_ptr<File> file_badnode, std::shared_ptr<File> file_badcell, std::shared_ptr<File> file_ov)
{
    Function u(_problem->trial_space());
    std::vector<double> value;
    value.resize(bad_dofs.size(),1);
    u.vector()->zero();
    u.vector()->set_local(value.data(), bad_dofs.size(),bad_dofs.data());
    u.vector()->apply("insert");
    *file_badnode << u;

    value.resize(bad_cell_dofs.size(),1);
    u.vector()->zero();
    u.vector()->set(value.data(), bad_cell_dofs.size(),bad_cell_dofs.data());
    u.vector()->apply("insert");
    *file_badcell << u;

    value.resize(bad_cell_ov_dofs.size(),1);
    u.vector()->zero();
    u.vector()->set(value.data(), bad_cell_ov_dofs.size(),bad_cell_ov_dofs.data());
    u.vector()->apply("insert");
    *file_ov << u;

}

void NonlinearVariationalSolver::NonlinearCoarseDiscreteProblem::construct_coarse_IS(const GenericVector& res, double tol, GenericMatrix& _J, dolfin::la_index overlap)
{
    std::shared_ptr<const Mesh> _mesh = _problem->trial_space()->mesh();
    std::shared_ptr<const GenericDofMap> _dofmap = _problem->trial_space()->dofmap();

    clear_dofs_values();
    PETScVector is_dirichlet_dof(res.down_cast<PETScVector>());
    is_dirichlet_dof = 1;
    std::pair<std::size_t, std::size_t> range = is_dirichlet_dof.local_range();
    //info("vector local range %d-%d", range.first, range.second);
    // Iterate over cells
    for (CellIterator cell(*_mesh); !cell.end(); ++cell)
    {
        //info("iterating cells");
        const ArrayView<const dolfin::la_index> cell_dofs = _dofmap->cell_dofs(cell->index());
        bool added=false;
        for(std::size_t i =0; i<cell_dofs.size(); ++i)
        {
            //info("iterating cell dofs");
            double value, zero=0;
            //info("check dof(%d) in cell(%d) vs range (%d - %d)", cell_dofs[i], cell->index(), range.first, range.second);
            //if(cell_dofs[i]>= range.second- range.first)
            //    continue;
            res.get_local(&value, 1, &(cell_dofs[i]));
            //info("current value (%5.2f) vs tol (%5.2f)", value, tol);
            if(tol < std::abs(value)) 
            {
                bad_dofs.push_back(cell_dofs[i]);
                //info("found bad");
                if(!added)// found bad dof
                {
                    //info("adding");
                    //info("setting bad dofs");
                    for(std::size_t j =0; j <cell_dofs.size();++j)
                    {
                        is_dirichlet_dof.set_local(&zero, 1, &cell_dofs[j]);
                    }
                    added = true;
                }
            }
        }
    }
    is_dirichlet_dof.apply("insert");
    //info("add dofs");

    for(dolfin::la_index dof = range.first; dof < range.second; ++dof)
    { 
        double value;
        is_dirichlet_dof.get(&value, 1, &dof);
        if(std::abs(value) > 1./2) 
            dirichlet_dofs0.push_back(dof);
        else 
            bad_cell_dofs.push_back(dof);
    }

    info("#bad_cell_dofs: %d, #dirichlet_dofs %d", bad_cell_dofs.size(), dirichlet_dofs0.size());

    if(overlap)
    {
        PETScMatrix& J= _J.down_cast<PETScMatrix>();
        IS is;
        ISCreateGeneral(PETSC_COMM_WORLD, bad_cell_dofs.size(),bad_cell_dofs.data(),PETSC_COPY_VALUES,&is);
        MatIncreaseOverlap(J.mat(), 1, &is, overlap);
        PetscInt n;
        const PetscInt * nidx;
        ISGetLocalSize(is, &n);
        ISGetIndices(is, &nidx);

        for(PetscInt i = 0; i <n; i++)
        {
            double zero=0;
            is_dirichlet_dof.set(&zero, 1, &nidx[i]);
        }
        ISRestoreIndices(is,&nidx);
        ISDestroy(&is);
        info("#coarse_dofs with %d ov: %d", overlap, n);
    }
    is_dirichlet_dof.apply("insert");

    for(dolfin::la_index dof = range.first; dof < range.second; ++dof)
    { 
        double value;
        is_dirichlet_dof.get(&value,1,&dof);
        if(std::abs(value)>1./2) 
            add_dirichlet_dof(dof, 0.0);
        else bad_cell_ov_dofs.push_back(dof);
    }

}
