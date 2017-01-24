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

#include <dolfin/common/NoDeleter.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/function/Function.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/fem/Form.h>
#include <dolfin/la/PETScMatrix.h>
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
       File r_hist("r_hist.pvd");
       info("begin snes solve with nonlinear elimination");
       Parameters para_coarse = parameters("snes_solver");
       para_coarse["absolute_tolerance"] = (double)(parameters["NL_atol"]);
       para_coarse["relative_tolerance"] = (double)(parameters["NL_rtol"]);
       para_coarse["report"] = (bool)(parameters["NL_report"]);
       double rho0 = parameters["NL_res_r"];
       double rho1 = parameters["NL_infty_r"];
       double rho2 = parameters["NL_size_r"];
       int overlap = parameters["NL_overlap"];

       Function res(u->function_space());
       double norm_res0,norm_res;
       nonlinear_problem->F(*res.vector(), *u->vector());
       norm_res0 =res.vector()->norm("l2");

       bool converged = false;
       std::size_t iter = 0, max_iter = parameters("snes_solver")["maximum_iterations"] ;
       parameters("snes_solver")["maximum_iterations"] = 1;
       while(iter< max_iter)
       {
           snes_solver->parameters.update(parameters("snes_solver"));
           ret = snes_solver->solve(*nonlinear_problem, *u->vector());
           ++iter;
           converged = std::get<1>(ret);
           nonlinear_problem->F(*res.vector(), *u->vector());
           norm_res =res.vector()->norm("l2");
           info("%d iterations of global nonlinear iteration, last residual:%5.2e, vs current residual:%5.2e",
                   iter, norm_res0, norm_res);
           if(converged)
               break;
           begin("Check Nonlinearity.");
           if(norm_res < norm_res0*rho0)
           {
               info("Reduction of residual %.2f < rho0(%.2f), No Need for nonlinear elimination",
                       norm_res/norm_res0, rho0);
               norm_res0 = norm_res;
               end();
               continue;
           }
           info("Reduction of residual %.2f > rho0(%.2f), Need for nonlinear elimination",
                   norm_res/norm_res0, rho0);
           // find all good dof
           double infty = res.vector()->norm("linf");
           info("finding bad dofs [ > %f*%.2f (infty norm * rho1)]", infty, rho1);
           PETScMatrix J;
           nonlinear_problem->J(J,*u->vector());
           nonlinear_coarse_problem->construct_coarse_IS(*res.vector(), infty*rho1, J, overlap);
           std::size_t size_total = res.vector()->size();
           std::size_t size_good = nonlinear_coarse_problem->size_dirichlet_dofs();

           if( size_total - size_good > rho2*size_total)
           {
               info("Size of coarse problem %d > %d*%.2f (total_size * rho2), No Need for nonlinear elimination",
                       size_total-size_good, size_total, rho2);
               norm_res0 = norm_res;
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

           r_hist << res ;
           nonlinear_problem->F(*res.vector(), *u->vector());
           r_hist << res ;
           double norm_res_c = res.vector()->norm("l2");
           if(norm_res_c<norm_res)
           {
               info("Residual norm with correction %f < %f, accept nonlinear elimination",
                       norm_res_c,norm_res);
               norm_res0 = norm_res_c;
           }
           else {
               info("Residual norm with correction %f > %f, not accept nonlinear elimination",
                       norm_res_c,norm_res);
               //norm_res0 =norm_res;
               //*u->vector() = x_copy;
               //norm_res0 = norm_res_c;
           }
           info("Solved Coarse Problem");
           end();
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
  b.set_local(values.data(), values.size(), dirichlet_dofs.data());
  b.apply("insert");

  // Print vector
  dolfin_assert(_solver);
  const bool print_rhs = _solver->parameters["print_rhs"];
  if (print_rhs)
    info(b, true);
}
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
  A.ident_local(dirichlet_dofs.size(),dirichlet_dofs.data());
  A.apply("insert");

  // Print matrix
  dolfin_assert(_solver);
  const bool print_matrix = _solver->parameters["print_matrix"];
  if (print_matrix)
    info(A, true);
}
//-----------------------------------------------------------------------------

