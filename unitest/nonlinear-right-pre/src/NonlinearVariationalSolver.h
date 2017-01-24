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
// Modified by Marie E. Rognes, 2011.
// Modified by Corrado Maurini, 2013.
//
// First added:  2011-01-14 (2008-12-26 as VariationalProblem.h)
// Last changed: 2013-11-20

#ifndef __NONLINEAR_VARIATIONAL_SOLVER_H
#define __NONLINEAR_VARIATIONAL_SOLVER_H

#include <dolfin/nls/NonlinearProblem.h>
#include <dolfin/nls/NewtonSolver.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/GenericDofMap.h>
#include <petscmat.h>
#include "PETScSNESSolver.h"
#include "NonlinearVariationalProblem.h"
#include "ColorIO.h"

namespace dolfin
{

  /// This class implements a solver for nonlinear variational problems.

  class NonlinearVariationalSolver : public Variable
  {
  public:

    /// Create nonlinear variational solver for given problem
    NonlinearVariationalSolver(NonlinearVariationalProblem& problem);

    /// Create nonlinear variational solver for given problem (shared
    /// pointer version)
    NonlinearVariationalSolver(std::shared_ptr<NonlinearVariationalProblem> problem);

    /// Solve variational problem with bound constraints defined by
    /// GenericVectors
    ///
    /// *Arguments*
    ///     lb (_GenericVector_)
    ///         The linear solver.
    ///     ub (_GenericVector_)
    ///         The factory.
    /// *Returns*
    ///     std::pair<std::size_t, bool>
    ///         Pair of number of Newton iterations, and whether
    ///         iteration converged)
    std::pair<std::size_t, bool> solve(const GenericVector& lb,
                                       const GenericVector& ub);

    /// Solve variational problem with bound constraints defined by
    /// GenericVectors (shared pointer version)
    ///
    /// *Arguments*
    ///     lb (_std::shared_ptr<const GenericVector>_)
    ///         The linear solver.
    ///     ub (_std::shared_ptr<const GenericVector>_)
    ///         The factory.
    /// *Returns*
    ///     std::pair<std::size_t, bool>
    ///         Pair of number of Newton iterations, and whether
    ///         iteration converged)
    std::pair<std::size_t, bool>
      solve(std::shared_ptr<const GenericVector> lb,
            std::shared_ptr<const GenericVector> ub);

    /// Solve variational problem with bound constraints defined by Functions
    ///
    /// *Arguments*
    ///     lb (_Function_)
    ///         The linear solver.
    ///     ub (_Function_)
    ///         The factory.
    /// *Returns*
    ///     std::pair<std::size_t, bool>
    ///         Pair of number of Newton iterations, and whether
    ///         iteration converged)
    std::pair<std::size_t, bool> solve(const Function& lb,
                                       const Function& ub);

    /// Solve variational problem with bound constraints defined by
    /// Functions (shared pointer version)
    ///
    /// *Arguments*
    ///     lb (_std::shared_ptr<const Function>_)
    ///         The linear solver.
    ///     ub (_std::shared_ptr<const Function>_)
    ///         The factory.
    /// *Returns*
    ///     std::pair<std::size_t, bool>
    ///         Pair of number of Newton iterations, and whether
    ///         iteration converged)
    std::pair<std::size_t, bool> solve(std::shared_ptr<const Function> lb,
                                       std::shared_ptr<const Function> ub);

    /// Solve variational problem
    ///
    /// *Returns*
    ///     std::pair<std::size_t, bool>
    ///         Pair of number of Newton iterations, and whether
    ///         iteration converged)
    std::pair<std::size_t, bool> solve();

    /// Default parameter values
    static Parameters default_parameters()
    {
      Parameters p("nonlinear_variational_solver");

      p.add("symmetric", false);
      p.add("print_rhs", false);
      p.add("print_matrix", false);
      p.add("NL_atol", 1e-8);
      p.add("NL_rtol", 1e-2);
      p.add("NL_res_r", 1e-1);
      p.add("NL_infty_r", 1e-1);
      p.add("NL_size_r", 1e-1);
      p.add("NL_overlap", 0);
      p.add("NL_report", true);
      p.add("dt", 1e-1);

      std::set<std::string> nonlinear_solvers = {"newton"};
      std::string default_nonlinear_solver = "newton";
      p.add(NewtonSolver::default_parameters());

      #ifdef HAS_PETSC
      p.add(PETScSNESSolver::default_parameters());
      nonlinear_solvers.insert("snes");
      #endif

      p.add("nonlinear_solver", default_nonlinear_solver, nonlinear_solvers);

      return p;
    }

  private:

    // Nonlinear (algebraic) problem
    class NonlinearDiscreteProblem : public NonlinearProblem
    {
    public:

      // Constructor
      NonlinearDiscreteProblem(
        std::shared_ptr<NonlinearVariationalProblem> problem,
        std::shared_ptr<NonlinearVariationalSolver> solver);

      // Destructor
      ~NonlinearDiscreteProblem();

      // Compute F at current point x
      virtual void F(GenericVector& b, const GenericVector& x);

      // Compute J = F' at current point x
      virtual void J(GenericMatrix& A, const GenericVector& x);

    private:

      // Problem and solver objects
      std::shared_ptr<NonlinearVariationalProblem> _problem;
      std::shared_ptr<NonlinearVariationalSolver> _solver;

    };

    // Nonlinear coarse (algebraic) problem
    class NonlinearCoarseDiscreteProblem : public NonlinearProblem
    {
    public:

      // Constructor
      NonlinearCoarseDiscreteProblem(
        std::shared_ptr<NonlinearVariationalProblem> problem,
        std::shared_ptr<NonlinearVariationalSolver> solver);

      // Destructor
      ~NonlinearCoarseDiscreteProblem();

      // Compute F at current point x
      virtual void F(GenericVector& b, const GenericVector& x);

      // Compute J = F' at current point x
      virtual void J(GenericMatrix& A, const GenericVector& x);

      void add_dirichlet_dof(dolfin::la_index dof, double value)
      {
          dirichlet_dofs.push_back(dof); values.push_back(value);
      }


      void clear_dofs_values()
      {
          dirichlet_dofs.clear(); 
          values.clear(); 
          bad_dofs.clear(); 
          dirichlet_dofs0.clear();
          bad_cell_dofs.clear();
          bad_cell_ov_dofs.clear();
      }


      void set_all_dofs(GenericVector& u)
      {
          clear_dofs_values();
          std::pair<std::size_t, std::size_t> range = u.local_range();
          for(dolfin::la_index dof = range.first; dof<range.second; dof++)
          {
              double value;
              u.get(&value, 1, &dof);
              //add_good_dof(dof, value);
              add_dirichlet_dof(dof, 0.0);
          }
          dirichlet_dofs0 = dirichlet_dofs;
      }

      void construct_coarse_IS(const GenericVector& res, double tol, GenericMatrix& _J, dolfin::la_index overlap=0)
      {
          std::shared_ptr<const Mesh> _mesh = _problem->trial_space()->mesh();
          std::shared_ptr<const GenericDofMap> _dofmap = _problem->trial_space()->dofmap();

          clear_dofs_values();
          std::vector<bool> is_dirichlet_dof(_dofmap->global_dimension(), true);
          // Iterate over cells
          for (CellIterator cell(*_mesh); !cell.end(); ++cell)
          {
              //info("iterating cells");
              const ArrayView<const dolfin::la_index> cell_dofs = _dofmap->cell_dofs(cell->index());
              bool added=false;
              for(std::size_t i =0; i<cell_dofs.size(); ++i)
              {
                  //info("iterating cell dofs");
                  double value;
                  res.get(&value, 1, &(cell_dofs[i]));
                  if(tol < std::abs(value)) 
                  {
                      bad_dofs.push_back(cell_dofs[i]);
                      if(!added)// found bad dof
                      {
                          //info("setting bad dofs");
                          for(std::size_t j =0; j <cell_dofs.size();++j)
                          {is_dirichlet_dof[cell_dofs[j]] = false;}
                          added = true;
                      }
                  }
              }
          }
          //info("add dofs");

          for(dolfin::la_index dof = 0; dof <is_dirichlet_dof.size(); ++dof)
          { if(is_dirichlet_dof[dof]) 
                dirichlet_dofs0.push_back(dof);
            else 
                bad_cell_dofs.push_back(dof);
          }
          info("#bad_cell_dofs: %d, #dirichlet_dofs %d", bad_cell_dofs.size(), dirichlet_dofs0.size());

          if(overlap)
          {
              PETScMatrix& J= _J.down_cast<PETScMatrix>();
              IS is;
              ISCreateGeneral(PETSC_COMM_SELF, bad_cell_dofs.size(),bad_cell_dofs.data(),PETSC_COPY_VALUES,&is);
              MatIncreaseOverlap(J.mat(), 1, &is, overlap);
              PetscInt n;
              const PetscInt * nidx;
              ISGetLocalSize(is, &n);
              ISGetIndices(is, &nidx);

              for(PetscInt i = 0; i <n; i++)
              {
                  is_dirichlet_dof[nidx[i]] = false;
              }
              ISRestoreIndices(is,&nidx);
              ISDestroy(&is);
              info("#coarse_dofs with %d ov: %d", overlap, n);
          }

          for(dolfin::la_index dof = 0; dof <is_dirichlet_dof.size(); ++dof)
          { if(is_dirichlet_dof[dof]) 
              add_dirichlet_dof(dof, 0.0);
          else bad_cell_ov_dofs.push_back(dof);}

      }

      std::size_t size_dirichlet_dofs()
      {
          return dirichlet_dofs.size();
      }

      std::size_t size_bad_dofs()
      {
          return bad_dofs.size();
      }

      void restricted_update(GenericVector& u, const GenericVector& u_pre)
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

      void save_coarse(std::shared_ptr<File> file_badnode, std::shared_ptr<File> file_badcell, std::shared_ptr<File> file_ov)
      {
          Function u(_problem->trial_space());
          std::vector<double> value;
          value.resize(bad_dofs.size(),1);
          u.vector()->set(value.data(), bad_dofs.size(),bad_dofs.data());
          *file_badnode << u;

          value.resize(bad_cell_dofs.size(),1);
          u.vector()->set(value.data(), bad_cell_dofs.size(),bad_cell_dofs.data());
          *file_badcell << u;

          value.resize(bad_cell_ov_dofs.size(),1);
          u.vector()->set(value.data(), bad_cell_ov_dofs.size(),bad_cell_ov_dofs.data());
          *file_ov << u;

      }

    private:

      // Problem and solver objects
      std::shared_ptr<NonlinearVariationalProblem> _problem;
      std::shared_ptr<NonlinearVariationalSolver> _solver;
      std::vector<dolfin::la_index> dirichlet_dofs;
      std::vector<dolfin::la_index> dirichlet_dofs0;
      std::vector<dolfin::la_index> bad_dofs;
      std::vector<dolfin::la_index> bad_cell_dofs;
      std::vector<dolfin::la_index> bad_cell_ov_dofs;
      std::vector<double> values;

    };

    // The nonlinear problem
    std::shared_ptr<NonlinearVariationalProblem> _problem;

    // The nonlinear discrete problem
    std::shared_ptr<NonlinearDiscreteProblem> nonlinear_problem;
    std::shared_ptr<NonlinearCoarseDiscreteProblem> nonlinear_coarse_problem;

    // The Newton solver
    std::shared_ptr<NewtonSolver> newton_solver;

#ifdef HAS_PETSC
    // Or, alternatively, the SNES solver
    std::shared_ptr<PETScSNESSolver> snes_solver;
#endif

  };

}

#endif
