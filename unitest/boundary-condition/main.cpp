// Copyright (C) 2009 Harish Narayanyan
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
// Modified by Anders Logg, 2011
//
// First added:  2009-09-29
// Last changed: 2012-11-12
//
// This demo program solves a hyperelastic problem

// Begin demo

#include <dolfin.h>
#include "HyperElasticity.h"

using namespace dolfin;

#define R 0.01
#define Thk_med 0.00132 
#define Thk_adv 0.00096
#define Depth 0.005

// Sub domain for clamp at left end
class LeftPoint : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    if( //(std::abs(x[2]) < DOLFIN_EPS) &&
           (std::abs(x[0]+R+Thk_med+Thk_adv) < DOLFIN_EPS) 
      )
    {   
        std::cout<<"found left point"<<std::endl;
        std::cout<<x[0]<<"  "<<x[1]<<"  "<<x[2]<<std::endl;
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
        std::cout<<"found right point"<<std::endl;
        std::cout<<x[0]<<"  "<<x[1]<<"  "<<x[2]<<std::endl;
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
        std::cout<<"found bottom point"<<std::endl;
        std::cout<<x[0]<<"  "<<x[1]<<"  "<<x[2]<<std::endl;
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

    const double pressure = 0.5;
    values[0] = -pressure*n[0];
    values[1] = -pressure*n[1];
    values[2] = -pressure*n[2];
  }

private:

  const Mesh& mesh;

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

  // Define solution function
  Function u(V);

  // Set material parameters
  const double E  = 10.0;
  const double nu = 0.3;
  Constant mu(E/(2*(1 + nu)));
  Constant lambda(E*nu/((1 + nu)*(1 - 2*nu)));

  // Create (linear) form defining (nonlinear) variational problem
  HyperElasticity::ResidualForm F(V);
  F.mu = mu; F.lmbda = lambda; F.u = u;
  F.B = B; F.T = pressure_normal;
  F.ds = boundary_mark;

  // Create jacobian dF = F' (for use in nonlinear solver).
  HyperElasticity::JacobianForm J(V, V);
  J.mu = mu; J.lmbda = lambda; J.u = u;

  // Solve nonlinear variational problem F(u; v) = 0
  solve(F == 0, u, bcs, J);

  // Save solution in VTK format
  File file("displacement.pvd");
  file << u;

  // Plot solution
  plot(u);
  interactive();

  
  return 0;
}
