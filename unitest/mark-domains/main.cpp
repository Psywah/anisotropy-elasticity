// Copyright (C) 2007-2008 Anders Logg
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
// First added:  2007-04-24
// Last changed: 2011-01-25
//
// This demo program demonstrates how to mark sub domains
// of a mesh and store the sub domain markers as a mesh
// function to a DOLFIN XML file.
//
// The sub domain markers produced by this demo program
// are the ones used for the Stokes demo programs.

#include <dolfin.h>

using namespace dolfin;

int main()
{
  set_log_level(1);

#define R 0.01
#define Thk_med 0.00132
#define Thk_adv 0.00096
#define Depth 0.002
#define Thk_pla 0.006
#define H_ca ((Thk_pla-.002)/2.)
#define W_ca (H_ca*2.5)

  // adventitia domain (3)
  class Adventitia : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return sqrt(x[0]*x[0]+x[1]*x[1]) >= R+Thk_med - DOLFIN_EPS;
    }
  };

  // calcification domain (4)
  class Calcification : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        double cy = -R + Thk_pla/2.;
        double cxr = sqrt(W_ca*W_ca - H_ca*H_ca);
        double cxl = -cxr;
        return (sqrt((x[0]- cxl)*(x[0]-cxl)+(x[1]-cy)*(x[1]-cy))
               + sqrt((x[0]- cxr)*(x[0]-cxr)+(x[1]-cy)*(x[1]-cy)))
                <= (2*W_ca + DOLFIN_EPS);
    }
  };

  // Media domain (2)
  class Media : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return (sqrt(x[0]*x[0]+x[1]*x[1]) < R+Thk_med + DOLFIN_EPS) &&
               (sqrt(x[0]*x[0]+x[1]*x[1]) >= R -DOLFIN_EPS);
    }
  };

  // outer boundary (3)
  class OuterBnd : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return on_boundary && (sqrt(x[0]*x[0]+x[1]*x[1]) + DOLFIN_EPS*10 >=  R+Thk_med); 
    }
  }; 

  // Cross boundary boundary (2)
  class CrossBnd : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return (( x[2] + DOLFIN_EPS*10 >=  Depth) ||
               ( x[2] <= DOLFIN_EPS*10)) && on_boundary; 
    }
  }; 

  // onboundary boundary (2)
  class OnBnd : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return  on_boundary; 
    }
  }; 


  // Read mesh
  Mesh mesh("../../mesh/tube-4components-short.xml");

  // Create mesh functions over the cells and acets
  MeshFunction<std::size_t> sub_domains_mark(mesh, mesh.topology().dim() );
  MeshFunction<std::size_t> boundary_mark(mesh, mesh.topology().dim() - 1);

  // Mark all cells  as sub domain 1
  sub_domains_mark = 1;
  boundary_mark = 0;

  Calcification calcification;
  calcification.mark(sub_domains_mark, 4);
  Adventitia adventitia;
  adventitia.mark(sub_domains_mark, 3);
  Media media;
  media.mark(sub_domains_mark, 2);

  OnBnd on_bnd;
  on_bnd.mark(boundary_mark,1);
  OuterBnd out_bnd;
  out_bnd.mark(boundary_mark,3);
  CrossBnd cross_bnd;
  cross_bnd.mark(boundary_mark, 2);


  
  // Save sub domains to file
  File file("tube-4components-short-domains-marker.xml");
  file << sub_domains_mark;

  // Save sub domains to file
  File file_bnd("tube-4components-short-boundary-marker.xml");
  file_bnd << boundary_mark;
  

  
  

  
  
  /*
  sub_domains_mark = 0;
  boundary_mark = 0;

  // Read sub domains from file
  File file_read("../../mesh/tube-4components_physical_region.xml");
  file_read >> sub_domains_mark;

  // Save sub domains to file
  File file_bnd_read("../../mesh/tube-4components_facet_region.xml");
  file_bnd_read >> boundary_mark;
  */
  
  


  plot(mesh);
  plot(sub_domains_mark);
  plot(boundary_mark);
  interactive();


  return 0;


}
