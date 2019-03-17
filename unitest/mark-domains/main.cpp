// Copyright (C) 2007-2008 Anders Logg
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
//  DOLFIN is distributed in the hope that it will be useful,
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
#include <dolfin/io/HDF5File.h>
//#include "Cell.h"

using namespace dolfin;

int main()
{
  set_log_level(1);
/*
#define R 0.008
#define Thk_med 0.00132
#define Thk_adv 0.00096
//#define Depth 0.03
#define Depth 0.02
*/

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
        return (sqrt(x[0]*x[0]+x[1]*x[1]) <= R+Thk_med + DOLFIN_EPS) &&
               (sqrt(x[0]*x[0]+x[1]*x[1]) >= R -DOLFIN_EPS);
    }
  };

  // Plaque domain (1)
  class Plaque : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return (sqrt(x[0]*x[0]+x[1]*x[1]) <= R + DOLFIN_EPS);// &&
               //(sqrt(x[0]*x[0]+x[1]*x[1]) >= R -DOLFIN_EPS);
    }
  };

  // outer boundary (3)
  class OuterBnd : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return on_boundary && (sqrt(x[0]*x[0]+x[1]*x[1]) >=  R+Thk_med+Thk_adv - DOLFIN_EPS); 
    }
  }; 

  // inner boundary (3)
  class InnerBnd : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return on_boundary && (sqrt(x[0]*x[0]+x[1]*x[1]) <=  R+ DOLFIN_EPS); 
    }
  }; 

  // inlet boundary (2)
  class InletBnd : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return ( x[2] <= DOLFIN_EPS) && on_boundary;
    }
  };

  // outlet boundary (2)
  class OutletBnd : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return ( x[2] + DOLFIN_EPS >=  Depth)&& on_boundary;
    }
  }; 

  // Cross boundary boundary (2)
  class CrossBnd : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return (( x[2] + DOLFIN_EPS >=  Depth) ||
               ( x[2] <= DOLFIN_EPS)) && on_boundary; 
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

  HDF5File filew(MPI_COMM_WORLD,"tube-4components.h5","w");

//  for(int i=1;i<5;i++){
      // Read mesh
      //Mesh mesh("../../mesh/tube-2layerL2.xml");

  int i=4;
      Mesh mesh(MPI_COMM_WORLD, std::string("../../mesh/tube-4components") + std::to_string(i) + std::string(".xml"));

      info("mesh h:%f", mesh.hmax());
      info(mesh);
      //mesh.move(Constant(1.0,1.0,1.0));

      //return 0;
      //
      // Create mesh functions over the cells and acets
      MeshFunction<std::size_t> sub_domains_mark(reference_to_no_delete_pointer(mesh), mesh.topology().dim() );
      MeshFunction<std::size_t> boundary_mark(reference_to_no_delete_pointer(mesh), mesh.topology().dim() - 1);

      // Mark all cells  as sub domain 1
      //sub_domains_mark = 3;
      boundary_mark = 0;

      Media media;
      media.mark(sub_domains_mark, 2,false);
      Plaque plaque;
      plaque.mark(sub_domains_mark, 4,false);
      Calcification calcification;
      calcification.mark(sub_domains_mark, 1,false);
      Adventitia adventitia;
      adventitia.mark(sub_domains_mark, 3,false);

      InnerBnd inner_bnd;
      inner_bnd.mark(boundary_mark,3,false);
      OuterBnd outer_bnd;
      outer_bnd.mark(boundary_mark,4,false);
      InletBnd inlet_bnd;
      inlet_bnd.mark(boundary_mark,1,false);
      OutletBnd outlet_bnd;
      outlet_bnd.mark(boundary_mark,2,false);

      //HDF5File filew(MPI_COMM_WORLD,"tube-2layerL2.h5","w");
      //filew.write(mesh,"/mesh2");
      //filew.write(sub_domains_mark,"/subdomains_mark2");
      //filew.write(boundary_mark,"/facet_mark2");
/*      std::string s1("/mesh");
      s1 = s1 + std::to_string(i);
      filew.write(mesh,s1.c_str());

      std::string s2("/subdomains_mark");
      s2 = s2 + std::to_string(i);
      filew.write(sub_domains_mark,s2.c_str());

      std::string s3("/facet_mark");
      s3 = s3 + std::to_string(i);
      filew.write(boundary_mark,s3.c_str());

  }
  filew.close();
  */




  /*Plaque plaque;
    plaque.mark(sub_domains_mark, 1,false);
    Calcification calcification;
    calcification.mark(sub_domains_mark, 4,false);


    Mesh _mesh;
#define INDEX "2"
for(int i=0; i< 2;i ++)
{
_mesh  = refine(mesh);
mesh = _mesh;
}
*/


  // Save sub domains to file
  /*File file_mesh("tube-4components-short-refine" INDEX ".xml");
    file_mesh << mesh;
    File file("tube-4components-short-refine" INDEX "-domains-marker.xml");
    file << sub_domains_mark;

  // Save sub domains to file
  File file_bnd("tube-4components-short-refine" INDEX "-boundary-marker.xml");
  file_bnd << boundary_mark;
  */



  /* 
     Mesh mesh;
  //HDF5File filer(mesh.mpi_comm(),"carotidHII.h5","r");
  HDF5File filer(mesh.mpi_comm(),"tube-4components.h5","r");
  filer.read(mesh,"/mesh2",false);
  //filer.flush();

  MeshFunction<std::size_t> sub_domains_mark(reference_to_no_delete_pointer(mesh), mesh.topology().dim());
  MeshFunction<std::size_t> boundary_mark(reference_to_no_delete_pointer(mesh),mesh.topology().dim()-1);

  filer.read(sub_domains_mark,"/subdomains_mark2");
  std::cout<< "finished reading subdomains mark\n"<<std::flush<<std::endl;

  filer.read(boundary_mark,"facet_mark2");
  std::cout<< "finished reading\n"<<std::flush<<std::endl;
  */
  
  

  
  
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
