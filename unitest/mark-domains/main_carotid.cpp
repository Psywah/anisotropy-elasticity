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
  class OnBnd : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return  on_boundary; 
    }
  }; 

class LeftSide : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return (std::abs(x[0]) <  DOLFIN_EPS);
    }
  };

class RightSide : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return (std::abs(x[0] - Length) <  DOLFIN_EPS);
    }
  };
  */

  // Read mesh
  Mesh mesh("../../mesh/carotid_HII.xml");

  //return 0;
  // Create mesh functions over the cells and acets
  MeshFunction<std::size_t> sub_domains_mark(reference_to_no_delete_pointer(mesh), mesh.topology().dim(), 1 );
  MeshFunction<std::size_t> boundary_mark(reference_to_no_delete_pointer(mesh), mesh.topology().dim() - 1);

  // Mark all cells  as sub domain 1
  //sub_domains_mark = 1;
  //boundary_mark = 0;

  // Save sub domains to file
  //File file("carotid_HII_physical_region.xml");
  //file << sub_domains_mark;

  // Save sub domains to file
  //File file_bnd("carotid_HII_facet_region.xml");
  //file_bnd << boundary_mark;
  

  
  

  
  
  
  sub_domains_mark = 0;
  boundary_mark = 0;

  // Read sub domains from file
  File file_read("../../mesh/carotid_HII_physical_region.xml");
  file_read >> sub_domains_mark;

  // Save sub domains to file
  File file_bnd_read("../../mesh/carotid_HII_facet_region.xml");
  file_bnd_read >> boundary_mark;

  HDF5File filew(MPI_COMM_WORLD,"carotidHII.h5","w");
  filew.write(mesh,"mesh");
  filew.write(sub_domains_mark,"subdomains_mark");
  filew.write(boundary_mark,"facet_mark");
  filew.close();

  
  HDF5File filer(MPI_COMM_WORLD,"carotidHII.h5","r");
  filer.read(mesh,"mesh",false);
  filer.read(sub_domains_mark,"subdomains_mark");
  filer.read(boundary_mark,"facet_mark");
  filer.close();

  
  
  


  plot(mesh);
  plot(sub_domains_mark);
  plot(boundary_mark);
  interactive();


  return 0;


}
