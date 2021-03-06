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
/*
  // Read mesh
  Mesh mesh(MPI_COMM_WORLD, "/Users/shihua/BaiduYun/Note/anisotropic_hyperelasticity/code/data/CARE-566-left/carotid-artery/carotid_HII1.xml");

  //return 0;
  // Create mesh functions over the cells and acets
  MeshFunction<std::size_t> sub_domains_mark(reference_to_no_delete_pointer(mesh), mesh.topology().dim(), 1 );
  MeshFunction<std::size_t> boundary_mark(reference_to_no_delete_pointer(mesh), mesh.topology().dim() - 1);
  
  // Read sub domains from file
  File file_read("/Users/shihua/BaiduYun/Note/anisotropic_hyperelasticity/code/data/CARE-566-left/carotid-artery/carotid_HII1_physical_region.xml");
  file_read >> sub_domains_mark;

  // Save sub domains to file
  File file_bnd_read("/Users/shihua/BaiduYun/Note/anisotropic_hyperelasticity/code/data/CARE-566-left/carotid-artery/carotid_HII1_facet_region.xml");
  file_bnd_read >> boundary_mark;
  */

 
  // Mark all cells  as sub domain 1
  //sub_domains_mark = 1;
  //boundary_mark = 0;

  // Save sub domains to file
  //File file("carotid_HII_physical_region.xml");
  //file << sub_domains_mark;

  // Save sub domains to file
  //File file_bnd("carotid_HII_facet_region.xml");
  //file_bnd << boundary_mark;
  
  // serial write a patch of meshes
  
  //HDF5File filew(MPI_COMM_WORLD,"carotidHII.h5","w");
  /*
  HDF5File filew(MPI_COMM_WORLD,"pipe.h5","w");
  for(int i =1;i<3;i++)
  {
      // Read mesh
      std::string prefix("/Users/shihua/BaiduYun/Note/anisotropic_hyperelasticity/code/data/pipe/pipe");
      //std::string prefix("/Users/shihua/BaiduYun/Note/anisotropic_hyperelasticity/code/data/CARE-566-left/carotid-artery/carotid_HII");
      Mesh mesh(MPI_COMM_WORLD, prefix + std::to_string(i) + std::string(".xml"));

      // Create mesh functions over the cells and acets
      MeshFunction<std::size_t> sub_domains_mark(reference_to_no_delete_pointer(mesh), mesh.topology().dim(), 1 );
      MeshFunction<std::size_t> boundary_mark(reference_to_no_delete_pointer(mesh), mesh.topology().dim() - 1);

      // Read sub domains from file
      File file_read(prefix + std::to_string(i) + std::string("_physical_region.xml"));
      file_read >> sub_domains_mark;

      // Save sub domains to file
      File file_bnd_read(prefix + std::to_string(i) + std::string("_facet_region.xml"));
      file_bnd_read >> boundary_mark;


      std::string s1("/mesh");
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






// serial write single mesh
 /*
  HDF5File filew(MPI_COMM_WORLD,"carotidHII.h5","w");
  filew.write(mesh,"/mesh");
  filew.write(sub_domains_mark,"/subdomains_mark");
  filew.write(boundary_mark,"/facet_mark");
  filew.close();
  */




// parallel write

  Mesh mesh;
  //HDF5File filer(mesh.mpi_comm(),"carotidHII.h5","r");
  HDF5File filer(mesh.mpi_comm(),"pipe.h5","r");
  filer.read(mesh,"/mesh2",false);
  //filer.flush();

  MeshFunction<std::size_t> sub_domains_mark(reference_to_no_delete_pointer(mesh), mesh.topology().dim());
  MeshFunction<std::size_t> boundary_mark(reference_to_no_delete_pointer(mesh),mesh.topology().dim()-1);

  filer.read(sub_domains_mark,"/subdomains_mark2");
  std::cout<< "finished reading subdomains mark\n"<<std::flush<<std::endl;

  filer.read(boundary_mark,"facet_mark2");
  std::cout<< "finished reading\n"<<std::flush<<std::endl;
  info("has PETSc");
  

  //filer.flush();
/*
  HDF5File filew(MPI_COMM_WORLD,"carotidHIIdist.h5","w");
  filew.write(mesh,"mesh");
  std::cout<< "wrote mesh\n"<<std::endl;
  filew.write(sub_domains_mark,"subdomains_mark");
  std::cout<< "wrote subdomains mark\n"<<std::endl;
  filew.write(boundary_mark,"facet_mark");
  std::cout<< "wrote facet mark\n"<<std::endl;

  //filew.close();
  //filer.close();
*/





  
     plot(mesh);
     plot(sub_domains_mark);
     plot(boundary_mark);
     interactive();
     
     
     


  return 0;


}
