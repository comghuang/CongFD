#include "data.hpp"

void Data::cgnsoutputInit1D()
{
   cgsize_t isize[3][1];
   int ni,nj,nk,i,j,k;
   int index_file,icelldim,iphysdim,index_base;
   int index_zone,index_coord;
   char basename[33],zonename[33];

/* create gridpoints for simple example: */
   real* x=new double [n];
   for (ind i=0;i<n;i++) x[i]=data[i];

/* WRITE X, Y, Z GRID POINTS TO CGNS FILE */
/* open CGNS file for write */
   if (cg_open("grid_c.cgns",CG_MODE_WRITE,&index_file)) cg_error_exit();
/* create base (user can give any name) */
   strcpy(basename,"Base");
   icelldim=1;
   iphysdim=1;
   cg_base_write(index_file,basename,icelldim,iphysdim,&index_base);
/* define zone name (user can give any name) */
   strcpy(zonename,"Zone  1");
/* vertex size */
   isize[0][0]=n;
/* cell size */
   isize[1][0]=isize[0][0]-1;
/* boundary vertex size (always zero for structured grids) */
   isize[2][0]=0;
/* create zone */
   cg_zone_write(index_file,index_base,zonename,*isize,CGNS_ENUMV(Structured),&index_zone);
/* write grid coordinates (user must use SIDS-standard names here) */
   cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),"CoordinateX",
       x,&index_coord);
/* close CGNS file */
   delete[] x;
   cg_close(index_file);
}

void Data::oneDsolOutput(real t,std::string label)
{

    int ni,nj,nk,i,j,k,index_file,index_base,index_zone,index_flow,index_field;

    real* u=new real[n];

/* WRITE FLOW SOLUTION TO EXISTING CGNS FILE */
/* open CGNS file for modify */
    if (cg_open("grid_c.cgns",CG_MODE_MODIFY,&index_file)) cg_error_exit();
/* we know there is only one base (real working code would check!) */
    index_base=1;
/* we know there is only one zone (real working code would check!) */
    index_zone=1;
/* define flow solution node name (user can give any name) */
   std::string solname=label+" t="+std::to_string(t);
/* create flow solution node */
    int a=cg_sol_write(index_file,index_base,index_zone,solname.c_str(),CGNS_ENUMV(Vertex),&index_flow);
/* write flow solution (user must use SIDS-standard names here) */
   for(ind ivar=0;ivar<nVar;ivar++)
   {
      std::string varName;
      switch (ivar)
      {
      case 0: 
         varName="r";
         break;
      case 1: 
         varName="u";
         break;
      case 2: 
         varName="e";
         break;
      
      default:
         break;
      }

      for(ind i=0;i<n;i++)
      {
         u[i]=(*this)(i,ivar);
      }
      int a=cg_field_write(index_file,index_base,index_zone,index_flow,
        CGNS_ENUMV(RealDouble),varName.c_str(),u,&index_field);
   }

/* close CGNS file */
   delete[] u;
    cg_close(index_file);
}

void Data::cgnsoutputInit2D()
{
   cgsize_t isize[3][2];
   int ni,nj,nk,i,j,k;
   int index_file,icelldim,iphysdim,index_base;
   int index_zone,index_coord;
   char basename[33],zonename[33];

/* create gridpoints for simple example: */
   real* x=new double [n];
   for (ind i=0;i<n;i++) x[i]=data[i];

/* WRITE X, Y, Z GRID POINTS TO CGNS FILE */
/* open CGNS file for write */
   if (cg_open("grid_c.cgns",CG_MODE_WRITE,&index_file)) cg_error_exit();
/* create base (user can give any name) */
   strcpy(basename,"Base");
   icelldim=2;
   iphysdim=2;
   cg_base_write(index_file,basename,icelldim,iphysdim,&index_base);
/* define zone name (user can give any name) */
   strcpy(zonename,"Zone  1");
/* vertex size */
   isize[0][0]=n;
/* cell size */
   isize[1][0]=isize[0][0]-1;
/* boundary vertex size (always zero for structured grids) */
   isize[2][0]=0;
/* create zone */
   cg_zone_write(index_file,index_base,zonename,*isize,CGNS_ENUMV(Structured),&index_zone);
/* write grid coordinates (user must use SIDS-standard names here) */
   cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),"CoordinateX",
       x,&index_coord);
/* close CGNS file */
   delete[] x;
   cg_close(index_file);
}