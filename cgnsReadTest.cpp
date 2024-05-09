/*   Program read_grid2zn_str   */
/*
Reads simple 3-D structured 2-zone grid from CGNS file
(companion program to write_grid2zn_str.c)

The CGNS grid file 'grid_c.cgns' must already exist.

Example compilation for this program is (change paths if needed!):

cc -I ../.. -c read_grid2zn_str.c
cc -o read_grid2zn_str_c read_grid2zn_str.o -L ../../lib -lcgns

(../../lib is the location where the compiled
library libcgns.a is located)
*/

#include <stdio.h>
/* cgnslib.h file must be located in directory specified by -I during compile: */
#include "cgnslib.h"
#include "data.hpp"

#if CGNS_VERSION < 3100
# define cgsize_t int
#endif

int main()
{
   int ifile,index_base,nzone,i,j,k;
   int n1to1=0;
   char zonename[33],basename[33];

   std::string a;
   int cell_dim=2,phys_dim=2;
   std::vector<Data> xs,ys;


    if (cg_open("testgrid.cgns",CG_MODE_READ,&ifile)) cg_error_exit();
   
   cg_base_read(ifile, index_base, basename,&cell_dim, &phys_dim);

   cgsize_t isize[3][cell_dim];
   cgsize_t irmin[cell_dim],irmax[cell_dim];

   int itrans[cell_dim];
   char connectname[33],donorname[33];
   cgsize_t ipnts[2][cell_dim],ipntsdonor[2][cell_dim];

   index_base=1;
    cg_nzones(ifile,index_base,&nzone);
    xs.resize(nzone);
    ys.resize(nzone);

    for (int index_zone=1; index_zone <= nzone; index_zone++)
    {

      std::cout<<"\n\n";
      if(cg_zone_read(ifile,index_base,index_zone,zonename,isize[0])) cg_error_exit();
      a=zonename;
      std::cout<<a<<'\n';

      irmin[0]=1;
      irmin[1]=1;

      irmax[0]=isize[0][0];
      irmax[1]=isize[0][1];
      
      int len=irmax[0]*irmax[1];
      //real *x=new real[len];
      //real *y=new real[len];
      real x[len*2],y[len*2];
      cg_coord_read(ifile,index_base,index_zone,"CoordinateX",CGNS_ENUMV(RealDouble),irmin,irmax,x);
      cg_coord_read(ifile,index_base,index_zone,"CoordinateY",CGNS_ENUMV(RealDouble),irmin,irmax,y);
      xs[index_zone-1].init(len,1);
      xs[index_zone-1].setValue(x,len);
      ys[index_zone-1].init(len,1);
      ys[index_zone-1].setValue(y,len);

      if(cg_n1to1(ifile, index_base, index_zone, &n1to1)) cg_error_exit();
      std::cout<<"n1to1=  "<<n1to1<<'\n';
      for (int i1to1 = 1; i1to1 <= n1to1; i1to1++)
      {
         if(cg_1to1_read(ifile,index_base,index_zone,i1to1,connectname,donorname,ipnts[0],ipntsdonor[0],itrans)) cg_error_exit();
         std::cout<<connectname<<'\n';
         std::cout<<ipnts[0][0]<<' '<<ipnts[0][1]<<" and "<<ipnts[1][0]<<' '<<ipnts[1][1]<<'\n';
      }
      
      std::cout<<"cg_1to1_read finish \n";

      int nbocos,normallist,ndataset;
      int normalindex[cell_dim];
      char boconame[33];
      CGNS_ENUMT(BCType_t) ibocotype;
      CGNS_ENUMT(PointSetType_t) iptset;
      CGNS_ENUMT(DataType_t) normaldatatype;
      cgsize_t npts,normallistflag;


      cg_nbocos(ifile,index_base,index_zone,&nbocos);
      std::cout<<"nbocos== "<<nbocos<<'\n';
      for(int ib=1;ib<=nbocos;ib++)
      {
         cg_boco_info(ifile,index_base,index_zone,ib,boconame,&ibocotype,
                   &iptset,&npts,normalindex,&normallistflag,&normaldatatype,&ndataset);
         std::cout<<BCTypeName[ibocotype]<<"   "<<PointSetTypeName[iptset]<<'\n';


         cgsize_t ipnts[npts][cell_dim];
         cg_boco_read(ifile,index_base,index_zone,ib,ipnts[0],&normallist);
         std::cout<<ipnts[0][0]<<' '<<ipnts[0][1]<<" and ";
         std::cout<<ipnts[1][0]<<' '<<ipnts[1][1]<<'\n';
      }

      
      
    }

    cg_close(ifile);
    return 0;
}
