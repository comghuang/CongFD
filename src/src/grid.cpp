#include"block.hpp"

void Block::initUniform(std::array<int,3> iMax_,std::array<double,6> calZone)
{
    iMax=iMax_;
    dim=0;
    for (int i = 0; i < 3; i++)
    {
        if (iMax[i]>2) dim++;
        else 
        {
            iMax[i]=2;
        }
        icMax[i]=iMax[i]-1;
    }
    nVer=1,nCel=1;
    for (int i = 0; i < 3; i++)
    {
        nVer*=iMax[i];
        nCel*=icMax[i];
    }
    coorVer.init(nVer,dim);
    coorCel.init(nCel,dim);
    for (int idim = 0; idim < dim; idim++)
    {
        double cmin=calZone[idim*2];
        double cmax=calZone[idim*2+1];
        double interval=(cmax-cmin)/(iMax[idim]-1);

        int l,m,n;
        int* onedIndex=((idim == 0 ) ? &l : ( idim == 1 ? &m : &n));
        
        for (l = 0; l < iMax[0]; l++)
        for (m = 0; m < iMax[1]; m++)
        for (n = 0; n < iMax[2]; n++)
        {
            int globalIndex=l+m*iMax[0]+n*iMax[0]*iMax[1];
            coorVer(globalIndex,idim)=cmin+(*onedIndex)*interval;
        }
    }
    
    //for cellcenter
    for (int idim = 0; idim < dim; idim++)
    {
        int l,m,n,iLen=(dim==1?2:dim==2? 4:8);
        std::vector<int> index;
        index.resize(iLen);
        for (l = 0; l < icMax[0]; l++)
        for (m = 0; m < icMax[1]; m++)
        for (n = 0; n < icMax[2]; n++)
        {
            int iVerGlobal=l+m*iMax[0]+n*iMax[0]*iMax[1];
            int iCelGlobal=l+m*icMax[0]+n*icMax[0]*icMax[1];
            double temp=0;
            index[0]=iVerGlobal;
            index[1]=iVerGlobal+1;
            if (dim>=2)
            {
                index[2]=iVerGlobal+iMax[0];
                index[3]=iVerGlobal+iMax[0]+1;
            }
            if (dim>=3)
            {
                index[4]=iVerGlobal+iMax[0]*iMax[1];
                index[5]=iVerGlobal+1+iMax[0]*iMax[1];
                index[6]=iVerGlobal+iMax[0]+iMax[0]*iMax[1];
                index[7]=iVerGlobal+iMax[0]+1+iMax[0]*iMax[1];
            }

            for(auto iver:index) 
            temp+=coorVer(iver,idim);
            
            
            coorCel(iCelGlobal,idim)=temp/iLen;
            
        }
    }
    


    
    
}

void Block::outputCgns()
{
   cgsize_t isize[3][dim];
   int ni,nj,nk,i,j,k;
   int index_file,icelldim,iphysdim,index_base;
   int index_zone,index_coord;
   char basename[33],zonename[33];

/* create gridpoints for simple example: */
   std::vector<real> x;
   x.reserve(nVer);

/* WRITE X, Y, Z GRID POINTS TO CGNS FILE */
/* open CGNS file for write */
   if (cg_open("grid_c.cgns",CG_MODE_WRITE,&index_file)) cg_error_exit();
/* create base (user can give any name) */
   strcpy(basename,"Base");
   icelldim=dim;
   iphysdim=dim;
   cg_base_write(index_file,basename,icelldim,iphysdim,&index_base);
/* define zone name (user can give any name) */
   strcpy(zonename,"Zone  1");
   for (int idim = 0; idim < dim; idim++)
   {
        /* vertex size */
        isize[0][idim]=iMax[idim];
        /* cell size */
        isize[1][idim]=isize[0][idim]-1;
        /* boundary vertex size (always zero for structured grids) */
        isize[2][idim]=0;
   }
   

/* create zone */
   cg_zone_write(index_file,index_base,zonename,*isize,CGNS_ENUMV(Structured),&index_zone);
/* write grid coordinates (user must use SIDS-standard names here) */
    for (int idim = 0; idim < dim; idim++)
    {
        std::string name=(idim==0)? "CoordinateX":(idim==1)? "CoordinateY":"CoordinateZ";
        for (int i=0 ; i<nVer ; i++ ) x[i]=coorVer(i,idim);
        cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),
        name.c_str(),
        x.data(),&index_coord);
    }
    
   
/* close CGNS file */
   cg_close(index_file);
}


real Block::operator()(int ic,int idim)
{
    return coorCel(ic,idim);
}