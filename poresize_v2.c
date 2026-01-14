/**
   @Program : poresizedistribution.c
   @Author  : Himangsu Bhaumik hb685@cam.ac.uk
   @Note    : PSD calculation Gelb and Gubbins Langmuir 1999, 15, 305-308
   @Compilation : gcc poresizedistribution.c  -lm
   
   @executuon  : ./a.out restart-data-file
   @output: 

   @version: Here I will find dr void cell to atom position
             Explicitely check if void is within atom? by checking dr<radii
             grid size can be made smaller.
   
   @Note : Be careful about the number of atoms defined at the
   begining, if you change the definition of neighbour the formula Rc=...
   should be changed!
   
   Reading lammps restart file: it reads lammps restart_data file
   which has some specific structure, for other kind of file one
   should modify the reading section. 
   
   
   


**/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#define N 10000
#define maxcell 80




int main(int argc, char *argv[])     
{
  FILE *fp,*fp1;
  static char FILE[1024],FILE1[1024];

  static char *e0;

  static double L,iL,Lf,Li,binw;

  static long i,j,k,m,set;
  
  static double x_cor,y_cor,z_cor,vx,vy,vz,tild,tild0,tild1;
  static double xr,yr,zr,dr,Rc=1.568;

  static int Nr,ipart,ineighbour,ilist,id,itype,nx,ny,nz;
  static double dx1,dx2;

  double x[N][3];
  int type[N],nlist[N];
  int indx[3],cell[maxcell][maxcell][maxcell];
  double cellD[maxcell][maxcell][maxcell];

  double radii[7]={0.88,0.92,0.96,1.0,1.04,1.08,1.12};
  
  double Dbin[1000];
  for(i=0;i<1000;i++)Dbin[i]=0;
  
  if(argc!=2){
    printf("Give filename as arguments\n");
    return 0;
    }
    
  int maxidx=0;double totcount=0;
  double sumD=0,Dcount=0;
  
    
    sprintf(FILE,"%s",argv[1]);
    fp=fopen(FILE,"r");
    printf("file open %s\n",FILE);
    
     while (!feof(fp)){
  	fscanf(fp,"%*[^\n]\n");
  	fscanf(fp,"%d\n",&Nr);
  	fscanf(fp,"%*[^\n]\n");
  	fscanf(fp,"%*[^\n]\n");
  	fscanf(fp,"%lf %lf \n",&Li,&Lf);L=Lf-Li;
  	for(i=0;i<21;i++)fscanf(fp,"%*[^\n]\n");
  	
  	for(i=0;i<N;i++){
  	  fscanf(fp,"%d %d %lf %lf %lf %d %d %d\n",&id,&itype,&x_cor,&y_cor,&z_cor,&nx,&ny,&nz);
  	  id-=1;
  	  x[id][0]=x_cor+(nx*L);
	  x[id][1]=y_cor+(ny*L);
	  x[id][2]=z_cor+(nz*L);
	  type[id]=itype;
  	}
  	for(i=0;i<=N;i++)fscanf(fp,"%*[^\n]\n");
      }
      fclose(fp);

      if(N!=Nr){
	printf("Something wrong with defined N and N that read from file\n");
	return 0;
      }
	
      
    double dl=L/maxcell;

    //i=518;
    //printf("%lf %lf %lf \n",x[i][0],x[i][1],x[i][2]);
    
    printf("File read completed\n");
    printf("L=%f\n",L);
    //=======================================================================
    
      iL=1.0/L;

      //==========Transform co-ordinate to wrapped one========== 
      for(i=0;i<N;i++){
	for(j=0;j<3;j++){
	 while(x[i][j]<0)x[i][j]+=L;
	 while(x[i][j]>L)x[i][j]-=L;
       }
      }
       //========================================================= 
      printf("Wrapped coordinate\n");

      for(i=0;i<N;i++)nlist[i]=0;
      
      for(ipart=0;ipart<N-1;ipart++){
	
	for(ineighbour=ipart+1;ineighbour<N;ineighbour++){

	  xr=x[ipart][0]-x[ineighbour][0];
	  yr=x[ipart][1]-x[ineighbour][1];
	  zr=x[ipart][2]-x[ineighbour][2];
	  xr=xr-tild0*round(yr*iL);
  	  xr=xr-L*round(xr*iL);
  	  yr=yr-L*round(yr*iL);
  	  zr=zr-L*round(zr*iL);
  	  dr=sqrt(xr*xr+yr*yr+zr*zr);
	  
	  Rc=(0.88+(type[ipart]-1)*0.02+(type[ineighbour]-1)*0.02)*1.4;
	  if(dr<Rc){
	    nlist[ipart]++;
	    nlist[ineighbour]++;
	  }
	}
      }//ipart
      //=========================================
      printf("Neighbourlist done\n");
      
      // discretize in cell: w/ particle 1, w/o 0

      for(i=0;i<maxcell;i++)
	for(j=0;j<maxcell;j++)
	  for(k=0;k<maxcell;k++){
	    cell[i][j][k]=0;
	    cellD[i][j][k]=0;
	  }
      
      for(i=0;i<N;i++){
	if(nlist[i]>0)// avoid lonely atom
	  {
	    for(j=0;j<3;j++)
	      indx[j]=x[i][j]/dl;
	    
	    cell[indx[0]][indx[1]][indx[2]]=1;
	    
	  }
      }
      printf("Discretize cell and avoid lonely atom\n");

      




      
      for(i=0;i<maxcell;i++){
	for(j=0;j<maxcell;j++){
	  for(k=0;k<maxcell;k++){

	    //printf("loop over cell %d %d %d\n",i,j,k);
	    
	    if(cell[i][j][k]==0){
	      double xi=i*dl,yi=j*dl,zi=k*dl;

	       double mindr=L;
	      for(ipart=0;ipart<N;ipart++){
		if(nlist[ipart]){ //avoid lonely atom flying in space
		  
		  xr=xi-x[ipart][0];
		  yr=yi-x[ipart][1];
		  zr=zi-x[ipart][2];
		  
		  xr=xr-L*round(xr*iL);
		  yr=yr-L*round(yr*iL);
		  zr=zr-L*round(zr*iL);
		  
		  dr=sqrt(xr*xr+yr*yr+zr*zr);
		  

		  if(dr<radii[type[ipart]-1]){
		    cell[i][j][k]=1;
		    goto skipcell;
		  }
		  if(dr<mindr)
		    mindr=dr;
		}
	      }



	      dr=mindr;//funny tranformation, huh!??


	      int ir=dr/dl;
	      //for(int ir=1;ir<maxcell/2;ir++){
	      //printf("ir assiging loop dr=%e\n",dr);
		for(int ii=i-ir;ii<i+ir;ii++){
		  int id=(ii+maxcell)%maxcell;
		  for(int jj=j-ir;jj<j+ir;jj++){
		    int jd=(jj+maxcell)%maxcell;
		    for(int kk=k-ir;kk<k+ir;kk++){
		      int kd=(kk+maxcell)%maxcell;

		      //Here we need to check the distance of each cell
		      //if they are within dr we will take them
		      xr=(i-id)*dl;
		      yr=(j-jd)*dl;
		      zr=(k-kd)*dl;
		      xr=xr-L*round(xr*iL);
		      yr=yr-L*round(yr*iL);
		      zr=zr-L*round(zr*iL);
		      double ds=sqrt(xr*xr+yr*yr+zr*zr);
		      if(ds<dr)
		      
		      if(cellD[id][jd][kd]<dr)
			cellD[id][jd][kd]=dr;
		      
		      		      
		    }//kk
		  }//jj
		}//ii
		//}//ir
	      
	      
	      
	      
	      
	      
	      
	    }//if cell==0
	  skipcell:;
	  }//k
	}//j
      }//i
      
      
      binw=2*dl;
      for(i=0;i<maxcell;i++)
	for(j=0;j<maxcell;j++)
	  for(k=0;k<maxcell;k++)
	    if(cell[i][j][k]==0){
	      //printf("%e\n",cellD[i][j][k]*dl);
	      double Dia=2*cellD[i][j][k];
	      sumD+=Dia;
	      Dcount++;
	      int idx=Dia/binw; //factor 2 for diameter
	      if(idx>maxidx)maxidx=idx;
	      Dbin[idx]++;
	      totcount++;
	    }
      

      sprintf(FILE,"hist_D.dat");
      fp=fopen(FILE,"w");
      for(i=0;i<=maxidx;i++)
        if(Dbin[i]>0)fprintf(fp,"%e %e\n",i*binw,Dbin[i]);
      fclose(fp);
  
      
  

      //sprintf(FILE,"cellconf.dat");
      //fp=fopen(FILE,"w");
      //fprintf(fp,"%d\n",maxcell*maxcell*maxcell);
      //fprintf(fp,"Lattice=\" %d 0.0 0.0 %d %d 0.0 0.0 0.0 %d \" Properties=type:I:1:pos:I:3\n",maxcell,0,maxcell,maxcell);
      //for(i=0;i<maxcell;i++)
      // 	for(j=0;j<maxcell;j++)
      // 	  for(k=0;k<maxcell;k++)
      // 	    fprintf(fp,"%d %d %d %d\n",cell[i][j][k],i,j,k);
      // 
      // 
      //fclose(fp);




      
      
      return 0;
} //end main

      





