/* Copyright (C) T.Hishinuma, All rights reserved.*/
#include<DD-AVX.hpp>
#include<stdio.h>
#include<string.h>
#include<ctype.h>
#define MM_BANNER		"%%MatrixMarket"
#define MM_MTX			"matrix"
#define MM_VEC			"vector"
#define MM_FMT			"coordinate"
#define MM_TYPE_REAL		"real"
#define MM_TYPE_GENERAL		"general"
#define MM_TYPE_SYMM		"symmetric"

int flg = 0;
///////////// DD_Vector input ////////////////
void D_Matrix::input(const char* filename){
#ifdef ddavx_debug
	 printf("D_Matrix::input()\n");
#endif
   char	buf[256],banner[128];
   FILE*	file;
   char mtx[64], fmt[64], dtype[64], dstruct[64];
   char	*p;
   int	row, col, nz;

   //file open
   if( filename==NULL )
   {
      printf("DD-AVX system: filename is NULL\n");
      abort();
   }
   file = fopen(filename, "r");
   if( file==NULL )
   {
      printf("DD-AVX system: can not open file %s\n",filename);
      abort();
   }
   if( fgets(buf, 256, file) == NULL )
   {
      printf("DD-AVX system: blank file %s\n",filename);
      abort();
   }

   //check banner
   rewind(file);
   fgets(buf, 1024, file);
   sscanf(buf, "%s %s %s %s %s", banner, mtx, fmt, dtype, dstruct);

   for(p=mtx;*p!='\0';p++)     *p = (char)tolower(*p);
   for(p=fmt;*p!='\0';p++)     *p = (char)tolower(*p);
   for(p=dtype;*p!='\0';p++)   *p = (char)tolower(*p);
   for(p=dstruct;*p!='\0';p++) *p = (char)tolower(*p);

   if( strncmp(banner, MM_BANNER, strlen(MM_BANNER)) != 0)
   {
      printf("DD-AVX system: Not Matrix Market banner, banner is %s\n",banner);
      abort();
   }
   if( strncmp(fmt, MM_FMT, strlen(MM_FMT))!=0 )
   {
      printf("DD-AVX system: Not Coodinate format\n");
      abort();
   }
   if( strncmp(dtype, MM_TYPE_REAL, strlen(MM_TYPE_REAL))!=0 )
   {
      printf("DD-AVX system: Not real\n");
      abort();
   }
   if( strncmp(dstruct, MM_TYPE_GENERAL, strlen(MM_TYPE_GENERAL))!=0 )
   {
     if( strncmp(dstruct, MM_TYPE_SYMM, strlen(MM_TYPE_SYMM))==0 )
     {
       flg = 1;
     } else {
       printf("DD-AVX system: Not general or symmetric\n");
       abort();
     }
   }
   do
   {
      if( fgets(buf, 1024, file) == NULL )
      {
	 printf("DD-AVX system: check size error\n");
	 abort();
      }
   }while( buf[0]=='%' );


   //check size
   if( sscanf(buf, "%d %d %d", &col, &row, &nz ) != 3 )
   {
      printf("DD-AVX system: matrix size in unknown\n");
      abort();
      if(col != row){
	 printf("DD-AVX system: matrix is not square matrix (row!=col)\n");
	 abort();
      }
   }
   N = row;
   nnz = nz;

   //input MM data as selected matrix storage format
   switch(format){
      case 0:
   	input_crs(file);
	 break;
      case 1:
	input_crs(file);
	 break;
      case 2:
	input_crs(file);
	 break;
      case 3:
	input_crs(file);
	 break;
      default:
	 printf("DD-AVX system: error, format is %d",format);
   }

   rewind(file);
   fclose(file);
}

void D_Matrix::input_coo(FILE *file){
#ifdef ddavx_debug
	 printf("    D_Matrix::input_coo()\n");
#endif
   char	buf[1024];
   int	i;
   int idx, jdx;
   double value;
   //rewind(file);

   /* read data */
   row = new int[nnz];
   col = new int[nnz];
   val = new double[nnz];//double
   printf("%d %d\n",N,nnz);

   for(i=0;i<nnz;i++)
   {
     if( fgets(buf, 1024, file) == NULL )
      {
	 printf("can't read data, [row col value]\n");
	 printf("bar\n");
	 abort();

	 printf("DD-AVX system: cant read data, [row col value]\n"); abort();
      }
     if( sscanf(buf,"%d %d %lf",&idx, &jdx, &value) != 3 )
      {
	 printf("DD-AVX system: not data, [col=%d,row=%d,val=%f]\n",idx,jdx,value); abort();
      }
      col[i] = idx-1;
      row[i] = jdx-1;
      val[i] = value;
   }
}

void D_Matrix::input_crs(FILE* file){
	char	buf[1024];
	int	i,j;
	int idx, jdx, count = 0, jb = 1;
	double value;
	int nnz2 = nnz;

	if(!flg)
	  {
	    /* read data */
	    row = new int[N+1];
	    col = new int[nnz];
	    val = new double[nnz];
	    row[0] = 0;
	    //	printf("DD-AVX system: format CRS, N = %d, nnz = %d, filesize = %.3f KB\n",N,nnz,((double)nnz*2+N)/1000);
	  for(i=0;i<nnz;i++)
	    {
	      if( fgets(buf, 1024, file) == NULL )
		{
		  printf("can't read data, [row col value]\n");
		  abort();
		}
	  
	      if( sscanf(buf,"%d %d %lf",&idx, &jdx, &value) != 3 )
		{
		  printf("not data, [col=%d,row=%d,val=%f]\n",idx,jdx,value);
		  abort();
		}
	      //printf("%d %d %e %d\n",idx,jdx,value, jb);
	  
	      if(jb != jdx){
		row[jdx-1] = i-1;
		//printf("%d,%d\n",jdx-1,row[jdx-1]);
	      }
	      jb = jdx;
	      count++;
	  
	      col[i] = idx-1;
	      val[i] = value;
		
	    }
	  row[N] = nnz;
	  for(j=1;j<N;j++)
	    row[j] += 1;
	}else{
	    /* read data */
	    row = new int[N+1];
	    col = new int[nnz];
	    val = new double[nnz];
	    row[0] = 0;
	  for(i=0;i<nnz;i++)
	    {
	      
	      if( fgets(buf, 1024, file) == NULL )
		{
		  printf("can't read data, [row col value]\n");
		  abort();
		}
	      
	      if( sscanf(buf,"%d %d %lf",&idx, &jdx, &value) != 3 )
		{
		  printf("not data, [col=%d,row=%d,val=%f]\n",idx,jdx,value);
		  abort();
		}

	      if(idx != jdx){
		nnz2 += 1;
	      }
	      
	      if( jb != jdx ){
		row[jdx-1] = i-1;
	      }
	      
	      jb = jdx;
	      count++;
	      col[i] = idx-1;
	      val[i] = value;
	    }
	  row[N] = nnz2;
	  for(j=1;j<N;j++)
	    row[j]++;
	  }
}

void D_Matrix::output_plane(const char *filename)
{
#ifdef ddavx_debug
	 printf("D_Matrix::output_plane()\n");
#endif
   char	buf[1024];
   int	i;
   FILE *file;

   file = fopen(filename, "w");
   if( file==NULL )
   {
      printf("can not open file %s",filename);
      abort();
   }
   if(N==0){
      printf("vector size N is 0\n");
      abort();
   }
   int count = 0;
   fprintf(file,"%%%MatrixMarket matrix coordinate real general\n");
   fprintf(file,"%d %d %d\n",N,N,nnz);
   for(int i=0;i<nnz;i++){
     if(row[count+1] <= i)
       count++;
     fprintf(file,"%d %d %1.15f\n",col[i]+1,count+1,val[i]);
   }
   /*
   for(i=0;i<N;i++)
   {
      fprintf(file,"%20.20e\n",hi[i]);
   }
   */
   fclose(file);
}
