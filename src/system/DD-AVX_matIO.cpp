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
      printf("filename is NULL\n");
      abort();
   }
   file = fopen(filename, "r");
   if( file==NULL )
   {
      printf("can not open file %s\n",filename);
      abort();
   }
   if( fgets(buf, 256, file) == NULL )
   {
      printf("blank file %s\n",filename);
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
      printf("Not Matrix Market banner, banner is %s\n",banner);
      abort();
   }
   if( strncmp(fmt, MM_FMT, strlen(MM_FMT))!=0 )
   {
      printf("Not Coodinate format\n");
      abort();
   }
   if( strncmp(dtype, MM_TYPE_REAL, strlen(MM_TYPE_REAL))!=0 )
   {
      printf("Not real\n");
      abort();
   }
   if( strncmp(dstruct, MM_TYPE_GENERAL, strlen(MM_TYPE_GENERAL))!=0 )
   {
      printf("Not general\n");
      abort();
   }
   do
   {
      if( fgets(buf, 1024, file) == NULL )
      {
	 printf("check size error\n");
	 abort();
      }
   }while( buf[0]=='%' );


   //check size
   if( sscanf(buf, "%d %d %d", &col, &row, &nz ) != 3 )
   {
      printf("matrix size in unknown\n");
      abort();
      if(col != row){
	 printf("matrix is not square matrix (row!=col)\n");
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
	 printf("error, format is %d",format);
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
   val = new double[nnz];
   printf("%d %d\n",N,nnz);

   for(i=0;i<nnz;i++)
   {
      if( fgets(buf, 1024, file) == NULL )
      {
	 printf("can't read data, [row col value]\n");
	 printf("bar\n");
	 abort();
      }
      if( sscanf(buf,"%d %d %lf",&idx, &jdx, &value) != 3 )
      {
	 printf("not data, [col=%d,row=%d,val=%f]\n",idx,jdx,value); abort();
      }
      col[i] = idx-1;
      row[i] = jdx-1;
      val[i] = value;
   }
}

void D_Matrix::input_crs(FILE* file){
	char	buf[1024];
	int	i,j;
	int idx, jdx, count=0, jb = 1;
	double value;

	/* read data */
	row = new int[N+1];
	col = new int[nnz];
	val = new double[nnz];
	row[0] = 0;
	printf("format CRS, N = %d, nnz = %d, filesize = %.3f KB\n",N,nnz,((double)nnz*2+N)/1000);

	for(i=0;i<nnz;i++)
	{
	  if( fgets(buf, 1024, file) == NULL )
	    {
	      printf("can't read data, [row col value]\n");
	      printf("hoge\n");
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
	/*
	for(j=0;j<nnz;j++)
	  printf("%d  %d\n",j,col[j]);
	printf("----------------------こる--------------------------------\n");
	for(j=0;j<N+1;j++)
	  printf("%d  %d\n",j,row[j]);
	printf("---------------------ろお--------------------------------------\n");
	for(j=0;j<nnz;j++)
	  printf("%d  %d\n",j,val[j]);
	printf("---------------------------------val----------------------------\n");
	*/
}
