/* Copyright (C) T.Hishinuma, All rights reserved.*/
#include<DD-AVX.hpp>
#include<stdio.h>
#include<string.h>
#include <immintrin.h>

void DD_AVX_SpMV_CRS_D(D_Matrix A, D_Vector vx, D_Vector vy)
{
#ifdef ddavx_debug
   printf("\tcompute CRS SpMV in double\n");
#endif
   int i,j,n,is,ie,nprocs,my_rank;
   DD_AVX_DECLAR;
   double *x,*y;
   double k;

   n   = vx.N;
   x   = vx.hi;
   y   = vy.hi;
#pragma omp parallel for
   for(i=0;i<n;i++){
      y[i] = 0;
   }

#pragma omp parallel for private(j,k)
   for(i=0;i<n;i++){
     k = y[i];
     for(j=A.row[i];j<A.row[i+1];j++){
       k += A.val[j] * x[A.col[j]];
     }
     y[i] = k;
   }
}

void DD_AVX_SpMV_CRS_DD(D_Matrix A, DD_Vector vx, DD_Vector vy)
{
#ifdef ddavx_debug
   printf("\tcompute CRS SpMV in DD\n");
#endif
   int	i,j,n;
   int	is,ie;
   int	j0,j1;
#if defined(USE_AVX)
   int j2,j3,ij;
   DD_AVX_AVX_TYPE	tt_hi, tt_lo;
   double		s_hi, s_lo;
   DD_AVX_AVX_TYPE	tmph,tmpl;
#endif
   
   int				*jj0;
   double		*vv0;
   double		*x,*y,*xl,*yl;
   DD_AVX_DECLAR;
   n     = A.N;
   x     = vx.hi;
   y     = vy.hi;
   xl    = vx.lo;
   yl    = vy.lo;
   //mp3_matvec
   jj0 = A.col;
   vv0 = A.val;

#if defined(USE_AVX)
#pragma omp parallel for schedule(guided) private(j,is,ie,j0,j1,j2,j3,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh,t3,tt_hi,tt_lo,s_hi,s_lo)
#elif defined(USE_SSE2)
#pragma omp parallel for schedule(guided) private(j,is,ie,j0,j1,tt,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
#else
#pragma omp parallel for schedule(guided) private(j,is,ie,j0,j1,tt,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
#endif

   for(i=0;i<n+1;i++)
   {
#if defined(USE_AVX)
      DD_AVX_AVX_TYPE tt_hi = DD_AVX_AVX_FUNC(setzero_pd)();
      DD_AVX_AVX_TYPE tt_lo = DD_AVX_AVX_FUNC(setzero_pd)();

      is = A.ptr[i];
      ie = A.ptr[i+1];
      for(j=is;j<ie-(DD_AVX_AVX_SIZE-1);j+=DD_AVX_AVX_SIZE)
      {
	      
	 j0 = jj0[j+0]; j1 = jj0[j+1];
	 j2 = jj0[j+2]; j3 = jj0[j+3];
	 DD_AVX_FMAD4_AVX_LDSD(tt_hi,tt_lo,tt_hi,tt_lo,x[j0],xl[j0],x[j1],xl[j1],x[j2],xl[j2],x[j3],xl[j3],vv0[j]);
      }

#define padding 1
#if padding == 1
      if(j==ie-3)
      {
	 j0 = jj0[j+0]; j1 = jj0[j+1];
	 j2 = jj0[j+2]; 
	 DD_AVX_FMAD4_AVX_LDSD(tt_hi,tt_lo,tt_hi,tt_lo,x[j0],xl[j0],x[j1],xl[j1],x[j2],xl[j2],0,0,vv0[j]);
	 j+=3;
      }
      else if(j==ie-2)
      {
	 j0 = jj0[j+0]; j1 = jj0[j+1];
	 DD_AVX_FMAD4_AVX_LDSD(tt_hi,tt_lo,tt_hi,tt_lo,x[j0],xl[j0],x[j1],xl[j1],0,0,0,0,vv0[j]);
	 j+=2;
      }
      else if(j==ie-1)
      {
	 j0 = jj0[j+0];
	 DD_AVX_FMAD4_AVX_LDSD(tt_hi,tt_lo,tt_hi,tt_lo,x[j0],xl[j0],0,0,0,0,0,0,vv0[j]);
	 j+=1;
      }

      DD_AVX_HADDALL_AVX(y[i],yl[i],tt_hi,tt_lo);
#endif

#elif defined(USE_SSE2)
      tt.hi[0] = tt.hi[1] = tt.lo[0] = tt.lo[1] = 0.0;

      is = A->ptr[i];
      ie = A->ptr[i+1];

      for(j=is;j<ie-1;j+=2)
      {
#if 1
	 j0 = jj0[j+0];
	 j1 = jj0[j+1];
#ifdef USE_SSE2
	 DD_AVX_QUAD_FMAD2_SSE2_LDSD(tt.hi[0],tt.lo[0],tt.hi[0],tt.lo[0],x[j0],xl[j0],x[j1],xl[j1],vv0[j]);
#endif
#else
	 j0 = jj0[j];
	 DD_AVX_QUAD_FMAD2_SSE2(tt.hi[0],tt.lo[0],tt.hi[0],tt.lo[0],x[j0],xl[j0],vv0[j]);
#endif
      }
#ifdef USE_SSE2
      DD_AVX_QUAD_ADD_SSE2(y[i],yl[i],tt.hi[0],tt.lo[0],tt.hi[1],tt.lo[1]);
#endif
      for(;j<ie;j++)
      {
	 j0 = jj0[j+0];
#ifdef USE_SSE2
	 DD_AVX_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],vv0[j]);
#endif
      }
      #else
      #error no implementation
#endif
   }

}


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/

void DD_AVX_TSpMV_CRS_D(D_Matrix A, D_Vector vx, D_Vector vy)
{
#ifdef ddavx_debug
   printf("\tcompute CRS SpMV in double\n");
#endif
   int i,j,n,is,ie,nprocs,my_rank;
   DD_AVX_DECLAR;
   double *x,*y;
   double k;

   n   = vx.N;
   x   = vx.hi;
   y   = vy.hi;
#pragma omp parallel for
   for(i=0;i<n;i++){
      y[i] = 0;
   }

   for(i=0;i<n;i++){
     k = x[i];
     for(j=A.row[i];j<A.row[i+1];j++){
       y[A.col[j]] += A.val[j] * k;
     }
   }
}

void DD_AVX_TSpMV_CRS_DD(D_Matrix A, DD_Vector vx, DD_Vector vy)
{
#ifdef ddavx_debug
   printf("\tcompute CRS SpMV in DD\n");
#endif
   int	i,j,n;
   int	is,ie;
   int	j0,j1;
#if defined(USE_AVX)
   int j2,j3,ij;
   DD_AVX_AVX_TYPE	tt_hi, tt_lo;
   double		s_hi, s_lo;
   DD_AVX_AVX_TYPE	tmph,tmpl;
#endif
   
   int				*jj0;
   double		*vv0;
   double		*x,*y,*xl,*yl;
   DD_AVX_DECLAR;
   n     = A.N;
   x     = vx.hi;
   y     = vy.hi;
   xl    = vx.lo;
   yl    = vy.lo;
   //mp3_matvec
   jj0 = A.col;
   vv0 = A.val;
   
#if defined(USE_AVX)
#pragma omp parallel for schedule(guided) private(j,is,ie,j0,j1,j2,j3,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh,t3,tt_hi,tt_lo,s_hi,s_lo) 
#elif defined(USE_SSE2)
#pragma omp parallel for schedule(guided) private(j,is,ie,j0,j1,tt,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
#else
#pragma omp parallel for schedule(guided)  private(j,is,ie,j0,j1,tt,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
#endif
  
   for(i=0;i<n+1;i++)
   {
#if defined(USE_AVX)
      DD_AVX_AVX_TYPE tt_hi = DD_AVX_AVX_FUNC(setzero_pd)();
      DD_AVX_AVX_TYPE tt_lo = DD_AVX_AVX_FUNC(setzero_pd)();

      is = A.row[i];
      ie = A.row[i+1];
      for(j=is;j<ie-(DD_AVX_AVX_SIZE-1);j+=DD_AVX_AVX_SIZE)
      {
	      
	 j0 = jj0[j+0]; j1 = jj0[j+1];
	 j2 = jj0[j+2]; j3 = jj0[j+3];
	 DD_AVX_FMAD4_AVX_LDSD(tt_hi,tt_lo,tt_hi,tt_lo,x[j0],xl[j0],x[j1],xl[j1],x[j2],xl[j2],x[j3],xl[j3],vv0[j]);
      }

#define padding 1
#if padding == 1
      if(j==ie-3)
      {
	 j0 = jj0[j+0]; j1 = jj0[j+1];
	 j2 = jj0[j+2]; 
	 DD_AVX_FMAD4_AVX_LDSD(tt_hi,tt_lo,tt_hi,tt_lo,x[j0],xl[j0],x[j1],xl[j1],x[j2],xl[j2],0,0,vv0[j]);
	 j+=3;
      }
      else if(j==ie-2)
      {
	 j0 = jj0[j+0]; j1 = jj0[j+1];
	 DD_AVX_FMAD4_AVX_LDSD(tt_hi,tt_lo,tt_hi,tt_lo,x[j0],xl[j0],x[j1],xl[j1],0,0,0,0,vv0[j]);
	 j+=2;
      }
      else if(j==ie-1)
      {
	 j0 = jj0[j+0];
	 DD_AVX_FMAD4_AVX_LDSD(tt_hi,tt_lo,tt_hi,tt_lo,x[j0],xl[j0],0,0,0,0,0,0,vv0[j]);
	 j+=1;
      }

      DD_AVX_HADDALL_AVX(y[i],yl[i],tt_hi,tt_lo);
#endif

#elif defined(USE_SSE2)
      tt.hi[0] = tt.hi[1] = tt.lo[0] = tt.lo[1] = 0.0;

      is = A->ptr[i];
      ie = A->ptr[i+1];

      for(j=is;j<ie-1;j+=2)
      {
#if 1
	 j0 = jj0[j+0];
	 j1 = jj0[j+1];
#ifdef USE_SSE2
	 DD_AVX_QUAD_FMAD2_SSE2_LDSD(tt.hi[0],tt.lo[0],tt.hi[0],tt.lo[0],x[j0],xl[j0],x[j1],xl[j1],vv0[j]);
#endif
#else
	 j0 = jj0[j];
	 DD_AVX_QUAD_FMAD2_SSE2(tt.hi[0],tt.lo[0],tt.hi[0],tt.lo[0],x[j0],xl[j0],vv0[j]);
#endif
      }
#ifdef USE_SSE2
      DD_AVX_QUAD_ADD_SSE2(y[i],yl[i],tt.hi[0],tt.lo[0],tt.hi[1],tt.lo[1]);
#endif
      for(;j<ie;j++)
      {
	 j0 = jj0[j+0];
#ifdef USE_SSE2
	 DD_AVX_QUAD_FMAD_SSE2(y[i],yl[i],y[i],yl[i],x[j0],xl[j0],vv0[j]);
#endif
      }
#else
#error no implementation
#endif
   }

}
