/* Copyright (C) T.Hishinuma, All rights reserved.*/
#include<DD-AVX.hpp>
#include<stdio.h>

////////////////////////////////////////////////////////////
void DD_AVX_scale(D_Scalar alpha, D_Vector vx){
#ifdef ddavx_debug
   printf("scale(D Scalar, D Vector)\n");
#endif
   DD_AVX_scale_D(alpha,vx);
}
////////////////////////////////////////////////////////////

void DD_AVX_scale(DD_Scalar alpha, D_Vector vx){
#ifdef ddavx_debug
   printf("scale(DD Scalar, D Vector)\n");
#endif
   DD_AVX_scale_DD(DD(alpha),DD(vx));
}
////////////////////////////////////////////////////////////

void DD_AVX_scale(D_Scalar alpha, DD_Vector vx){
#ifdef ddavx_debug
   printf("scale(D Scalar, DD Vector)\n");
#endif
   DD_AVX_scale_DD(DD(alpha),DD(vx));
}
////////////////////////////////////////////////////////////

void DD_AVX_scale(DD_Scalar alpha, DD_Vector vx){
#ifdef ddavx_debug
   printf("scale(DD Scalar, DD Vector)\n");
#endif
   DD_AVX_scale_DD(DD(alpha),DD(vx));
}
////////////////////////////////////////////////////////////

void DD_AVX_scale_D(D_Scalar alpha, D_Vector vx)
{
#ifdef ddavx_debug
   printf("    compute scale in double\n");
#endif

   int i,n,is,ie,nprocs,my_rank;
   DD_AVX_DECLAR;
   double *x;

   n   = vx.N;
   x   = vx.hi;

#pragma omp parallel for
   for(i=0; i<n; i++)
   {
      x[i] = alpha.hi * x[i];
   }
}

void DD_AVX_scale_DD(DD_Scalar alpha, DD_Vector vx)
{
#ifdef ddavx_debug
   printf("    compute scale in DD\n");
#endif

   int i,n,is,ie,nprocs,my_rank;
   DD_AVX_DECLAR;
   double *x,*xl;

   n   = vx.N;
   x   = vx.hi;
   xl  = vx.lo;

#if USE_AVX==1
   nprocs = omp_get_max_threads();
   DD_AVX_AVX_TYPE alpha_hi, alpha_lo;
   alpha_hi = DD_AVX_AVX_FUNC(broadcast_sd)(&alpha.hi);
   alpha_lo = DD_AVX_AVX_FUNC(broadcast_sd)(&alpha.lo);
#pragma omp parallel private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh,t3,is,ie,my_rank)
{
      my_rank = omp_get_thread_num();

      DD_AVX_GET_ISIE(my_rank,nprocs,n,is,ie);
      for(i=is;i<ie-(DD_AVX_AVX_SIZE-1);i+=DD_AVX_AVX_SIZE)
      {
	 DD_AVX_MULN_AVX(x[i],xl[i],x[i],xl[i],alpha_hi,alpha_lo);
      }
      for(;i<ie;i++)
      {
	 DD_AVX_MUL_SSE2(x[i],xl[i],x[i],xl[i],alpha.hi,alpha.lo);
      }
}

#elif USE_SSE2==1
   nprocs = omp_get_max_threads();
   DD_Scalar bx,aa,ax;
   aa.hi[0] = aa.hi[1] = alpha.hi;
   aa.lo[0] = aa.lo[1] = alpha.lo;
#pragma omp parallel private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh,is,ie,my_rank)
{
      my_rank = omp_get_thread_num();
      DD_AVX_GET_ISIE(my_rank,nprocs,n,is,ie);
      for(i=is;i<ie-1;i+=2)
      {
	 DD_AVX_MUL2_SSE2(x[i],xl[i],x[i],xl[i],aa.hi[0],aa.lo[0]);
      }
      for(;i<ie;i++)
      {
	 DD_AVX_MUL_SSE2(x[i],xl[i],x[i],xl[i],aa.hi[0],aa.lo[0]);
      }
}
#else
//   DD_AVX_DECLARS;
#pragma omp parallel for private(i,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
   for(i=0; i<n; i++)
   {
      DD_AVX_MUL(x[i],xl[i],x[i],xl[i],alpha.hi,alpha.lo);
   }
#endif
}
