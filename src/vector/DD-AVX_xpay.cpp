/* Copyright (C) T.Hishinuma, All rights reserved.*/
#include<DD-AVX.hpp>
#include<stdio.h>
////////////////////////////////////////////////////////////
extern void DD_AVX_xpay(D_Vector vx, D_Scalar alpha, D_Vector vy){
#ifdef ddavx_debug
   printf("xpay(D Vector, D Scalar, D Vector)\n");
#endif
   DD_AVX_xpay_D(vx,alpha,vy);
}
////////////////////////////////////////////////////////////

extern void DD_AVX_xpay(DD_Vector vx, D_Scalar alpha, D_Vector vy){
#ifdef ddavx_debug
   printf("xpay(DD Vector, D Scalar, D Vector)\n");
#endif
   DD_AVX_xpay_DD(DD(vx),DD(alpha),DD(vy));
}
////////////////////////////////////////////////////////////

extern void DD_AVX_xpay(D_Vector vx, D_Scalar alpha, DD_Vector vy){
#ifdef ddavx_debug
   printf("xpay(D Vector, D Scalar, DD Vector)\n");
#endif
   DD_AVX_xpay_DD(DD(vx),DD(alpha),DD(vy));
}
////////////////////////////////////////////////////////////

extern void DD_AVX_xpay(DD_Vector vx, D_Scalar alpha, DD_Vector vy){
#ifdef ddavx_debug
   printf("xpay(DD Vector, D Scalar, DD Vector)\n");
#endif
   DD_AVX_xpay_DD(DD(vx),DD(alpha),DD(vy));
}
////////////////////////////////////////////////////////////

extern void DD_AVX_xpay(D_Vector vx, DD_Scalar alpha, D_Vector vy){
#ifdef ddavx_debug
   printf("xpay(D Vector, DD Scalar, D Vector)\n");
#endif
   DD_AVX_xpay_DD(DD(vx),DD(alpha),DD(vy));
}
////////////////////////////////////////////////////////////

extern void DD_AVX_xpay(DD_Vector vx, DD_Scalar alpha, D_Vector vy){
#ifdef ddavx_debug
   printf("xpay(DD Vector, DD Scalar, D Vector)\n");
#endif
   DD_AVX_xpay_DD(DD(vx),DD(alpha),DD(vy));
}
////////////////////////////////////////////////////////////

extern void DD_AVX_xpay(D_Vector vx, DD_Scalar alpha, DD_Vector vy){
#ifdef ddavx_debug
   printf("xpay(D Vector, DD Scalar, DD Vector)\n");
#endif
   DD_AVX_xpay_DD(DD(vx),DD(alpha),DD(vy));
}
////////////////////////////////////////////////////////////

extern void DD_AVX_xpay(DD_Vector vx, DD_Scalar alpha, DD_Vector vy){
#ifdef ddavx_debug
   printf("xpay(DD Vector, DD Scalar, DD Vector)\n");
#endif
   DD_AVX_xpay_DD(DD(vx),DD(alpha),DD(vy));
}
////////////////////////////////////////////////////////////

void DD_AVX_xpay_D(D_Vector vx, D_Scalar alpha, D_Vector vy)
{
#ifdef ddavx_debug
   printf("    compute xpay in double\n");
#endif

   int i,n,is,ie,nprocs,my_rank;
   DD_AVX_DECLAR;
   double *x,*y;

   n   = vx.N;
   x   = vx.hi;
   y   = vy.hi;

#pragma omp parallel for
   for(i=0; i<n; i++)
   {
      y[i] = alpha.hi * y[i] + x[i];
   }
}

void DD_AVX_xpay_DD(DD_Vector vx, DD_Scalar alpha, DD_Vector vy)
{
#ifdef ddavx_debug
   printf("    compute axpy in DD\n");
#endif

   int i,n,is,ie,nprocs,my_rank;
   DD_AVX_DECLAR;
   double *x,*xl,*y,*yl;

   n   = vx.N;
   x   = vx.hi;
   y   = vy.hi;
   xl  = vx.lo;
   yl  = vy.lo;

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
	 DD_AVX_FMAN_AVX(y[i],yl[i],x[i],xl[i],alpha_hi,alpha_lo,y[i],yl[i]);
      }
      for(;i<ie;i++)
      {
	 DD_AVX_FMA_SSE2(y[i],yl[i],x[i],xl[i],alpha.hi,alpha.lo,y[i],yl[i]);
      }
}

#elif USE_SSE2==1
   nprocs = omp_get_max_threads();
   DD_Scalar bx,aa,ax;
   aa.hi[0] = aa.hi[1] = alpha.hi[0];
   aa.lo[0] = aa.lo[1] = alpha.lo[0];
#pragma omp parallel private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh,is,ie,my_rank)
{
      my_rank = omp_get_thread_num();
      DD_AVX_GET_ISIE(my_rank,nprocs,n,is,ie);
      for(i=is;i<ie-1;i+=2)
      {
	 DD_AVX_FMA2_SSE2(y[i],yl[i],x[i],xl[i],aa.hi[0],aa.lo[0],y[i],yl[i]);
      }
      for(;i<ie;i++)
      {
	 DD_AVX_FMA_SSE2(y[i],yl[i],x[i],xl[i],alpha.hi[0],alpha.lo[0],y[i],yl[i]);
      }
}
#else
//   DD_AVX_DECLARS;
#pragma omp parallel for private(i,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
   for(i=0; i<n; i++)
   {
      DD_AVX_FMA(y[i],yl[i],x[i],xl[i],alpha.hi,alpha.lo,y[i],yl[i]);
   }
#endif
}
