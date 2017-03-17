/* Copyright (C) T.Hishinuma, All rights reserved.*/
#include<DD-AVX.hpp>
#include<stdio.h>

//1//////////////////////////////////////////////////////////
void DD_AVX_axpyz(D_Scalar alpha, D_Vector vx, D_Vector vy, D_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(D Scalar, D Vector, D Vector, D_Vector)\n");
#endif
   DD_AVX_axpyz_D(alpha,vx,vy,vz);
}
//2//////////////////////////////////////////////////////////

void DD_AVX_axpyz(D_Scalar alpha, DD_Vector vx, D_Vector vy, D_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(D Scalar, DD Vector, D Vector, D_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}
//3//////////////////////////////////////////////////////////

void DD_AVX_axpyz(D_Scalar alpha, D_Vector vx, DD_Vector vy, D_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(D Scalar, D Vector, DD Vector, D_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}
//4//////////////////////////////////////////////////////////

void DD_AVX_axpyz(D_Scalar alpha, DD_Vector vx, DD_Vector vy, D_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(D Scalar, DD Vector, DD Vector, D_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}
//5//////////////////////////////////////////////////////////

void DD_AVX_axpyz(DD_Scalar alpha, D_Vector vx, D_Vector vy, D_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(DD Scalar, D Vector, D Vector, D_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}
//6//////////////////////////////////////////////////////////

void DD_AVX_axpyz(DD_Scalar alpha, DD_Vector vx, D_Vector vy, D_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(DD Scalar, DD Vector, D Vector, D_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}
//7//////////////////////////////////////////////////////////

void DD_AVX_axpyz(DD_Scalar alpha, D_Vector vx, DD_Vector vy, D_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(DD Scalar, D Vector, DD Vector, D_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}
//8//////////////////////////////////////////////////////////

void DD_AVX_axpyz(DD_Scalar alpha, DD_Vector vx, DD_Vector vy, D_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(DD Scalar, DD Vector, DD Vector, D_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}

//9//////////////////////////////////////////////////////////
//9//////////////////////////////////////////////////////////
//9//////////////////////////////////////////////////////////

void DD_AVX_axpyz(D_Scalar alpha, D_Vector vx, D_Vector vy, DD_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(D Scalar, D Vector, D Vector, DD_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}
//10//////////////////////////////////////////////////////////

void DD_AVX_axpyz(D_Scalar alpha, DD_Vector vx, D_Vector vy, DD_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(D Scalar, DD Vector, D Vector, DD_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}
//11//////////////////////////////////////////////////////////

void DD_AVX_axpyz(D_Scalar alpha, D_Vector vx, DD_Vector vy, DD_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(D Scalar, D Vector, DD Vector, DD_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}
//12//////////////////////////////////////////////////////////

void DD_AVX_axpyz(D_Scalar alpha, DD_Vector vx, DD_Vector vy, DD_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(D Scalar, DD Vector, DD Vector, DD_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}
//13//////////////////////////////////////////////////////////

void DD_AVX_axpyz(DD_Scalar alpha, D_Vector vx, D_Vector vy, DD_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(DD Scalar, D Vector, D Vector, DD_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}
//14//////////////////////////////////////////////////////////

void DD_AVX_axpyz(DD_Scalar alpha, DD_Vector vx, D_Vector vy, DD_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(DD Scalar, DD Vector, D Vector, DD_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}
//15//////////////////////////////////////////////////////////

void DD_AVX_axpyz(DD_Scalar alpha, D_Vector vx, DD_Vector vy, DD_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(DD Scalar, D Vector, DD Vector, DD_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}
//16//////////////////////////////////////////////////////////

void DD_AVX_axpyz(DD_Scalar alpha, DD_Vector vx, DD_Vector vy, DD_Vector vz){
#ifdef ddavx_debug
   printf("axpyz(DD Scalar, DD Vector, DD Vector, DD_Vector)\n");
#endif
   DD_AVX_axpyz_DD(DD(alpha),DD(vx),DD(vy),DD(vz));
}

////////////////////////////////////////////////////////////

void DD_AVX_axpyz_D(D_Scalar alpha, D_Vector vx, D_Vector vy, D_Vector vz)
{
#ifdef ddavx_debug
   printf("    compute axpyz in double\n");
#endif

   int i,n,is,ie,nprocs,my_rank;
   DD_AVX_DECLAR;
   double *x,*y,*z;

   n   = vx.N;
   x   = vx.hi;
   y   = vy.hi;
   z   = vz.hi;

#pragma omp parallel for
   for(i=0; i<n; i++)
   {
      z[i] = alpha.hi * x[i] + y[i];
   }
}

void DD_AVX_axpyz_DD(DD_Scalar alpha, DD_Vector vx, DD_Vector vy, DD_Vector vz)
{
#ifdef ddavx_debug
   printf("    compute axpyz in DD\n");
#endif

   int i,n,is,ie,nprocs,my_rank;
   DD_AVX_DECLAR;
   double *x,*xl,*y,*yl,*z,*zl;

   n   = vx.N;
   x   = vx.hi;
   y   = vy.hi;
   z   = vz.hi;

   xl  = vx.lo;
   yl  = vy.lo;
   zl  = vz.lo;

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
	 DD_AVX_FMAN_AVX(z[i],zl[i],y[i],yl[i],alpha_hi,alpha_lo,x[i],xl[i]);
      }
      for(;i<ie;i++)
      {
	 DD_AVX_FMA_SSE2(z[i],zl[i],y[i],yl[i],alpha.hi,alpha.lo,x[i],xl[i]);
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
	 DD_AVX_FMA2_SSE2(z[i],zl[i],y[i],yl[i],aa.hi[0],aa.lo[0],x[i],xl[i]);
      }
      for(;i<ie;i++)
      {
	 DD_AVX_FMA_SSE2(z[i],zl[i],y[i],yl[i],alpha.hi[0],alpha.lo[0],x[i],xl[i]);
      }
}
#else
//   DD_AVX_DECLARS;
#pragma omp parallel for private(i,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
   for(i=0; i<n; i++)
   {
      DD_AVX_FMA(z[i],zl[i],y[i],yl[i],alpha.hi,alpha.lo,x[i],xl[i]);
   }
#endif
}
