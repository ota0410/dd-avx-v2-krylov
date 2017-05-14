/* Copyright (C) T.Hishinuma, All rights reserved.*/
#include<DD-AVX.hpp>
#include<stdio.h>

//1//////////////////////////////////////////////////////////
void DD_AVX_nrm2(D_Vector vx, D_Scalar* val){
#ifdef ddavx_debug
   printf("D scalar = nrm2(D vector)\n");
#endif
   D_Scalar tval;
   DD_AVX_nrm2_D(vx, &tval);
   val->hi = tval.hi;
}
//2//////////////////////////////////////////////////////////

void DD_AVX_nrm2(DD_Vector vx, D_Scalar* val){
#ifdef ddavx_debug
   printf("D scalar = nrm2(DD vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_nrm2_DD(DD(vx), &tval);
   val->hi = tval.hi;
}
//3//////////////////////////////////////////////////////////

void DD_AVX_nrm2(D_Vector vx, DD_Scalar* val){
#ifdef ddavx_debug
   printf("DD scalar = nrm2(D vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_nrm2_DD(DD(vx), &tval);
   val->hi = tval.hi;
}
//4//////////////////////////////////////////////////////////

void DD_AVX_nrm2(DD_Vector vx, DD_Scalar* val){
#ifdef ddavx_debug
   printf("DD scalar = nrm2(DD vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_nrm2_DD(DD(vx), &tval);
   val->hi = tval.hi;
}
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

//1//////////////////////////////////////////////////////////
D_Scalar D_Scalar::nrm2(D_Vector vx){
#ifdef ddavx_debug
   printf("D scalar = nrm2(D vector)\n");
#endif
   D_Scalar tval;
   DD_AVX_nrm2_D(vx, &tval);
   this->hi = tval.hi;
   return *this;
}
//2//////////////////////////////////////////////////////////

D_Scalar D_Scalar::nrm2(DD_Vector vx){
#ifdef ddavx_debug
   printf("D scalar = nrm2(DD vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_nrm2_DD(DD(vx), &tval);
   this->hi = tval.hi;
   return *this;
}
//3//////////////////////////////////////////////////////////

DD_Scalar DD_Scalar::nrm2(D_Vector vx){
#ifdef ddavx_debug
   printf("DD scalar = nrm2(D vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_nrm2_DD(DD(vx), &tval);
   this->hi = tval.hi;
   this->lo = tval.lo;
   return *this;
}
//4//////////////////////////////////////////////////////////

DD_Scalar DD_Scalar::nrm2(DD_Vector vx){
#ifdef ddavx_debug
   printf("DD scalar = nrm2(DD vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_nrm2_DD(DD(vx), &tval);
   this->hi = tval.hi;
   this->lo = tval.lo;
   return *this;
}
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
void DD_AVX_nrm2_D(D_Vector vx, D_Scalar *val){
#ifdef ddavx_debug
   printf("    compute nrm2 in double\n");
#endif

   int i,n;
   double *x, dot;

   n = vx.N;
   x = vx.hi;
   dot = val->hi;

#pragma omp parallel for reduction(+: dot)
      for(i=0; i<n; i++)
      {
	 dot += x[i] * x[i];
      }
   val->hi = sqrt(dot);
}


void DD_AVX_nrm2_DD(DD_Vector vx, DD_Scalar* val){
#ifdef ddavx_debug
   printf("   compute nrm2 in DD\n");
#endif
   int i,j,n;
   int is,ie,nprocs,my_rank;
   double *gt,*gtavx;
   double *x,*xl;
   DD_AVX_DECLAR;
   nprocs = omp_get_max_threads();

#ifdef USE_AVX
   DD_Scalar dotm;
   __m256d dotm_v_hi0, dotm_v_lo0;
   gtavx = (double*)malloc(nprocs*8*sizeof(double));
#elif defined USE_SSE2
   __m128d dotm_v_hi0, dotm_v_lo0;
   gt = (double*)malloc(nprocs*4*sizeof(double));
#endif

   n  = vx.N;

   x  = vx.hi;
   xl = vx.lo;

#pragma omp parallel private(i,j,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh,t3,is,ie,my_rank)
   {
      my_rank = omp_get_thread_num();
      DD_AVX_GET_ISIE(my_rank,nprocs,n,is,ie);
#if defined(USE_AVX)
      dotm_v_hi0 = DD_AVX_AVX_FUNC(setzero_pd)();
      dotm_v_lo0 = DD_AVX_AVX_FUNC(setzero_pd)();
      double a[8]={0,0,0,0,0,0,0,0};

      for(i=is;i<ie-(DD_AVX_AVX_SIZE-1);i+=DD_AVX_AVX_SIZE)
      {
	 DD_AVX_FMAN_AVX(a[0],a[4],a[0],a[4],x[i],xl[i],x[i],xl[i]);
      }

      for(;i<ie;i++)
      {
	 DD_AVX_FMA_SSE2(a[0],a[4],a[0],a[4],x[i],xl[i],x[i],xl[i]);
      }

	  DD_AVX_HADDALL_AVX(a[0],a[4],a[0],a[4]);
      for(j=0;j<8;j++){
	 gtavx[my_rank*8+j]=a[j];
      }
#elif defined(USE_FMA2_SSE2)
      gt[my_rank*2] = gt[my_rank*2+1] = 0.0;
      gt[my_rank*2+2] = gt[my_rank*2+3] = 0.0;
      for(i=is;i<ie-1;i+=2)
      {
	 DD_AVX_FMA2_SSE2(gt[my_rank*2],gt[my_rank*2+2],gt[my_rank*2],gt[my_rank*2+2],\
	       x[i],xl[i],x[i],xl[i]);
      }
      DD_AVX_ADD_SSE2(gt[my_rank*2],gt[my_rank*2+1],gt[my_rank*2],gt[my_rank*2+1],gt[my_rank*2+2],gt[my_rank*2+3]);
      for(;i<ie;i++)
      {
	 DD_AVX_FMA_SSE2(gt[my_rank*2],gt[my_rank*2+1],gt[my_rank*2],gt[my_rank*2+1],\
	       x[i],xl[i],x[i],xl[i]);
      }
#else
      gt[my_rank*2] = gt[my_rank*2+1] = 0.0;
      for(i=is;i<ie;i++)
      {
	 DD_AVX_FMA(gt[my_rank*2],gt[my_rank*2+1],gt[my_rank*2],gt[my_rank*2+1],\
	       x[i],xl[i],x[i],xl[i]);
      }
#endif
   }

#ifdef USE_AVX
   dotm_v_hi0 = DD_AVX_AVX_FUNC(setzero_pd)();
   dotm_v_lo0 = DD_AVX_AVX_FUNC(setzero_pd)();
#endif

   for(i=0;i<nprocs;i++)
   {
#ifdef USE_AVX
      DD_AVX_ADD_SSE2(dotm.hi,dotm.lo,dotm.hi,dotm.lo,gtavx[i*8],gtavx[i*8+4]);
#else
      DD_AVX_ADD_SSE2(dotm.hi,dotm.lo,dotm.hi,dotm.lo,gt[i*2],gt[i*2+1]);
#endif	
   }
#ifndef USE_SSE2
   DD_AVX_SQRT(val->hi,val->lo,dotm.hi,dotm.lo);
#else
   DD_AVX_SQRT_SSE2(val->hi,val->lo,dotm.hi,dotm.lo);
#endif
}
