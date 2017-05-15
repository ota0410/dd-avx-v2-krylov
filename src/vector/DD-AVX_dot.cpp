/* Copyright (C) T.Hishinuma, All rights reserved.*/
#include<DD-AVX.hpp>
#include<stdio.h>

//1//////////////////////////////////////////////////////////
void DD_AVX_dot(D_Vector vx, D_Vector vy, D_Scalar* val){
#ifdef ddavx_debug
   printf("D scalar = dot(D vector, D vector)\n");
#endif
   D_Scalar tval;
   DD_AVX_dot_D(vx,vy, &tval);
   val->hi = tval.hi;
}
//2//////////////////////////////////////////////////////////

void DD_AVX_dot(DD_Vector vx, D_Vector vy, D_Scalar* val){
#ifdef ddavx_debug
   printf("D scalar = dot(DD vector, D vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(DD(vx),DD(vy), &tval);
   val->hi = tval.hi;
}
//3//////////////////////////////////////////////////////////

void DD_AVX_dot(D_Vector vx, DD_Vector vy, D_Scalar* val){
#ifdef ddavx_debug
   printf("D scalar = dot(D vector, DD vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(DD(vx),DD(vy), &tval);
   val->hi = tval.hi;
}
//4//////////////////////////////////////////////////////////

void DD_AVX_dot(DD_Vector vx, DD_Vector vy, D_Scalar* val){
#ifdef ddavx_debug
   printf("D scalar = dot(DD vector, DD vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(DD(vx),DD(vy), &tval);
   val->hi = tval.hi;
}
//5//////////////////////////////////////////////////////////

void DD_AVX_dot(D_Vector vx, D_Vector vy, DD_Scalar* val){
#ifdef ddavx_debug
   printf("DD scalar = dot(D vector, D vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(DD(vx),DD(vy), &tval);
   val->hi = tval.hi;
   val->lo = tval.lo;
}
//6//////////////////////////////////////////////////////////

void DD_AVX_dot(DD_Vector vx, D_Vector vy, DD_Scalar* val){
#ifdef ddavx_debug
   printf("DD scalar = dot(DD vector, D vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(DD(vx),DD(vy), &tval);
   val->hi = tval.hi;
   val->lo = tval.lo;
}
//7//////////////////////////////////////////////////////////

void DD_AVX_dot(D_Vector vx, DD_Vector vy, DD_Scalar* val){
#ifdef ddavx_debug
   printf("DD scalar = dot(D vector, DD vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(DD(vx),DD(vy), &tval);
   val->hi = tval.hi;
   val->lo = tval.lo;
}
//8//////////////////////////////////////////////////////////

void DD_AVX_dot(DD_Vector vx, DD_Vector vy, DD_Scalar* val){
#ifdef ddavx_debug
   printf("DD scalar = dot(DD vector, DD vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(DD(vx),DD(vy), &tval);
   val->hi = tval.hi;
   val->lo = tval.lo;
}
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

//1//////////////////////////////////////////////////////////
D_Scalar D_Scalar::dot(D_Vector vx, D_Vector vy){
#ifdef ddavx_debug
   printf("D scalar = dot(D vector, D vector)\n");
#endif
   D_Scalar tval;
   DD_AVX_dot_D(vx,vy, &tval);
   this->hi = tval.hi;
   return *this;
}
//2//////////////////////////////////////////////////////////

D_Scalar D_Scalar::dot(DD_Vector vx, D_Vector vy){
#ifdef ddavx_debug
   printf("D scalar = dot(DD vector, D vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(DD(vx),DD(vy), &tval);
   this->hi = tval.hi;
   return *this;
}
//3//////////////////////////////////////////////////////////

D_Scalar D_Scalar::dot(D_Vector vx, DD_Vector vy){
#ifdef ddavx_debug
   printf("D scalar = dot(D vector, DD vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(DD(vx),DD(vy), &tval);
   this->hi = tval.hi;
   return *this;
}
//4//////////////////////////////////////////////////////////

D_Scalar D_Scalar::dot(DD_Vector vx, DD_Vector vy){
#ifdef ddavx_debug
   printf("D scalar = dot(DD vector, DD vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(DD(vx),DD(vy), &tval);
   this->hi = tval.hi;
   return *this;
}
//5//////////////////////////////////////////////////////////

DD_Scalar DD_Scalar::dot(D_Vector vx, D_Vector vy){
#ifdef ddavx_debug
   printf("DD scalar = dot(D vector, D vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(DD(vx),DD(vy), &tval);
   this->hi = tval.hi;
   this->lo = tval.lo;
   return *this;
}
//6//////////////////////////////////////////////////////////

DD_Scalar DD_Scalar::dot(DD_Vector vx, D_Vector vy){
#ifdef ddavx_debug
   printf("DD scalar = dot(DD vector, D vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(DD(vx),DD(vy), &tval);
   this->hi = tval.hi;
   this->lo = tval.lo;
   return *this;
}
//7//////////////////////////////////////////////////////////

DD_Scalar DD_Scalar::dot(D_Vector vx, DD_Vector vy){
#ifdef ddavx_debug
   printf("DD scalar = dot(D vector, DD vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(DD(vx),DD(vy), &tval);
   this->hi = tval.hi;
   this->lo = tval.lo;
   return *this;
}
//8//////////////////////////////////////////////////////////
DD_Scalar DD_Scalar::dot(DD_Vector vx, DD_Vector vy){
#ifdef ddavx_debug
   printf("DD scalar = dot(DD vector, DD vector)\n");
#endif
   DD_Scalar tval;
   DD_AVX_dot_DD(vx, vy, &tval);
   this->hi = tval.hi;
   this->lo = tval.lo;
   return *this;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
void DD_AVX_dot_D(D_Vector vx, D_Vector vy, D_Scalar *val){
#ifdef ddavx_debug
   printf("    compute dot in double\n");
#endif

   int i,n;
   double *x,*y, dot;

   n = vx.N;

   x = vx.hi;
   y = vy.hi;
   dot = val->hi;

#pragma omp parallel for reduction(+: dot)
      for(i=0; i<n; i++)
      {
	 dot += x[i] * y[i];
      }
   val->hi = dot;
}


void DD_AVX_dot_DD(DD_Vector vx, DD_Vector vy, DD_Scalar* val){
#ifdef ddavx_debug
	printf("   compute dot in DD\n");
#endif
	int i,j,n;
	int is,ie,nprocs,my_rank;
	double *gt,*gtavx;
	double *x,*y,*xl,*yl;
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
	y  = vy.hi;
	xl = vx.lo;
	yl = vy.lo;

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
			DD_AVX_FMAN_AVX(a[0],a[4],a[0],a[4],y[i],yl[i],x[i],xl[i]);
		}

// 		DD_AVX_ADD_SSE2(a[0], a[4], a[0], a[4], a[0], a[4]);
// 		DD_AVX_ADD_SSE2(a[0], a[4], a[1], a[5], a[1], a[5]);
// 		DD_AVX_ADD_SSE2(a[0], a[4], a[2], a[6], a[2], a[6]);
// 		DD_AVX_ADD_SSE2(a[0], a[4], a[3], a[7], a[3], a[7]);
		DD_AVX_HADDALL_AVX(a[0],a[4],a[0],a[4]);
		for(;i<ie;i++)
		{
			DD_AVX_FMA_SSE2(a[0],a[4],a[0],a[4],y[i],yl[i],x[i],xl[i]);
		}

 		for(j=0;j<8;j++){
 			gtavx[my_rank*8+j]=a[j];
 		}
#elif defined(USE_FMA2_SSE2)
		gt[my_rank*2] = gt[my_rank*2+1] = 0.0;
		gt[my_rank*2+2] = gt[my_rank*2+3] = 0.0;
		for(i=is;i<ie-1;i+=2)
		{
			DD_AVX_FMA2_SSE2(gt[my_rank*2],gt[my_rank*2+2],gt[my_rank*2],gt[my_rank*2+2],\
					y[i],yl[i],x[i],xl[i]);
		}
		DD_AVX_ADD_SSE2(gt[my_rank*2],gt[my_rank*2+1],gt[my_rank*2],gt[my_rank*2+1],gt[my_rank*2+2],gt[my_rank*2+3]);
		for(;i<ie;i++)
		{
			DD_AVX_FMA_SSE2(gt[my_rank*2],gt[my_rank*2+1],gt[my_rank*2],gt[my_rank*2+1],\
					y[i],yl[i],x[i],xl[i]);
		}
#else
		gt[my_rank*2] = gt[my_rank*2+1] = 0.0;
		for(i=is;i<ie;i++)
		{
			DD_AVX_FMA(gt[my_rank*2],gt[my_rank*2+1],gt[my_rank*2],gt[my_rank*2+1],\
					y[i],yl[i],x[i],xl[i]);
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
	val->hi = dotm.hi;
	val->lo = dotm.lo;
}
