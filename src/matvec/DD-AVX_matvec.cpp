/* Copyright (C) T.Hishinuma, All rights reserved.*/
#include<DD-AVX.hpp>
#include<stdio.h>

////////////////////////////////////////////////////////////
void DD_AVX_SpMV(D_Matrix A, D_Vector vx, D_Vector vy){
#ifdef ddavx_debug
	printf("SpMV(D Matrix, D Vector, D Vector)\n");
#endif
	switch(A.format){
		case 0:
			DD_AVX_SpMV_CRS_D(A,vx,vy);
			break;
		case 1:
			DD_AVX_SpMV_CRS_D(A,vx,vy);
			break;
		case 2:
			DD_AVX_SpMV_CRS_D(A,vx,vy);
			break;
		case 3:
			DD_AVX_SpMV_CRS_D(A,vx,vy);
			break;
		default:
			printf("error, format is %d",A.format);
	}
}
////////////////////////////////////////////////////////////

void DD_AVX_SpMV(D_Matrix A, DD_Vector vx, D_Vector vy){
#ifdef ddavx_debug
	printf("SpMV(D Matrix, DD Vector, D Vector)\n");
#endif
	D_Vector ty;
	switch(A.format){
		case 0:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 1:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 2:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 3:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),D(vy));
			break;
		default:
			printf("error, format is %d",A.format);
	}
}
////////////////////////////////////////////////////////////

void DD_AVX_SpMV(D_Matrix A, D_Vector vx, DD_Vector vy){
#ifdef ddavx_debug
	D_Vector ty;
	switch(A.format){
		case 0:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 1:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 2:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 3:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),D(vy));
			break;
		default:
			printf("error, format is %d",A.format);
	}
	printf("SpMV(D Scalar, D Vector, DD Vector)\n");
#endif
}
////////////////////////////////////////////////////////////

void DD_AVX_SpMV(D_Matrix A, DD_Vector vx, DD_Vector vy){
#ifdef ddavx_debug
	printf("SpMV(D Scalar, DD Vector, DD Vector)\n");
#endif
	D_Vector ty;
	switch(A.format){
		case 0:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 1:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 2:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 3:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),D(vy));
			break;
		default:
			printf("error, format is %d",A.format);
	}
}
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
void DD_AVX_TSpMV(D_Matrix A, D_Vector vx, D_Vector vy){
#ifdef ddavx_debug
	printf("SpMV(D Matrix, D Vector, D Vector)\n");
#endif
	switch(A.format){
		case 0:
			DD_AVX_TSpMV_CRS_D(A,vx,vy);
			break;
		case 1:
			DD_AVX_TSpMV_CRS_D(A,vx,vy);
			break;
		case 2:
			DD_AVX_TSpMV_CRS_D(A,vx,vy);
			break;
		case 3:
			DD_AVX_TSpMV_CRS_D(A,vx,vy);
			break;
		default:
			printf("error, format is %d",A.format);
	}
}
////////////////////////////////////////////////////////////

void DD_AVX_TSpMV(D_Matrix A, DD_Vector vx, D_Vector vy){
#ifdef ddavx_debug
	printf("SpMV(D Matrix, DD Vector, D Vector)\n");
#endif
	D_Vector ty;
	switch(A.format){
		case 0:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 1:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 2:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 3:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		default:
			printf("error, format is %d",A.format);
	}
}
////////////////////////////////////////////////////////////

void DD_AVX_TSpMV(D_Matrix A, D_Vector vx, DD_Vector vy){
#ifdef ddavx_debug
	printf("SpMV(D Scalar, D Vector, DD Vector)\n");
#endif
	D_Vector ty;
	switch(A.format){
		case 0:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 1:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 2:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 3:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		default:
			printf("error, format is %d",A.format);
	}
}
////////////////////////////////////////////////////////////

void DD_AVX_TSpMV(D_Matrix A, DD_Vector vx, DD_Vector vy){
#ifdef ddavx_debug
	printf("SpMV(D Scalar, DD Vector, DD Vector)\n");
#endif
	D_Vector ty;
	switch(A.format){
		case 0:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 1:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 2:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		case 3:
			ty.malloc(vy.N);
			ty.hi = vy.hi;
			DD_AVX_SpMV_CRS_D(A,D(vx),ty);
			break;
		default:
			printf("error, format is %d",A.format);
	}
}
////////////////////////////////////////////////////////////
