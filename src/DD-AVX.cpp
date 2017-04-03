/* Copyright (C) T.Hishinuma, All rights reserved.*/
#include<DD-AVX.hpp>
#include<stdio.h>
#include<string.h>
#include<ctype.h>
void DD_AVX_version(){
   PRINT_VERSION;
}

///////////// Scalar Cast ////////////////
// DD -> D////////////////////////////////////////////
DD_Scalar DD(double a){
#ifdef ddavx_debug
	 printf("DD(): cast scalar double->DD\n");
#endif
   DD_Scalar b;
   b.hi = a;
   b.lo = 0.0;
   return b;
}

DD_Scalar DD(D_Scalar a){
#ifdef ddavx_debug
	 printf("DD(): cast scalar D->DD\n");
#endif
   DD_Scalar b;
   b.hi = a.hi;
   b.lo = 0.0;
   return b;
}

DD_Scalar DD(DD_Scalar a){
#ifdef ddavx_debug
	 printf("DD(): cast scalar DD->DD (nothing to do)\n");
#endif
   return a;
}

// D -> DD/////////////////////////////////////////////////
D_Scalar D(double a){
#ifdef ddavx_debug
	 printf("D(): cast scalar Double->D\n");
#endif
	 D_Scalar x;
	 x.hi = a;
	 return x;
}

D_Scalar D(D_Scalar a){
#ifdef ddavx_debug
	 printf("D(): cast scalar D->D (nothing to do)\n");
#endif
	 return a;
}

D_Scalar D(DD_Scalar a){
#ifdef ddavx_debug
	 printf("D(): cast scalar DD->D\n");
#endif
	 D_Scalar x;
	 x.hi = a.hi;
	 return x;
}

//Cast DD_Vector/////////////////////////////////////////////////
//
// DD -> D
DD_Vector DD(D_Vector a){
#ifdef ddavx_debug
	 printf("DD(): cast vector D->D\n");
#endif
   DD_Vector b;
   b.malloc(a.N);

   for(int i=0;i<a.N;i++){
   b.hi[i] = a.hi[i];
   b.lo[i] = 0.0;
   }
   return b;
}

DD_Vector DD(DD_Vector a){
#ifdef ddavx_debug
	 printf("DD(): cast vector DD->DD (nothing to do)\n");
#endif
   return a;
}

// D -> DD
D_Vector D(D_Vector a){
#ifdef ddavx_debug
	 printf("D(): cast vector D->D (nothing to do)\n");
#endif
   return a;
}

D_Vector D(DD_Vector a){
#ifdef ddavx_debug
	 printf("D(): cast vector DD->D\n");
#endif
   D_Vector b;
   b.malloc(a.N);
   for(int i=0;i<a.N;i++){
   b.hi[i] = a.hi[i];
   }
   return b;
}

///////////// Vector member function ////////////////
///////////// D_Vector ////////////////
// input and output functions are written in System.c
D_Vector D_Vector::operator=(const D_Vector& D) {
#ifdef ddavx_debug
   printf("D_Vector::operator=(D_Vector)\n");
#endif
   if(this->N != D.N){
      printf("error D.N != D.N in cast operator");
      abort();
   }
#pragma omp parallel for
   for(int i=0; i<N; i++){
      this->hi[i] = D.hi[i];
   }
   return *this;
}

D_Vector D_Vector::operator=(const DD_Vector& DD) {
#ifdef ddavx_debug
   printf("D_Vector::operator=(DD_Vector)\n");
#endif
   if(this->N != DD.N){
      printf("error D.N != DD.N in cast operator");
      abort();
   }
#pragma omp parallel for
   for(int i=0; i<N; i++){
      this->hi[i] = DD.hi[i];
   }
   return *this;
}

D_Vector D_Vector::copy(D_Vector D) {
#ifdef ddavx_debug
   printf("D_Vector::copy(D_Vector)\n");
#endif
   if(this->N != D.N){
      printf("error D.N != D.N in cast operator");
      abort();
   }
#pragma omp parallel for
   for(int i=0; i<N; i++){
      this->hi[i] = D.hi[i];
   }
   return *this;
}

D_Vector D_Vector::copy(DD_Vector DD) {
#ifdef ddavx_debug
   printf("D_Vector::operator=(DD_Vector)\n");
#endif
   if(this->N != DD.N){
      printf("error D.N != DD.N in cast operator");
      abort();
   }
#pragma omp parallel for
   for(int i=0; i<N; i++){
      this->hi[i] = DD.hi[i];
   }
   return *this;
}

void D_Vector::print(int n){
#ifdef ddavx_debug
	 printf("D_vector::print()\n");
#endif
   printf("%e\n", hi[n]);
}

void D_Vector::print_all(){
#ifdef ddavx_debug
	 printf("D_vector::print_all()\n");
#endif
   for(int i=0;i<N;i++){
   	printf("%e\n", hi[i]);
   }
}


int D_Vector::getsize (){
#ifdef ddavx_debug
	 printf("D_vector::getsize()\n");
#endif
   return N;
}

void D_Vector::malloc(int n){
#ifdef ddavx_debug
	 printf("D_vector::malloc()\n");
#endif
   hi = new double[n];
   N = n;
}

void D_Vector::free(){
   delete hi;
}

void D_Vector::broadcast(double val){
#ifdef ddavx_debug
	 printf("D_Vector::broadcast(Double)\n");
#endif
#pragma omp parallel for
   for(int i=0;i<N;i++){
	hi[i] = val;
   }

}

void D_Vector::broadcast(D_Scalar val){
#ifdef ddavx_debug
	 printf("D_Vector::broadcast(D)\n");
#endif
#pragma omp parallel for
   for(int i=0;i<N;i++){
	hi[i] = val.hi;
   }

}

void D_Vector::broadcast(DD_Scalar val){
#ifdef ddavx_debug
	 printf("D_Vector::broadcast(DD)\n");
#endif
#pragma omp parallel for
   for(int i=0;i<N;i++){
	hi[i] = val.hi;
   }
}
///////////// DD_Vector ////////////////
// input and output functions are written in System.c
DD_Vector DD_Vector::operator=(const D_Vector& D) {
#ifdef ddavx_debug
   printf("DD_Vector::operator=D_Vector\n");
#endif
   if(this->N != D.N){
      printf("error DD.N != D.N in cast operator");
      abort();
   }
#pragma omp parallel for
   for(int i=0; i<N; i++){
      this->hi[i] = D.hi[i];
      this->lo[i] = 0;
   }
   return *this;
}

DD_Vector DD_Vector::operator=(const DD_Vector& DD) {
#ifdef ddavx_debug
   printf("DD_Vector::operator=DD_Vector\n");
#endif
   if(this->N != DD.N){
      printf("error DD.N != DD.N in cast operator");
      abort();
   }
#pragma omp parallel for
   for(int i=0; i<N; i++){
      this->hi[i] = DD.hi[i];
      this->lo[i] = DD.lo[i];
   }
   return *this;
}

DD_Vector DD_Vector::copy(D_Vector D) {
#ifdef ddavx_debug
   printf("DD_Vector::operator=D_Vector\n");
#endif
   if(this->N != D.N){
      printf("error DD.N != D.N in cast operator");
      abort();
   }
#pragma omp parallel for
   for(int i=0; i<N; i++){
      this->hi[i] = D.hi[i];
      this->lo[i] = 0;
   }
   return *this;
}

DD_Vector DD_Vector::copy(DD_Vector DD) {
#ifdef ddavx_debug
   printf("DD_Vector::operator=DD_Vector\n");
#endif
   if(this->N != DD.N){
      printf("error DD.N != DD.N in cast operator");
      abort();
   }
#pragma omp parallel for
   for(int i=0; i<N; i++){
      this->hi[i] = DD.hi[i];
      this->lo[i] = DD.lo[i];
   }
   return *this;
}

void DD_Vector::print(int n){
#ifdef ddavx_debug
	 printf("DD_Vector::print()\n");
#endif
   printf("hi = %e lo = %e\n", hi[n], lo[n]);
}

void DD_Vector::print_all(){
#ifdef ddavx_debug
	 printf("DD_Vector::print_all()\n");
#endif
   for(int i=0;i<N;i++){
   	printf("hi = %e lo = %e\n", hi[i], lo[i]);
   }
}


int DD_Vector::getsize(){
#ifdef ddavx_debug
	 printf("DD_Vector::get_size()\n");
#endif
   return N;
}

void DD_Vector::malloc(int n){
#ifdef ddavx_debug
	 printf("DD_vector::malloc()\n");
#endif
   hi = new double[n*2];
   lo = &hi[n];
   N = n;
}

void DD_Vector::free(){
   delete hi;
}

void DD_Vector::broadcast(double val){
#ifdef ddavx_debug
	 printf("DD_vector::broadcast(double)\n");
#endif
#pragma omp parallel for
   for(int i=0;i<N;i++){
	hi[i] = val;
	lo[i] = 0;
   }

}

void DD_Vector::broadcast(D_Scalar val){
#ifdef ddavx_debug
	 printf("DD_Vector::broadcast(D)\n");
#endif
#pragma omp parallel for
   for(int i=0;i<N;i++){
	hi[i] = val.hi;
	lo[i] = 0;
   }

}

void DD_Vector::broadcast(DD_Scalar val){
#ifdef ddavx_debug
	 printf("DD_Vector::broadcast(DD)\n");
#endif
#pragma omp parallel for
   for(int i=0;i<N;i++){
	hi[i] = val.hi;
	lo[i] = val.lo;
   }
}

///////////// D_Matrix ////////////////
// input and output functions are written in System.c
/*
void D_Matrix::print_all(){
#ifdef ddavx_debug
	 printf("DD_Matrix::print_all()\n");
#endif
   for(int i=0;i<nnz;i++){
   	printf("%d\t%d\t%e\n",row[i],col[i],val[i]);
   }
}+
*/
