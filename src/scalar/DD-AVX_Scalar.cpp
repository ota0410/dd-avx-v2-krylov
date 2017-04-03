/* Copyright (C) T.Hishinuma, All rights reserved.*/
#include<DD-AVX.hpp>
#include<stdio.h>
#include<string.h>
///////////// DOUBLE ////////////////
///////////// DOUBLE ////////////////
///////////// DOUBLE ////////////////
//arithmetic operator
D_Scalar D_Scalar::operator=(double a) {
#ifdef ddavx_debug
   printf("D_Scalar::operator=double\n");
#endif
   this->hi = a;
   return *this;
}

D_Scalar D_Scalar::operator=(const D_Scalar& D) {
#ifdef ddavx_debug
   printf("D_Scalar::operator=D\n");
#endif
   this->hi = D.hi;
   return *this;
}

D_Scalar D_Scalar::operator=(const DD_Scalar& DD) {
#ifdef ddavx_debug
   printf("D_Scalar::operator=DD\n");
#endif
   this->hi = DD.hi;

   return *this;
}

D_Scalar D_Scalar::operator-(void) {
#ifdef ddavx_debug
   printf("D_Scalar::operator-\n");
#endif
   D_Scalar tmp;
   tmp.hi = -1 * this->hi;
   return tmp;
}

//arithmetic operator
/////////////////////////////////////////////////////////////////////////
D_Scalar D_Scalar::operator+(double rhv) {
#ifdef ddavx_debug
   printf("[D  ] D_Scalar D_Scalar::operator+(double)\n");
#endif

   this->hi = this->hi + rhv;
   return *this;
}

D_Scalar D_Scalar::operator+(D_Scalar rhv) {
#ifdef ddavx_debug
   printf("[D  ] D_Scalar D_Scalar::operator+(D)\n");
#endif

   this->hi = this->hi + rhv.hi;
   return *this;
}

DD_Scalar D_Scalar::operator+(DD_Scalar rhv) {
#ifdef ddavx_debug
   printf("[Mix] DD_Scalar D_Scalar::operator+(DD)\n");
#endif
   DD_Scalar r;
   DD_AVX_ADD(r.hi, r.lo, this->hi, 0.0, rhv.hi,rhv.lo);

   return r;
}
/////////////////////////////////////////////////////////////////////////
//
D_Scalar D_Scalar::operator-(double rhv) {
#ifdef ddavx_debug
   printf("[D  ] D_Scalar D_Scalar::operator-(Double)\n");
#endif

   this->hi = this->hi - rhv;
   return *this;
}

D_Scalar D_Scalar::operator-(D_Scalar rhv) {
#ifdef ddavx_debug
   printf("[D  ] D_Scalar D_Scalar::operator-(D)\n");
#endif

   this->hi = this->hi - rhv.hi;
   return *this;
}

DD_Scalar D_Scalar::operator-(DD_Scalar rhv) {
#ifdef ddavx_debug
   printf("[Mix] DD_Scalar D_Scalar::operator-(DD)\n");
#endif

   DD_Scalar r;
   DD_AVX_ADD(r.hi, r.lo, this->hi, 0.0, -rhv.hi,-rhv.lo);

   return r;
}
/////////////////////////////////////////////////////////////////////////

D_Scalar D_Scalar::operator*(double rhv) {
#ifdef ddavx_debug
   printf("[D  ] D_Scalar D_Scalar::operator*(double)\n");
#endif

   this->hi = this->hi * rhv;
   return *this;
}

D_Scalar D_Scalar::operator*(D_Scalar rhv) {
#ifdef ddavx_debug
   printf("[D  ] D_Scalar D_Scalar::operator*(D)\n");
#endif

   this->hi = this->hi * rhv.hi;
   return *this;
}

DD_Scalar D_Scalar::operator*(DD_Scalar rhv) {
#ifdef ddavx_debug
   printf("[Mix] DD_Scalar D_Scalar::operator*(DD)\n");
#endif

   DD_Scalar r;
   DD_AVX_MULD(r.hi, r.lo, rhv.hi,rhv.lo, this->hi);
   return r;
}
/////////////////////////////////////////////////////////////////////////

D_Scalar D_Scalar::operator/(double rhv) {
#ifdef ddavx_debug
   printf("[D  ] D_Scalar D_Scalar::operator/(double)\n");
#endif

   this->hi = this->hi / rhv;
   return *this;
}

D_Scalar D_Scalar::operator/(D_Scalar rhv) {
#ifdef ddavx_debug
   printf("[D  ] D_Scalar D_Scalar::operator/(D)\n");
#endif

   this->hi = this->hi / rhv.hi;
   return *this;
}

DD_Scalar D_Scalar::operator/(DD_Scalar rhv) {
#ifdef ddavx_debug
   printf("[Mix] DD_Scalar D_Scalar::operator/(DD)\n");
#endif

   DD_Scalar r;
   DD_AVX_DIV(r.hi, r.lo, this->hi, 0.0, rhv.hi,rhv.lo);
   return r;
}
/////////////////////////////////////////////////////////////////////////

D_Scalar::operator double() {
#ifdef ddavx_debug
   printf("D_Scalar::cast operator double()\n");
#endif
   double a;
   a = this->hi;
   return a;
}

D_Scalar::operator DD_Scalar() {
#ifdef ddavx_debug
   printf("D_Scalar::cast operator DD_Scalar()\n");
#endif
   DD_Scalar a;
   a.hi = this->hi;
   a.lo = 0;
   return a;
}

///////////// Double-Double ////////////////
///////////// Double-Double ////////////////
///////////// Double-Double ////////////////
///////////// Double-Double ////////////////
//!binary operator
DD_Scalar DD_Scalar::operator=(double a) {
#ifdef ddavx_debug
   printf("DD_Scalar::operator=double\n");
#endif
   this->hi = a;
   this->lo = 0;
   return *this;
}
DD_Scalar DD_Scalar::operator=(const D_Scalar& D) {
#ifdef ddavx_debug
   printf("DD_Scalar::operator=D\n");
#endif
   this->hi = D.hi;
   this->lo = 0;

   return *this;
}
DD_Scalar DD_Scalar::operator=(const DD_Scalar& DD) {
#ifdef ddavx_debug
   printf("DD_Scalar::operator=DD\n");
#endif
   this->hi = DD.hi;
   this->lo = DD.lo;

   return *this;
}
DD_Scalar DD_Scalar::operator-(void) {
#ifdef ddavx_debug
   printf("D_Scalar::operator-\n");
#endif
   DD_Scalar tmp;
   tmp.hi = -1 * this->hi;
   tmp.lo = -1 * this->lo;
   return tmp;
}
/////////////////////////////////////////////////////
//arithmetic operator

DD_Scalar DD_Scalar::operator+(double rhv) {
#ifdef ddavx_debug
   printf("[Mix] DD_Scalar DD_Scalar::operator+(double)\n");
#endif
   DD_Scalar r;
   DD_AVX_ADD(r.hi, r.lo, this->hi, this->lo, rhv, 0.0);
   return r;
}

DD_Scalar DD_Scalar::operator+(D_Scalar rhv) {
#ifdef ddavx_debug
   printf("[Mix] DD_Scalar DD_Scalar::operator+(D)\n");
#endif
   DD_Scalar r;
   DD_AVX_ADD(r.hi, r.lo, this->hi, this->lo, rhv.hi, 0.0);
   return r;
}

DD_Scalar DD_Scalar::operator+(DD_Scalar rhv) {
#ifdef ddavx_debug
   printf("[DD ] DD_Scalar DD_Scalar::operator+(DD)\n");
#endif
   DD_Scalar r;
   DD_AVX_ADD(r.hi, r.lo, this->hi, this->lo, rhv.hi,rhv.lo);

   return r;
}
/////////////////////////////////////////////////////////////////////////
 
DD_Scalar DD_Scalar::operator-(double rhv) {
#ifdef ddavx_debug
   printf("[Mix] DD_Scalar DD_Scalar::operator-(double)\n");
#endif

   DD_Scalar r;
   DD_AVX_ADD(r.hi, r.lo, this->hi, this->lo, -rhv,0.0);

   return r;
}

DD_Scalar DD_Scalar::operator-(D_Scalar rhv) {
#ifdef ddavx_debug
   printf("[Mix] DD_Scalar DD_Scalar::operator-(D)\n");
#endif

   DD_Scalar r;
   DD_AVX_ADD(r.hi, r.lo, this->hi, this->lo, -rhv.hi,0.0);

   return r;
}

DD_Scalar DD_Scalar::operator-(DD_Scalar rhv) {
#ifdef ddavx_debug
   printf("[DD ] DD_Scalar DD_Scalar::operator-(DD)\n");
#endif

   DD_Scalar r;
   DD_AVX_ADD(r.hi, r.lo, this->hi, this->lo, -rhv.hi,-rhv.lo);

   return r;
}

/////////////////////////////////////////////////////////////////////////

DD_Scalar DD_Scalar::operator*(double rhv) {
#ifdef ddavx_debug
   printf("[Mix] DD_Scalar DD_Scalar::operator*(double)\n");
#endif

   DD_Scalar r;
   DD_AVX_MULD(r.hi, r.lo, this->hi, this->lo, rhv);

   return r;
}

DD_Scalar DD_Scalar::operator*(D_Scalar rhv) {
#ifdef ddavx_debug
   printf("[Mix] DD_Scalar DD_Scalar::operator*(D)\n");
#endif

   DD_Scalar r;
   DD_AVX_MULD(r.hi, r.lo, this->hi, this->lo, rhv.hi);

   return r;
}

DD_Scalar DD_Scalar::operator*(DD_Scalar rhv) {
#ifdef ddavx_debug
   printf("[DD ] DD_Scalar DD_Scalar::operator*(DD)\n");
#endif

   DD_Scalar r;
   DD_AVX_MUL(r.hi, r.lo, this->hi, this->lo, rhv.hi,rhv.lo);
   return r;
}
/////////////////////////////////////////////////////////////////////////

DD_Scalar DD_Scalar::operator/(double rhv) {
#ifdef ddavx_debug
   printf("[Mix] DD_Scalar DD_Scalar::operator/(double)\n");
#endif

   DD_Scalar r;
   DD_AVX_DIV(r.hi, r.lo, this->hi, this->lo, rhv,0.0);

   return r;
}

DD_Scalar DD_Scalar::operator/(D_Scalar rhv) {
#ifdef ddavx_debug
   printf("[Mix] DD_Scalar DD_Scalar::operator/(D)\n");
#endif

   DD_Scalar r;
   DD_AVX_DIV(r.hi, r.lo, this->hi, this->lo, rhv.hi,0.0);

   return r;
}

DD_Scalar DD_Scalar::operator/(DD_Scalar rhv) {
#ifdef ddavx_debug
   printf("[DD ] DD_Scalar DD_Scalar::operator/(DD)\n");
#endif

   DD_Scalar r;
   DD_AVX_DIV(r.hi, r.lo, this->hi, this->lo, rhv.hi,rhv.lo);
   return r;
}

/////////////////////////////////////////////////////////////////////////

DD_Scalar::operator D_Scalar() {
#ifdef ddavx_debug
   printf("DD_Scalar::cast operator D_Scalar()\n");
#endif
   D_Scalar a;
   a.hi = this->hi;
   return a;
}
DD_Scalar::operator double() {
#ifdef ddavx_debug
   printf("DD_Scalar::cast operator double()\n");
#endif
   double a;
   a = this->hi;
   return a;
}

void D_Scalar::hello( D_Scalar alpha ){
#ifdef ddavx_debug
  printf("hello %f\n",alpha);
#endif
  printf("hello %f\n",alpha);
}
