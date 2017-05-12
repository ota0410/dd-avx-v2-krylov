#include <DD-AVX.hpp>
#include <time.h>
#include <iostream>
#include <quadmath.h>

using std::cout;
using std::endl;

#define N A.N

int main( int argc, char* argv[] )
{

  char filename1[256] = {'\0'};/*元の行列*/
  char filename2[256] = {'\0'};/*転置行列*/

  sprintf( filename1, "%s",argv[1] );
  sprintf( filename2, "%s",argv[2] );
 
  D_Matrix A,tA;
  A.input( filename1 );
  tA.input( filename2 );

  DD_Vector b;
  DD_Vector init;
  DD_Vector x;
  DD_Vector y;
  DD_Vector r;
  DD_Vector rr;
  DD_Vector p;
  DD_Vector pp;
  DD_Vector q;
  DD_Vector qq;
  //  __float128 r_nrm2;

  init.malloc( N );
  b.malloc( N );
  r.malloc( N );
  rr.malloc( N );
  x.malloc( N );
  y.malloc( N );
  p.malloc( N );
  pp.malloc( N );
  q.malloc( N );
  qq.malloc( N );

  DD_Scalar b_nrm2;
  DD_Scalar resid;
  DD_Scalar r_nrm2;
  DD_Scalar a_scl;
  DD_Scalar b_scl;
  DD_Scalar c_scl;
  DD_Scalar alpha;
  DD_Scalar beta;
  //  D_Scalar param;

  srand( ( unsigned )time( NULL ) );
  
  for (int i = 0; i < N ; i++) {
    x.hi[i] = ( DD_Scalar )( rand() % 10 + 2 );
    x.lo[i] = 0.0;
  }
  
  init.broadcast( 1 );
  b.broadcast( 0 );
  r.broadcast( 0 );
  rr.broadcast( 0 );
  p.broadcast( 0 );
  pp.broadcast( 0 );
  q.broadcast( 0 );
  qq.broadcast( 0 );
  y.broadcast( 0 );

  DD_AVX_SpMV( A,x,y );
  DD_AVX_SpMV( A,init,b );
  init.free();

  DD_AVX_xpay( b,( DD_Scalar )(-1.0),y );
  r.copy( y );
  y.free();

  for (int i = 0; i < N; i++) {
    rr.hi[i] = ( DD_Scalar )( rand() % 10 + 2 );
    rr.lo[i] = 0.0;
    r.lo[i] = 0.0;
  }

  p.copy( r );
  pp.copy( rr );
 
  for (int i = 0; i < N; i++) {
    r_nrm2.hi += (r.hi[i] * r.hi[i]) + 2 * r.hi[i] * r.lo[i];
    r_nrm2.lo += (r.lo[i] * r.lo[i]);
    b_nrm2.hi += (b.hi[i] * b.hi[i]) + 2 * b.hi[i] * b.lo[i];
    b_nrm2.lo += (b.lo[i] * b.lo[i]);
  }

  //  r.print_all();
  //  DD_AVX_nrm2(r,&r_nrm2);
  // r_nrm2.print();

  r_nrm2 = sqrt( r_nrm2 );
  b_nrm2 = sqrt( b_nrm2 );

  resid = (DD_Scalar)r_nrm2 / (DD_Scalar)b_nrm2 ;
  //  D_Scalar z = r_nrm2.hi / b_nrm2.hi;
  // z.print();

  cout << "---------------------------------相対残差--------------------------------------\n" << endl;

  int count = 1;
  while ( resid > pow( 10,-12 ) ) {
    cout << count << ":";
    resid.print();

    DD_AVX_SpMV( A,p,q );
    DD_AVX_SpMV( tA,pp,qq );
    
    a_scl = 0.0;
    b_scl = 0.0;
    for (int i = 0;i < N; i++) {
      a_scl.hi += rr.hi[i] * r.hi[i] + rr.hi[i] * r.lo[i] + rr.lo[i] * r.hi[i];
      a_scl.lo += rr.lo[i] * r.lo[i];
      b_scl.hi += pp.hi[i] * q.hi[i] + pp.hi[i] * q.lo[i] + pp.lo[i] * q.hi[i];
      b_scl.lo += pp.lo[i] * q.lo[i];
    }
    
    alpha = (DD_Scalar)a_scl / (DD_Scalar)b_scl;    
 
    DD_AVX_axpy( alpha,p,x );
    DD_AVX_axpy( -alpha,q,r );
    DD_AVX_axpy( -alpha,qq,rr );

    c_scl = 0.0;
    for (int i = 0;i < N; i++) {
      c_scl.hi += rr.hi[i] * r.hi[i] + rr.hi[i] * r.lo[i] + rr.lo[i] * r.hi[i];
      c_scl.lo += rr.lo[i] * r.lo[i];
    }

    beta =  (DD_Scalar)c_scl / (DD_Scalar)a_scl; 

    DD_AVX_xpay( r,beta,p );
    DD_AVX_xpay( rr,beta,pp );    
    
    r_nrm2 = 0.0;
    for (int i = 0; i < N; i++) {
      r_nrm2.hi += (r.hi[i] * r.hi[i]) + 2 * r.hi[i] * r.lo[i] ;
      r_nrm2.lo += (r.lo[i] * r.lo[i]);
    }

    r_nrm2 = sqrt( r_nrm2 );
    resid = (DD_Scalar)r_nrm2 / (DD_Scalar)b_nrm2;
    count++;

  }

  cout << count << ":" ;
  resid.print();

  cout << "----------------------------------解ベクトル------------------------------------\n" << endl;
  //  x.print_all();
  cout << "--------------------------------------------------------------------------------\n" << endl;
  
  A.free();
  tA.free(); 
  b.free();
  x.free();
  p.free();
  pp.free();
  r.free();
  rr.free();
  q.free();
  qq.free();

  return 0;
}
