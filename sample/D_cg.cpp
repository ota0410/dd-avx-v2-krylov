#include <DD-AVX.hpp>
#include <iostream>
#include <omp.h>
#include <cstring>

using std::cout;
using std::endl;

#define N A.N
#define TOL 1.0e-10
#define maxitr N
#define SOLVE 0

D_Scalar dot( D_Vector x, D_Vector y )
{
  
  D_Scalar scl = 0.0;
  int size = x.getsize();
  for ( int i = 0; i < size ; i++ ) {
    scl = scl + x.hi[i] * y.hi[i];
  }  
  return scl;
}

void cg( D_Matrix A )
{
  D_Vector x;
  D_Vector b;
  D_Vector r;
  D_Vector p;
  D_Vector q;
  D_Vector vec_tmp;
  D_Vector init;

  x.malloc( N );
  b.malloc( N );
  r.malloc( N );
  p.malloc( N );
  q.malloc( N );
  vec_tmp.malloc( N );
  init.malloc( N );
  
  x.broadcast( 0 );
  for(int i = 0; i < N; i++){
    init.hi[i] = pow(-1,i+1);
  }
  b.broadcast( 0 );

  D_Scalar resid;
  D_Scalar r_nrm2;
  D_Scalar b_nrm2;
  D_Scalar alpha;
  D_Scalar beta;
  D_Scalar a_scl;
  D_Scalar b_scl;
  D_Scalar c_scl;
  D_Scalar tmp;

  DD_AVX_SpMV( A,init,b );
  DD_AVX_SpMV( A,x,vec_tmp );
  DD_AVX_xpay( b,(D_Scalar)(-1.0),vec_tmp );
  r.copy( vec_tmp );
  p.copy( r );

  DD_AVX_nrm2( r, &r_nrm2 );
  DD_AVX_nrm2( b, &b_nrm2 );
  resid = r_nrm2 / b_nrm2;
  
  int count = 1;
  //  cout << "---------------------------------相対残差----------------------------------"  << endl;

  while ( count < maxitr ) {
    //    resid.print();
    DD_AVX_SpMV( A,p,q );
    
    a_scl = dot( r,r );
    b_scl = dot( p,q );
    alpha = a_scl / b_scl;
    
    DD_AVX_axpy( alpha,p,x );
    DD_AVX_axpy( -alpha,q,r );

    if ( resid < TOL )
      break;
    
    c_scl = dot( r,r );
    beta = c_scl / a_scl;
    
    DD_AVX_xpay( r,beta,p );
    DD_AVX_nrm2( r,&r_nrm2 );
    resid = r_nrm2 / b_nrm2 ;
    count++;
  }
  cout << count << ",";
  resid.print();
  if ( SOLVE ){
    cout << "----------------------------------解--------------------------------------\n" << endl;
    x.print_all();
    cout << "---------------------------------------------------------------------------"  << endl;
  }

  b.free();
  x.free();
  r.free();
  p.free();
  q.free();
}

int main( int argc, char* argv[] )
{
  char filename[256] = {'\0'};
  char *matrix;
  sprintf( filename,"%s",argv[1] );
  D_Matrix A;
  A.input( filename );
  matrix = strrchr(filename,'/');
  for(int i=1;;i++){
    if(matrix[i] == '.')
      break;
    cout << matrix[i];
  }
  cout << ",";
  
  double start,end;
  // start = omp_get_wtime();
  cg( A );
  // end = omp_get_wtime();
  //  cout << end - start << endl;
  return 0;
}
