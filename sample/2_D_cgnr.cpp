#include <DD-AVX.hpp>
#include <iostream>
#include <string.h>
#include <cstdlib>
#include <ctime>
#include <omp.h>

using std::cout;
using std::endl;

#define N A.N
#define maxitr N
#define TOL 1.0e-10
#define SOLVE 0
#define SWEEP 0
#define OMEGA 1.22

D_Scalar dot( D_Vector x, D_Vector y )
{
  int i;
  D_Scalar scl = 0.0;

  for ( i = 0; i < x.getsize() ; i++ ) 
    scl = scl + x.hi[i] * y.hi[i];
    
  return scl;
}

void cgnr( D_Matrix A )
{
  D_Vector x;
  D_Vector b;
  D_Vector bb;
  D_Vector r;
  D_Vector p;
  D_Vector q;
  D_Vector z;
  D_Vector init;

  x.malloc( N );
  b.malloc( N );
  bb.malloc( N );
  r.malloc( N );
  p.malloc( N );
  q.malloc( N );
  z.malloc( N );
  init.malloc( N );
  
  x.broadcast( 0 );
  //  init.broadcast( 1 );
  
  for (int i = 0; i < N; i++)
    init.hi[i] = pow( -1.0,i );
  
  /*    
  srand((unsigned)time(NULL));  
  for (int i = 0; i < N; i++) 
    init.hi[i] = (double)rand()/(double)RAND_MAX;  
  */
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
  D_Scalar e;
  D_Scalar y;

  DD_AVX_SpMV( A,init,b );
  DD_AVX_TSpMV( A,b,bb );
  DD_AVX_SpMV( A,x,z );
  DD_AVX_xpay( b,( D_Scalar )(-1.0),z );

  r.copy( z );
  /********************************/
  if ( SWEEP ){
    D_Scalar del;
    D_Vector e,f;
    e.malloc( N );
    f.malloc( N );
    for(int i = 0; i < N; i++){
      e.hi[i] = 1.0;
      DD_AVX_SpMV( A,e,f );
      del = OMEGA * ( dot( f,r ) / dot( f,f ) );
      DD_AVX_axpy( del,e,x );
      DD_AVX_axpy( -del,f,r );
      e.hi[i] = 0.0;
    }
  }
  /********************************/
  DD_AVX_TSpMV( A,r,z );
  p.copy( z );

  DD_AVX_nrm2( r, &r_nrm2 );
  DD_AVX_nrm2( b, &b_nrm2 );
  resid = r_nrm2 / b_nrm2;
  int count = 1;
  //  cout << "---------------------------------相対残差----------------------------------"  << endl;
  /*
  D_Vector error;
  D_Scalar serror;
  error.malloc( N );
  DD_AVX_axpyz( ( D_Scalar )(-1.0), init, x, error );
  DD_AVX_nrm2( error, &serror );
  serror.print();
  */
  while ( resid > TOL && count < maxitr ) {
    resid.print();
    //serror.print();
    DD_AVX_SpMV( A,p,q );
    
    a_scl = dot( z,z );//内積
    b_scl = dot( q,q );//内積

    alpha = a_scl / b_scl;

    DD_AVX_axpy( alpha,p,x );
    DD_AVX_axpy( -alpha,q,r );
    DD_AVX_TSpMV( A,r,z );
    
    c_scl = dot( z,z );//内積
    beta = c_scl / a_scl;

    DD_AVX_xpay( z,beta,p );
    DD_AVX_nrm2( r,&r_nrm2 );
    //    DD_AVX_nrm2( z,&r_nrm2 );
    resid = r_nrm2 / b_nrm2 ;

    // DD_AVX_axpyz( ( D_Scalar )(-1.0), init, x, error );
    // DD_AVX_nrm2( error, &serror );
    
    count++;
  }
 

  // cout << "ループの計測時間は";
  /*
  cout << end - start << endl; 
  */
  if ( SOLVE ) {
    cout << "----------------------------------解--------------------------------------\n" << endl;
    x.print_all();
    cout << "---------------------------------------------------------------------------"  << endl;
  }
  cout << count << "," ;
  resid.print();
  //  serror.print();
  
  b.free();
  x.free();
  r.free();
  p.free();
  q.free();
}

int main ( int argc, char* argv[] )
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
  start = omp_get_wtime();
  cgnr( A );
  end = omp_get_wtime();
  cout << "%実行時間は" << end - start << endl;
  
  return 0;
}
