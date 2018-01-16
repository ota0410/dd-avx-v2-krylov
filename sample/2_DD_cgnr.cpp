#include <DD-AVX.hpp>
#include <iostream>
#include <cstring>
#include <ctime>
#include <omp.h>

using std::cout;
using std::endl;

#define TOL 1.0e-10
#define N A.N
#define ITR 5000
#define SOLVE 0

DD_Vector TSpMV( D_Matrix A,DD_Vector x, DD_Vector y )
{
  DD_Scalar c1,c2;
  DD_Scalar temp;
  y.broadcast( 0.0 );
  int size = x.getsize();
  for (int i = 0 ; i < size; i++) {
    c1 = x.getelm(i);      
    for (int j = A.row[i]; j < A.row[i+1]; j++){
      c2 = y.getelm( A.col[j] );
      temp =  c2 + ( DD_Scalar )A.val[j] * c1;
      y.chgelm(  temp, A.col[j] );
    }
  }
  return y;
}

DD_Vector SpMV( D_Matrix A,DD_Vector x, DD_Vector y ) 
{
  DD_Scalar c;
  DD_Scalar temp;
  int size = x.getsize();
  int i,j;
  y.broadcast( 0.0 );
  for (i = 0; i < size; i++) {
    temp = 0.0;
    for (j = A.row[i]; j < A.row[i+1]; j++){
      c = x.getelm( A.col[j] );
      temp = temp + ( DD_Scalar )A.val[j] * c;
    }
    y.chgelm( temp,i );
  }
  return y;
}

void cgnr( D_Matrix A )
{

  DD_Vector x;
  DD_Vector b;
  DD_Vector r;
  DD_Vector p;
  DD_Vector q;
  DD_Vector z;
  DD_Vector init;

  x.malloc( N );
  b.malloc( N );
  r.malloc( N );
  p.malloc( N );
  q.malloc( N );
  z.malloc( N );
  init.malloc( N );
  
  x.broadcast( 0.0 );
  //  init.broadcast( 1.0 );
  
  srand((unsigned)time(NULL));
  for (int i = 0; i < N; i++) 
    init.hi[i] = (double)rand()/(double)RAND_MAX;  
  
  /*
  for (int i = 0; i < N; i++){ 
    init.hi[i] = pow(-1,i);
    init.lo[i] = 0.0;
  }
  */
  b.broadcast( 0.0 );
  D_Scalar resid;
  DD_Scalar r_nrm2;
  DD_Scalar b_nrm2;

  DD_Scalar alpha;
  DD_Scalar beta;
  DD_Scalar a_scl;
  DD_Scalar b_scl;
  DD_Scalar c_scl;

  SpMV( A,init,b );
  SpMV( A,x,z );

  DD_AVX_xpay( b,( DD_Scalar )( -1.0 ),z );

  r.copy( z );
  TSpMV( A,r,z );
  p.copy( z );

  DD_AVX_nrm2( r, &r_nrm2 );
  DD_AVX_nrm2( b, &b_nrm2 );
  resid = r_nrm2 / b_nrm2;

  int count = 1;

  while ( count < ITR ) {
    resid.print();
    SpMV( A,p,q );
    
    a_scl.dot( z,z );
    b_scl.dot( q,q );

    alpha = a_scl / b_scl;

    DD_AVX_axpy( alpha,p,x );
    DD_AVX_axpy( -alpha,q,r );

    TSpMV( A,r,z );    
    c_scl.dot( z,z );
    beta = c_scl / a_scl;

    if ( TOL > resid )
      break;

    DD_AVX_xpay( z,beta,p );
    DD_AVX_nrm2( r,&r_nrm2 );
    resid = r_nrm2 / b_nrm2 ;

    count++;
  }
  //itr_end
  if ( SOLVE ) {
    cout << "----------------------------------解--------------------------------------\n" << endl;
    x.print_all();
    cout << "--------------------------------------------------------------------------"  << endl;
  }

  DD_Vector error;
  DD_Scalar serror;
  DD_Scalar sinit;
  error.malloc( N );
  DD_AVX_axpyz( ( DD_Scalar )(-1.0), init, x, error );
  DD_AVX_nrm2( error, &serror );
  DD_AVX_nrm2( init, &sinit );
  serror = serror / sinit;
  
  cout << "error = ";
  serror.print();

  
  cout << count << ",";
  resid.print();
  /*
  if ( ( N - count ) > 0 )
    cout << "%収束しました" << endl;
  else 
    cout << "%収束しませんでした" << endl;
  */
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

  cout << "%計測時間は" << end - start << endl;
  return 0;
}
