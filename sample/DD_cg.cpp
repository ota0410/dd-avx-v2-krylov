#include <DD-AVX.hpp>
#include <iostream>
#include <omp.h>
#include <string.h>

using std::cout;
using std::endl;

#define SOLVE 0
#define TIME 1
#define N A.N
#define TOL 1.0e-10
#define maxitr N

DD_Vector SpMV( D_Matrix A,DD_Vector x, DD_Vector y ) 
{
  int i,j;
  DD_Scalar k;
  DD_Scalar temp;
  int size = x.getsize();
  y.broadcast( 0.0 );
  for (i = 0; i < size; i++) {
    temp = 0.0;
    for (j = A.row[i]; j < A.row[i+1]; j++){
      temp = temp + A.val[j] * x.getelm( A.col[j] );
    }
    y.chgelm( temp,i );
  }
  return y;
}

void cg( D_Matrix A )
{
  DD_Vector x;
  DD_Vector b;
  DD_Vector r;
  DD_Vector p;
  DD_Vector q;
  DD_Vector vec_tmp;
  DD_Vector init;

  x.malloc( N );
  b.malloc( N );
  r.malloc( N );
  p.malloc( N );
  q.malloc( N );
  vec_tmp.malloc( N );
  init.malloc( N );
  
  x.broadcast( 0 );
  
  for (int i = 0; i < N; i++){
    init.hi[i] = 1.0 * pow(-1.0,i);
    init.lo[i] = 0.0;
  }
  //  init.broadcast( 1.0 );

  b.broadcast( 0 );
  D_Scalar resid;
  D_Scalar r_nrm2;
  D_Scalar b_nrm2;

  DD_Scalar alpha;
  DD_Scalar beta;
  DD_Scalar a_scl;
  DD_Scalar b_scl;
  DD_Scalar c_scl;
  DD_Scalar temp;

  DD_AVX_SpMV( A,init,b );
  // b.print_all();
  DD_AVX_SpMV( A,x,vec_tmp );
  DD_AVX_xpay( b,( DD_Scalar )(-1.0),vec_tmp );
  r.copy( vec_tmp );
  p.copy( r );

  DD_AVX_nrm2( r, &r_nrm2 );
  DD_AVX_nrm2( b, &b_nrm2 );
  resid = r_nrm2 / b_nrm2;
  
  int count = 1;
  // cout << "---------------------------------相対残差----------------------------------"  << endl;

  while ( count < maxitr ) {

    //    resid.print();
    SpMV( A,p,q );
    a_scl.dot( r,r );
    b_scl.dot( p,q );

    alpha = a_scl / b_scl;
    
    DD_AVX_axpy( alpha,p,x );
    DD_AVX_axpy( -alpha,q,r );

    if( resid < TOL )
      break;
    
    c_scl.dot( r,r );
    beta = c_scl / a_scl;
    DD_AVX_xpay( r,beta,p );
    
    DD_AVX_nrm2( r,&r_nrm2 );
    resid = r_nrm2 / b_nrm2 ;
    count++;
  }
  cout << count << ",";
  resid.print();
  
  /*
  cout << "\n計測時間は";
  cout << end - start << endl;
  */
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
  cg( A );
  end = omp_get_wtime();
  if ( TIME )
    printf("%f\n",end-start);
  return 0;
}
