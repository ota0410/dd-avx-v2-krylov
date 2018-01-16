#include <DD-AVX.hpp>
#include <iostream>
#include <time.h>
#include <omp.h>

using std::cout;
using std::endl;

#define MAXITR M.N
#define TOL 1.0e-10
#define SOLVE 0

DD_Vector TSpMV( D_Matrix A,DD_Vector x, DD_Vector y )
{
  DD_Scalar c1;
  DD_Scalar c2;
  DD_Scalar temp;
  y.broadcast( 0.0 );
  int size = x.getsize();
  for (int i = 0 ; i < size; i++) {
    c1 = x.getelm( i );      
    for (int j = A.row[i]; j < A.row[i+1]; j++){
      c2 = y.getelm( A.col[j] );
      temp = c2 + ( DD_Scalar )A.val[j] * c1;
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

  for (int i = 0; i < size; i++) {
    temp = 0.0;
    for (int j = A.row[i]; j < A.row[i+1]; j++){
      c = x.getelm( A.col[j] );
      temp = temp + ( DD_Scalar )A.val[j] * c;
    }
    y.chgelm( temp,i );
  }
  return y;
}

D_Matrix SpMM( D_Matrix A )
{
  D_Matrix C;
  int N = A.N;
  C.malloc( N );

  DD_Vector x,y;
  x.malloc( N );
  y.malloc( N );

  int count = 0;
  int h,i,j,k;
  //start for
  for (k=0; k<N; k++) {
    x.broadcast( 0.0 );
    for (i=0; i<N; i++) { //make vector
      for (j=A.row[i]; j<A.row[i+1]; j++) {
	if ( A.col[j] == k ) {
	  x.hi[i] = A.val[j];
	  break;
	} else if ( A.col[j] > k ) {
	break;
	}
      }
    }
    TSpMV( A,x,y );
    //DD_AVX_TSpMV( A,x,y );
#pragma omp prallel for private( h )
    for (h=0; h<N; h++) {
      if ( y.hi[h] != 0 ){
	C.val[count] = y.hi[h];
	C.row[count] = h;
	count++;
      }
    }
    C.col[k+1] = count;
  }
  //End For

  C.nnz = count;
  C.col[N] = C.nnz+1;
  //  printf("%d\n",C.nnz);  
  #pragma omp parallel for
  for (int i=0;i<N+1; i++) 
    C.ptr[i] = C.col[i];
  #pragma omp parallel for
  for (int i=0;i<C.nnz; i++)
    C.col[i] = C.row[i];
  #pragma omp parallel for
  for (int i=0;i<N+1; i++)
    C.row[i] = C.ptr[i];
  
  free( C.ptr );
  A.free();
  return C;
}

void cgnr( D_Matrix A ) {

  D_Matrix M;
  DD_Vector x;
  DD_Vector b;
  DD_Vector bb;
  DD_Vector r;
  DD_Vector p;
  DD_Vector q;
  DD_Vector vec_tmp;
  DD_Vector init;

  x.malloc( A.N );
  b.malloc( A.N );
  bb.malloc( A.N );
  r.malloc( A.N );
  p.malloc( A.N );
  q.malloc( A.N );
  vec_tmp.malloc( A.N );
  init.malloc( A.N );

  x.broadcast( 0.0 );
  b.broadcast( 0.0 );
  init.broadcast( 1.0 );

  /*
  for( int i=0; i<A.N; i++ )
    init.hi[i] = pow( -1.0,i );
  */
  /*
  srand((unsigned)time(NULL));
  for (int i = 0; i < A.N; i++) 
    init.hi[i] = (double)rand()/(double)RAND_MAX;  
  */

  //D_Scalar
  D_Scalar resid;
  DD_Scalar r_nrm2;
  DD_Scalar b_nrm2;

  //DD_Scalar
  DD_Scalar alpha;
  DD_Scalar beta;
  DD_Scalar a_scl;
  DD_Scalar b_scl;
  DD_Scalar c_scl;

  M.malloc( A.N );
  M = SpMM( A );

  SpMV( M,x,vec_tmp );
  SpMV( M,init,b );

  DD_AVX_xpay( b,( DD_Scalar )( -1.0 ),vec_tmp );
  r.copy( vec_tmp );
  vec_tmp.free();
  p.copy( r );

  DD_AVX_nrm2( r, &r_nrm2 );
  DD_AVX_nrm2( b, &b_nrm2 );

  resid = r_nrm2 / b_nrm2;
  int count = 1;
  double start,end;

  start = omp_get_wtime();
  while ( count < MAXITR ) {
    resid.print();
    SpMV( M,p,q );
    
    a_scl.dot( r,r );
    b_scl.dot( p,q );
    alpha = a_scl / b_scl ;

    DD_AVX_axpy( alpha,p,x );
    DD_AVX_axpy( -alpha,q,r );
    
    c_scl.dot( r,r );

    if ( resid < TOL )
      break;
    beta = c_scl / a_scl ;
    DD_AVX_xpay( r,beta,p );
    DD_AVX_nrm2( r,&r_nrm2 );
    resid = r_nrm2 / b_nrm2 ;
    count++;
  }
  end = omp_get_wtime();

  if ( SOLVE ) {
    cout << "----------------------------------解--------------------------------------\n" << endl;
    x.print_all();
    cout << "--------------------------------------------------------------------------"  << endl;
  }
  cout << "\n" << "%反復回数は" << count-1 << "\n" << endl;
  cout << "ループの時間は" << end - start << endl;

  M.free();
  x.free();
  r.free();
  p.free();
  q.free();
}

int main ( int argc, char* argv[] )
{
  char filename[256] = {'\0'};
  sprintf( filename,"%s",argv[1] );
  D_Matrix A;
  A.input( filename );
  double start,end;
  
  start = omp_get_wtime();
  cgnr( A );
  end = omp_get_wtime();
  
  cout << "\n" << "%計測時間は";
  cout << end - start << endl;
}
