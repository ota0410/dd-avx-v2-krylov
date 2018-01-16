#include <DD-AVX.hpp>
#include <iostream>
#include <time.h>
#include <omp.h>

using std::cout;
using std::endl;

#define maxitr M.N
#define TOL 1.0e-10
#define SOLVE 0

D_Matrix SpMM( D_Matrix A )
{
  D_Matrix C;
  C.malloc( A.N );
  D_Vector x,y;
  x.malloc( A.N );
  y.malloc( A.N );
  C.N = x.getsize();

  int count = 0;
  for (int k=0; k<A.N; k++) {
    x.broadcast( 0.0 );
    for (int i=0; i<A.N; i++) {
      for (int j=A.row[i]; j<A.row[i+1]; j++) {
	if ( A.col[j] == k )
	  x.hi[i] = A.val[j];
      }
    }
    DD_AVX_TSpMV( A,x,y );
    for (int h=0; h<A.N; h++) {
      if ( y.hi[h] != 0 ){
	C.val[count] = y.hi[h];
	C.row[count] = h;
	count++;
      }
    }
    C.col[k+1] = count;
  }
  C.nnz = count;
  C.col[C.N] = C.nnz+1;
  
  for (int i=0;i<C.N+1; i++)
    C.ptr[i] = C.col[i];
  for (int i=0; i<C.nnz; i++)
    C.col[i] = C.row[i];
  for (int i=0;i<C.N+1; i++)
    C.row[i] = C.ptr[i];

  free( C.ptr );
  //  C.print_all();
  return C;
}

D_Scalar dot( D_Vector x, D_Vector y )
{
  D_Scalar scl = 0.0;
  D_Scalar tmp = 0.0;
  D_Scalar yy = 0.0;
  D_Scalar e = 0.0;
  int size = x.getsize();

  for( int i = 0;i < size; i++ )
    scl = scl + x.hi[i] * y.hi[i];
  
  return scl;
}

void cgnr( D_Matrix A )
{
  D_Matrix M;

  D_Vector x;
  D_Vector b;
  D_Vector bb;
  D_Vector r;
  D_Vector p;
  D_Vector q;
  D_Vector vec_tmp;
  D_Vector init;

  x.malloc( A.N );
  b.malloc( A.N );
  bb.malloc( A.N );
  r.malloc( A.N );
  p.malloc( A.N );
  q.malloc( A.N );
  vec_tmp.malloc( A.N );
  init.malloc( A.N );

  x.broadcast( 0 );
  b.broadcast( 0 );
  init.broadcast( 1 );
  /*
  for( int i=0; i<A.N; i++ )
    init.hi[i] = pow( -1.0,i );
  */
  /*
  srand((unsigned)time(NULL));
  for (int i = 0; i < A.N; i++) 
    init.hi[i] = (double)rand()/(double)RAND_MAX;  
  */
  D_Scalar resid;
  D_Scalar r_nrm2;
  D_Scalar b_nrm2;
  D_Scalar alpha;
  D_Scalar beta;
  D_Scalar a_scl;
  D_Scalar b_scl;
  D_Scalar c_scl;

  M.malloc( A.N );
  /*
  DD_AVX_SpMV( A,init,bb );
  init.free();

  DD_AVX_TSpMV( A,bb,b );
  bb.free();
  */
  M = SpMM( A );

  DD_AVX_SpMV( M,x,vec_tmp );
  DD_AVX_SpMV( A,init,b );
  DD_AVX_TSpMV( A,b,bb );
  DD_AVX_xpay( bb,(D_Scalar)(-1.0),vec_tmp );
  r.copy( vec_tmp );
  vec_tmp.free();
  p.copy( r );

  DD_AVX_nrm2( r, &r_nrm2 );
  DD_AVX_nrm2( bb, &b_nrm2 );

  resid = r_nrm2 / b_nrm2;
  int count = 1;
  //  cout << "---------------------------------相対残差----------------------------------"  << endl;

  double start,end;
  start = omp_get_wtime();
  while ( resid > TOL && count < maxitr ) {
    resid.print();
    DD_AVX_SpMV( M,p,q );
    
    a_scl = dot( r,r );
    b_scl = dot( p,q );

    alpha = a_scl / b_scl ;
    
    DD_AVX_axpy( alpha,p,x );
    DD_AVX_axpy( -alpha,q,r );
    
    c_scl = dot( r,r );
    beta = c_scl / a_scl;

    DD_AVX_xpay( r,beta,p );
    DD_AVX_nrm2( r,&r_nrm2 );
    resid = r_nrm2 / b_nrm2 ;

    count++;
  }

  end = omp_get_wtime();
  resid.print();
  /*
  cout << "ループの計測時間は";
  cout << end - start << endl; 
  */
  if ( SOLVE ) {
    cout << "----------------------------------解--------------------------------------\n" << endl;
    x.print_all();
    cout << "---------------------------------------------------------------------------"  << endl;
  }

  cout << "%反復回数は" << count -1 << "\n" << endl;
  M.free();
  b.free();
  x.free();
  r.free();
  p.free();
  q.free();
}

int main ( int argc, char* argv[] )
{
  char filename[256] = {'\0'};
  sprintf( filename,"%s",argv[1] );
  D_Matrix A,M;
  A.input( filename );
  double start,end;
  start = omp_get_wtime();
  cgnr( A );
  end = omp_get_wtime();
  cout << "%実行時間は" << end - start << endl; 
  return 0;
}
