#include <DD-AVX.hpp>
#include <iostream>
#include <string.h>
#include <omp.h>

using std::cout;
using std::endl;

#define TOL 1.0e-10
#define MAXITR A.N
#define SOLVE 0

DD_Vector TSpMV( D_Matrix A,DD_Vector x, DD_Vector y )
{
  int i,j;
  DD_Scalar k;
  DD_Scalar temp;
  y.broadcast( 0.0 );
  int size = x.getsize();
  for (i = 0 ; i < size; i++) {
    for (j = A.row[i]; j < A.row[i+1]; j++){
      temp = y.getelm( A.col[j] ) + A.val[j] * x.getelm( i );
      y.chgelm(  temp, A.col[j] );
    }
  }
  return y;
}

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

void bicg( D_Matrix A ) {  

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

  init.malloc( A.N );
  b.malloc( A.N );
  r.malloc( A.N );
  rr.malloc( A.N );
  x.malloc( A.N );
  y.malloc( A.N );
  p.malloc( A.N );
  pp.malloc( A.N );
  q.malloc( A.N );
  qq.malloc( A.N );

  //D_Scalar
  D_Scalar resid;
  D_Scalar b_nrm2;
  D_Scalar r_nrm2;
 
  //DD_Scalar
  DD_Scalar a_scl;
  DD_Scalar b_scl;
  DD_Scalar c_scl;
  DD_Scalar alpha;
  DD_Scalar beta;

  x.broadcast( 0.0 );
  
  //  init.broadcast( 1.0 );
  /* 
  srand((unsigned)time(NULL));
  for (int i = 0; i < A.N; i++) 
    init.hi[i] = (double)rand()/(double)RAND_MAX;
  */
  //  init.broadcast( 1.0 );
  //  init.output_plane("rand_vec.mtx");
         
  for (int i = 0; i < A.N; i++){
    init.hi[i] = 1.0 * pow(-1.0,i);
    init.lo[i] = 0.0;
  }
  
  
  r.broadcast( 0.0 );
  p.broadcast( 0.0 );
  pp.broadcast( 0.0 );
  q.broadcast( 0.0 );
  qq.broadcast( 0.0 );

  SpMV( A,x,y );
  SpMV( A,init,b );

  DD_AVX_xpay( b,( DD_Scalar )(-1.0),y );
  r.copy( y );
  rr.copy( r );
  p.copy( r );
  pp.copy( rr );

  DD_AVX_nrm2( r,&r_nrm2 );
  DD_AVX_nrm2( b,&b_nrm2 );
  b_nrm2 = (DD_Scalar)1.0 / b_nrm2;
  resid = r_nrm2 * b_nrm2 ;

  int count = 1;
  while ( resid > TOL ) {
    resid.print();

    SpMV( A,p,q );
    TSpMV( A,pp,qq ); 

    a_scl.dot( rr,r );
    b_scl.dot( pp,q );
    
    if ( count > MAXITR ) {
      /*
      cout << "\n%BiCGは与えられた許容誤差に収束しませんでした" << endl;
      cout << "%最後の相対残差は\n" << endl;*/
      break;
    }
    
    alpha = a_scl / b_scl ;

    DD_AVX_axpy( alpha,p,x );
    DD_AVX_axpy( -alpha,q,r );
    DD_AVX_axpy( -alpha,qq,rr );

    c_scl.dot( rr,r );
    beta = c_scl / a_scl ;

    DD_AVX_xpay( r,beta,p );
    DD_AVX_xpay( rr,beta,pp );

    DD_AVX_nrm2( r,&r_nrm2 );
    resid = r_nrm2 * b_nrm2 ;

    count++;
  }
  
  if ( SOLVE ) {
    cout << "----------------------------------解ベクトル------------------------------------\n" << endl;
    x.print_all();
    cout << "--------------------------------------------------------------------------------" << endl;
  }
  
  // cout << "%次元数は" << A.N << endl;
  cout << count << ",";
  resid.print();
  DD_Vector error;
  DD_Scalar serror;
  DD_Scalar sinit;
  error.malloc( A.N );
  DD_AVX_axpyz( ( DD_Scalar )(-1.0), init, x, error );
  DD_AVX_nrm2( error, &serror );
  DD_AVX_nrm2( init, &sinit );
  serror = serror / sinit;
  
  cout << "error = ";
  serror.print();
  
  /*
  if ( ( A.N - count ) > 0 )
    cout << "%収束しました" << endl;
  else 
    cout << "%収束しませんでした"  << endl;
  */
  A.free();
  b.free();
  x.free();
  p.free();
  pp.free();
  r.free();
  rr.free();
  q.free();
  qq.free();  
}

int main( int argc, char* argv[] )
{

  char filename[256] = {'\0'};
  char *matrix;
  sprintf( filename, "%s",argv[1] );

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
  bicg( A );
  end = omp_get_wtime();
  cout << end - start << endl;
  return 0;
}
