#include <DD-AVX.hpp>
#include <iostream>
#include <ctime>
#include <omp.h>
#include <cstring>

using std::cout;
using std::endl;

#define TOL 1.0e-10
#define MAXITR 5000
#define SOLVE 0
#define PRECOND 0
#define OMEGA 1.22

DD_Vector SpMV( D_Matrix A,DD_Vector x, DD_Vector y ) 
{
  DD_Scalar k,temp;
  int size = x.getsize();

  y.broadcast( 0.0 );
  for (int i = 0; i < size; i++) {
    temp = 0.0;
    for (int j = A.row[i]; j < A.row[i+1]; j++){
      k = x.getelm( A.col[j] );
      temp = temp + ( DD_Scalar )A.val[j] * k;
    }
    y.chgelm( temp,i );
  }
  return y;
}

DD_Vector TSpMV( D_Matrix A,DD_Vector x, DD_Vector y )
{
  DD_Scalar k,temp;
  y.broadcast( 0.0 );
  int size = x.getsize();
  for (int i = 0 ; i < size; i++) {
    k = x.getelm(i);      
    for (int j = A.row[i]; j < A.row[i+1]; j++){
      temp = y.getelm( A.col[j] ) + ( DD_Scalar )A.val[j] * k;
      y.chgelm(  temp, A.col[j] );
    }
  }
  return y;
}

void cgne( D_Matrix A )
{
  DD_Vector x;
  DD_Vector b;
  DD_Vector r;
  DD_Vector rr;
  DD_Vector p;
  DD_Vector q;
  DD_Vector vtmp;
  DD_Vector init;

  x.malloc( A.N );
  b.malloc( A.N );
  r.malloc( A.N );
  rr.malloc( A.N ); 
  p.malloc( A.N );
  q.malloc( A.N );
  vtmp.malloc( A.N );
  init.malloc( A.N );
      
  x.broadcast( 0 );
  //  init.broadcast( 1 );
  init.broadcast( 0.0 );
  
  srand((unsigned)time(NULL));
  for (int i = 0; i < A.N; i++)
    init.hi[i] = (double)rand()/(double)RAND_MAX;
  /*
  for (int i = 0; i < A.N; i++){
    init.hi[i] = pow(-1,i);
    init.lo[i] = 0.0;
  }
  */
  b.broadcast( 0 );
  r.broadcast( 0 );
  
  D_Scalar resid;
  D_Scalar r_nrm2;
  D_Scalar b_nrm2;
  
  DD_Scalar alpha;
  DD_Scalar beta;
  DD_Scalar a_scl;
  DD_Scalar b_scl;
  DD_Scalar c_scl;

  SpMV( A,init,b );
  //init.free();
  /****** Forward NE-SOR sweep*********/
  if ( PRECOND ){
    DD_Scalar bet,del,d,row;
    DD_Vector e,f;
    e.malloc( A.N );
    f.malloc( A.N );
    for(int i = 0; i < A.N; i++){
      bet = b.getelm(i);
      e.hi[i] = 1.0;
      
      TSpMV( A,e,f );
      
      d.dot( f,x );
      row.dot( f,f );
      
      del = OMEGA * ( bet - d ) / row;
      
      DD_AVX_axpy( del,f,x );
      e.hi[i] = 0.0;
    }
    e.free();
    f.free();
  }
  /***********************************/

  SpMV( A,x,vtmp );

  DD_AVX_xpay( b,( DD_Scalar )(-1.0),vtmp );
  r.copy( vtmp );
  vtmp.free();

  TSpMV( A,r,p );
  DD_AVX_nrm2( r, &r_nrm2 );
  DD_AVX_nrm2( b, &b_nrm2 );

  resid = r_nrm2 / b_nrm2;

  int count = 1;
  while ( resid > TOL && count < MAXITR ) {
    resid.print();
    SpMV( A,p,q );
    
    a_scl.dot( r,r );
    b_scl.dot( p,p );
 
    alpha = a_scl / b_scl ;
    
    DD_AVX_axpy( alpha,p,x );
    DD_AVX_axpy( -alpha,q,r );

    c_scl.dot( r,r );
    beta = c_scl / a_scl ;

    TSpMV( A,r,rr );

    DD_AVX_xpay( rr,beta,p );
    DD_AVX_nrm2( r,&r_nrm2 ); 

    resid = r_nrm2 / b_nrm2 ;
    count++;
  }
  if ( SOLVE ) {
    cout << "----------------------------------解--------------------------------------\n" << endl;
    x.print_all();
    cout << "---------------------------------------------------------------------------"  << endl;
  }
  
  cout << count << ",";
  resid.print();
  /*
  if ( A.N - count > 0 )
    cout << "%収束しました" << endl;
  else
    cout << "%収束しませんでした" << endl;
  */
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
  
  A.free();
  b.free();
  x.free();
  r.free();
  rr.free();
  p.free();
  q.free();
}

int main ( int argc, char* argv[] )
{
  char filename[256] = {'\0'};
  sprintf( filename,"%s",argv[1] );

  D_Matrix A;
  char *matrix;
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
  cgne( A );
  end = omp_get_wtime();
  cout << end - start << endl; 
  return 0;
}
