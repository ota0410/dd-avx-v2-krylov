#include <DD-AVX.hpp>
#include <cstring>
#include <omp.h>
#include <ctime>
#include <iostream>
#include <cstdlib>

using std::cout;
using std::endl;
using std::string;

#define TOL 1.0e-10
#define MAXITR A.N
#define SOLVE 0

D_Scalar dot( D_Vector x, D_Vector y )
{
  double tmp = 0.0;
  D_Scalar scl;
  int i,size = x.getsize();
#pragma omp parallel reduction(+:tmp)
  for ( i = 0; i < size ; i++ ) 
    tmp = tmp + x.hi[i] * y.hi[i];

  scl.hi = tmp;
  return scl;
}

void bicg( D_Matrix A )
{
  D_Vector init;
  D_Vector b;
  D_Vector x;
  D_Vector y;
  D_Vector r;
  D_Vector rr;
  D_Vector p;
  D_Vector pp;
  D_Vector q;
  D_Vector qq;

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

  D_Scalar resid;
  D_Scalar r_nrm2;
  D_Scalar b_nrm2;
  D_Scalar a_scl;
  D_Scalar b_scl;
  D_Scalar c_scl;
  D_Scalar alpha;
  D_Scalar beta;
  
  x.broadcast( 0 );
  y.broadcast( 0 );
  //  init.broadcast( 1 );
  
  for (int i = 0; i < A.N; i++) 
    init.hi[i] = pow(-1,i);
  
  //  init.broadcast( 1.0 );  
  p.broadcast( 0 );
  pp.broadcast( 0 );
  q.broadcast( 0 );
  qq.broadcast( 0 );
  b.broadcast( 0 );
  r.broadcast( 0 );

  DD_AVX_SpMV( A,x,y );
  DD_AVX_SpMV( A,init,b );
  
  for (int i = 0; i < A.N; i++) 
    r.hi[i] = b.hi[i] - y.hi[i];
  
  rr.copy( r );
  p.copy( r );
  pp.copy( rr );

  DD_AVX_nrm2( r, &r_nrm2 );
  DD_AVX_nrm2( b, &b_nrm2 );

  resid = r_nrm2 * ( 1.0 / b_nrm2 );

  int count = 1;

  D_Vector error;
  D_Scalar serror;
  error.malloc( A.N );
  DD_AVX_axpyz( ( D_Scalar )(-1.0), init, x, error );
  /*
  DD_AVX_nrm2( error, &serror );
  serror.print();
  */
  while ( 1 ) {
                          
    DD_AVX_SpMV( A,p,q );
    DD_AVX_TSpMV( A,pp,qq );

    a_scl = dot( r,rr );
    b_scl = dot( pp,q );

    /*           
    a_scl.dot( r,rr );
    b_scl.dot( pp,q );
    */

    if ( count > MAXITR ) {
      /*
      cout << endl;
      cout << "BiCGは許容誤差に収束することなく停止しました" << endl;
      cout << "最後の相対残差は\n" << endl;
      */
      break;
    }

    if ( ( a_scl == 0.0 || b_scl == 0.0 ) )
      cout << "break down" << ",";
    
    alpha = a_scl * ( 1.0 / b_scl );
    resid.print();
    //serror.print();
    
    DD_AVX_axpy( alpha,p,x );
    DD_AVX_axpy( -alpha,q,r );
    DD_AVX_axpy( -alpha,qq,rr );
    
    if ( resid < TOL ) {
      break;
    }

    c_scl = dot( r,rr );
    beta =  c_scl * ( 1.0 / a_scl );

    DD_AVX_xpay( r,beta,p );
    DD_AVX_xpay( rr,beta,pp );

    //residual
    DD_AVX_nrm2( r,&r_nrm2 );
    resid = r_nrm2  / b_nrm2;
    
    //error
    DD_AVX_axpyz( ( D_Scalar )(-1.0), init, x, error );
    DD_AVX_nrm2( error, &serror );
    
    count++;
  }
  /*
  if ( SOLVE ){
    cout << "----------------------------------解ベクトル------------------------------------\n" << endl;
    x.print_all();
    cout << "--------------------------------------------------------------------------------" << endl;
  }
  */
  cout << count << ",";
  //  resid.print();
  // serror.print();
  //  cout << "%次元数は" << A.N << endl;
  /*
  if ( ( A.N-count ) > 0 ) 
    cout << "%収束しました" << endl;
  else
    cout << "%収束しませんでした" << endl;
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
  cout << "%計測時間は";
  cout << end - start << endl;

  return 0;
}
