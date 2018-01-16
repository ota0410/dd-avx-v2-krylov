#include <DD-AVX.hpp>
#include <iostream>
#include <string.h>
#include <time.h>
#include <omp.h>

using std::cout;
using std::endl;

#define N A.N
#define TOL 1.0e-10
#define SOLVE 0
#define maxitr A.N
#define OMEGA 1.22
#define PRECOND 0

D_Scalar dot( D_Vector a, D_Vector b )
{
  double scl = 0.0;
  D_Scalar dot;
  int n = a.getsize();
  //#pragma omp parallel for reduction(+:scl)
  for( int i = 0;i < n; i++ )
    scl = scl + a.hi[i] * b.hi[i];
  dot.hi = scl;
  return dot;
}

void cgne( D_Matrix A )
{
  D_Vector x;
  D_Vector b;
  D_Vector r;
  D_Vector rr;
  D_Vector p;
  D_Vector pp;
  D_Vector q;
  D_Vector vec_tmp;
  D_Vector init;

  x.malloc( N );
  b.malloc( N );
  r.malloc( N );
  rr.malloc( N ); 
  p.malloc( N );
  q.malloc( N );
  vec_tmp.malloc( N );
  init.malloc( N );

  x.broadcast( 0 );
  /*  
  srand((unsigned)time(NULL));
  for (int i = 0; i < N; i++) 
    init.hi[i] = (double)rand()/(double)RAND_MAX;
  
  init.broadcast(1.0);
  */
  
  for (int i = 0; i < N; i++) 
    init.hi[i] = pow(-1.0,i);
  
  b.broadcast( 0 );
  
  D_Scalar resid;
  D_Scalar r_nrm2;
  D_Scalar b_nrm2;
  D_Scalar alpha;
  D_Scalar beta;
  D_Scalar a_scl;
  D_Scalar b_scl;
  D_Scalar c_scl;

  DD_AVX_SpMV( A,init,b );
  /****** Forward NE-SOR sweep*********/
  if ( PRECOND ){
    D_Scalar bet,del,d,row;
    D_Vector e,f;
    e.malloc( N );
    f.malloc( N );
    for(int i = 0; i < N; i++){
      bet = b.hi[i];
      e.hi[i] = 1.0;
      
      DD_AVX_TSpMV( A,e,f );
      
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
  DD_AVX_SpMV( A,x,vec_tmp );
  DD_AVX_xpay( b,( D_Scalar )(-1.0),vec_tmp );
  r.copy( vec_tmp );

  DD_AVX_TSpMV( A,r,p );
  DD_AVX_nrm2( r, &r_nrm2 );
  DD_AVX_nrm2( b, &b_nrm2 );
  resid = r_nrm2 / b_nrm2;
  /*
  D_Vector error;
  D_Scalar serror;
  error.malloc( N );
  DD_AVX_axpyz( ( D_Scalar )(-1.0), init, x, error );
  DD_AVX_nrm2( error, &serror );
  serror.print();
  */
  int count = 1;
  while ( count < N ) {

    resid.print();
    //serror.print();
    DD_AVX_SpMV( A,p,q );

    a_scl = dot( r,r );
    b_scl = dot( p,p );
    
    alpha = a_scl / b_scl;

    DD_AVX_axpy( alpha,p,x );
    DD_AVX_axpy( -alpha,q,r );    

    c_scl = dot( r,r );

    beta = c_scl / a_scl;
    if ( resid < TOL ){
      //      resid.print();
      break;
    }
    
    DD_AVX_TSpMV( A,r,rr );
    DD_AVX_xpay( rr,beta,p );
    
    DD_AVX_nrm2( r,&r_nrm2 ); 
    resid = r_nrm2 / b_nrm2;
    
    // DD_AVX_axpyz( ( D_Scalar )(-1.0), init, x, error );
    //  DD_AVX_nrm2( error, &serror );
    count++;
  }
  
  if ( SOLVE ){
    cout << "----------------------------------解--------------------------------------\n" << endl;
    x.print_all();
    cout << "---------------------------------------------------------------------------"  << endl;
  }
  cout << count  << ",";
  //  resid.print();
  
  A.free();
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
  cgne( A );
  end = omp_get_wtime();
    
  cout << "計測時間は";
  cout << end - start << endl;
  
  return 0;
}
