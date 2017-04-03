#include "../include/DD-AVX.hpp"
#include <math.h>
#include <time.h>
#include <iostream>
/*************************************************
DDの限界はN=13(程度),Dの限界はN=???
gammaの値とサイズNの間の関係がよく分からん
**************************************************/

using std::cout;
using std::endl;

int main( int argc, char* argv[] ){

  char filename1[256] = {'\0'};/*普通の行列*/
  char filename2[256] = {'\0'};/*転置行列*/

  sprintf(filename1, "%s",argv[1] );
  sprintf(filename2, "%s",argv[2] );
  D_Matrix A,tA;
  A.input(filename1);
  tA.input(filename2);
  A.format = 1;/*CRS形式で保存*/
  tA.format = 1;

  #define N A.N

  D_Vector init,b,x,y,r,rr,p,pp,q,qq;
  init.malloc( N );
  b.malloc( N );
  r.malloc( N );
  rr.malloc( N );
  x.malloc( N );
  y.malloc( N );
  p.malloc( N );
  pp.malloc( N );
  q.malloc( N );
  qq.malloc( N );

  D_Scalar num,param,alf,bet,a_scl,b_scl,c_scl,minus;
  D_Scalar dist0,dist1,rate;
  //srand( ( unsigned )time( NULL ) );
  num = ( double )( rand() % 10 + 2 );/*値を適当に*/
  x.broadcast( num );
  init.broadcast( 1 );
  b.broadcast( 0 );
  r.broadcast( 0 );
  rr.broadcast( 0 );
  p.broadcast( 0 );
  pp.broadcast( 0 );
  q.broadcast( 0 );
  qq.broadcast( 0 );
  y.broadcast( 0 );

  DD_AVX_SpMV( A,x,y );
  DD_AVX_SpMV( A,init,b );/* b = A * init */
  init.free();/*メモリリーク対策*/

  DD_AVX_xpay( b,(D_Scalar)(-1),y );
  r.copy( y );
  y.free();/*メモリリーク対策*/

  do { /*内積が非零になるまで*/
    num = ( rand() % 5 + 2 );
    rr.broadcast( num );
    param.dot(r,rr);
  } while ( param == 0 );

  p.copy( r );
  pp.copy( rr );
  DD_AVX_nrm2( r,&dist0 );
  DD_AVX_nrm2( b,&dist1 );
 
  int count = 1;
  rate = ( dist0 / dist1 );
  cout<< "---------------------------------相対残差--------------------------------------\n" << endl;
  while( ( dist0 / dist1 ) > pow( 10,-12 ) ){
   
    cout << count << ":";
    rate.print();/*相対残差*/

    DD_AVX_SpMV( A,p,q );/*q = A*p */
    DD_AVX_SpMV( tA,pp,qq );/* 行列の転置処理 q_k = A_t * p_k */
    
    alf = ( a_scl.dot( rr,r ) / b_scl.dot( pp,q ) );
    param = a_scl;

    DD_AVX_axpy( alf,p,x );
    DD_AVX_axpy( -alf,q,r );
    DD_AVX_axpy( -alf,qq,rr );

    q.broadcast( 0 );
    qq.broadcast( 0 );

    bet = ( c_scl.dot( rr,r ) / param );

    DD_AVX_xpay(r,bet,p);
    DD_AVX_xpay(rr,bet,pp);
    
    DD_AVX_nrm2( r,&dist0 );
    DD_AVX_nrm2( b,&dist1 );

    rate = ( dist0 / dist1 );
    count++;
  }
  cout << "-------------------------------------------------------------------------------\n" << endl;
  /*
 for(int i=0;i<N;i++)
    x.print( i );
  */
 A.free();/*free関数はやむなく自作*/
 tA.free();
 
 b.free();
 x.free();
 p.free();
 pp.free();
 r.free();
 rr.free();
 q.free();
 qq.free();
 return 0;
}
