#include <stdio.h>
#include <iostream>
#include <math.h>
#include <DD-AVX.hpp>

int main( int argc, char *argv[]){

  char filename1[256] = {'\0'};/*元の行列*/
  sprintf( filename1, "%s",argv[1] );
 
  D_Matrix A;
  A.input( filename1 );


  DD_Vector x;
  DD_Vector y;
  DD_Vector u;
  DD_Vector v;

  D_Scalar alpha;
  DD_Scalar beta = 5.0;
  DD_Scalar temp2,bb;
  
  x.malloc( A.N );
  y.malloc( A.N );
  u.malloc( A.N );
  v.malloc( A.N );

  x.broadcast( 2.0 );//D_vector
  y.broadcast( 0.0 );//D_vector


  
  return 0;
  
}
