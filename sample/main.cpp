#include <stdio.h>
#include <iostream>
#include <math.h>
#include <DD-AVX.hpp>

//N=12までは正常に作動します
#define N 12

int main(){

  D_Vector x;
  D_Vector y;
  DD_Vector u;
  DD_Vector v;

  D_Scalar alpha;
  DD_Scalar beta;
  DD_Scalar div;
  DD_Scalar a;
  DD_Scalar b;
  DD_Scalar c;
  DD_Scalar d;
  
  x.malloc( N );
  y.malloc( N );
  u.malloc( N );
  v.malloc( N );

  x.broadcast(3);//D_vector
  y.broadcast(3);//D_vector

  u.broadcast(3);//DD_vector
  v.broadcast(3);//DD_vector

  alpha.dot(x,y);//D_Scalar
  beta.dot(u,v); //DD_Scalar

  std::cout << "Dot" << std::endl;
  std::cout << "D:";
  alpha.print();
  std::cout << "DD: ";
  beta.print();
  std::cout << "-------------------------------------------------------------------" << std::endl;
  std::cout << "nrm2" << std::endl;
  DD_AVX_nrm2(x,&alpha);
  DD_AVX_nrm2(v,&beta);
  std::cout << "D: ";
  alpha.print();
  std::cout << "DD: ";
  beta.print();


  return 0;
  
}
