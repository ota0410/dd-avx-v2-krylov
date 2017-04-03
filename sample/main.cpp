/* Copyright (C) T.Hishinuma, All rights reserved.*/
/* 2016/12/06*/
#include "../include/DD-AVX.hpp"
#include <iostream>

int main(int argc, char* argv[])
{
   char filename1[256] = {'\0'};
   //   char filename2[256] = {'\0'};

   sprintf(filename1, "%s",argv[1] );
   //  sprintf(filename2, "%s",argv[2] );

   DD_AVX_version();

   D_Scalar da = 1.0, db = 2.0, dc = 3.0;
   DD_Scalar a = 4, b = 5, c = 6;
   
   D_Matrix A/*B*/;
   A.format = 1;/*crs*/
   //  B.format = 1;
   A.input(filename1);
   //  B.input(filename2);
   
   DD_Vector x,y;/*宣言のちmalloc*/

   x.malloc(A.N); /*vectorを行列の大きさ分領域確保*/
   x.broadcast(1.41421356); /*全要素を1で初期化*/
   y.malloc(A.N);
   y.broadcast(0);

   DD_AVX_axpy(a, x, y);
   y.print_all();
   std::cout << "------------------------------------------------" << std::endl;
   DD_AVX_SpMV(A, x ,y);
   // DD_AVX_SpMV(B, x ,y);
   //   DD_AVX_scale(a,x);
   //   a.dot(x,y);
   //   da.dot(x,y);
   //   std::cout << da.dot(x,y) << std::endl;
   //   A.getsize();
   //   A.print_all();
   y.print_all();
   std::cout << "------------------------------------------------" << std::endl;
   // a.print();
   // b.print();
   x.print_all();
   std::cout << "-----------------------------------------------" << std::endl;
   y.print_all();
   std::cout << "-----------------------------------------------" << std::endl;
   A.print_all();/*行->列->値の順番*/
   // std::cout << "-----------------------------------------------" << std::endl;
   //A.getsize();
   //   std::cout << A.N << std::endl;
   da.hello( da );
   //   std::cout << da << std::endl;
   a = (double)1;
   a = da = 1.0;

   b = a * da + b / dc + (a + b)/da + (da-db);
   //b.print();
}
