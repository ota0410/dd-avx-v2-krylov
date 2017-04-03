/* Copyright (C) T.Hishinuma, All rights reserved.*/
/* 2016/12/06*/
#include<DD-AVX.hpp>

int main(int argc, char* argv[])
{
   char file_name[256] = {'\0'};
   sprintf(file_name, "%s",argv[1] );

   DD_AVX_version();

   D_Scalar da = 1.0, db = 2.0, dc = 3.0;
   DD_Scalar a = 4, b = 5, c = 6;
   
   D_Matrix A;
   A.format = 1;
   A.input(file_name);
   
   DD_Vector x,y;
   x.malloc(A.N);
   x.broadcast(1);
   y.malloc(A.N);
   y.broadcast(0);

   DD_AVX_axpy(a, x, y);
   DD_AVX_SpMV(A, x, y);
   DD_AVX_TSpMV(A, x, y);
   a.dot(x,y);
   da.dot(x,y);
//   y.print_all();
   a.print();

   a = (double)1;
   a = da = 1.0;

   b = a * da + b / dc + (a + b)/da + (da-db);
   //b.print();

}
