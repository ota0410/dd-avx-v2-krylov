#include <stdio.h>
#include <iostream>
#include <math.h>
#include <DD-AVX.hpp>

int main( int argc, char *argv[]){

  char filename1[256] = {'\0'};/*元の行列*/
  sprintf( filename1, "%s",argv[1] );
 
  D_Matrix A;
  A.input( filename1 );
  //  A.print_all();
  printf("--------------col----------------\n");
  for(int i=0; i < A.nnz; i++)
    printf("%d\n",A.col[i]);
  printf("--------------row----------------\n");
  for(int i=0; i < A.N+1; i++)
    printf("%d\n",A.row[i]);
  printf("--------------val----------------\n");
  for(int i=0; i < A.nnz; i++)
    printf("%f\n",A.val[i]);
  printf("-----------------------------------\n");
  A.print_all();
  return 0;
  
}
