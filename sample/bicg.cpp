#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../include/DD-AVX.hpp"
#include <iostream>

#define N 100
#define GAMMA 0.8

DD_Scalar distanceRate( DD_Vector *vec0,DD_Vector *vec1 ){
  DD_Scalar dist_0 = 0, dist_1 = 0;
  /*
  for(int i=0;i<N;i++){
    square_dist_0 += (vec0[i] * vec0[i]);
    square_dist_1 += (vec1[i] * vec1[i]);
  }
  DD_Scalar dist_0 = sqrt(square_dist_0);
  DD_Scalar dist_1 = sqrt(square_dist_1);*/
  DD_AVX_nrm2( vec0,dist_0);
  DD_AVX_nrm2( vec1,dist_1);
  DD_Scalar rate = (dist_0/dist_1);
  //printf("%.12Lf\n",rate);/*相対残差*/
  // printf("%E\n",rate);/*相対残差*/
  //  std::cout << rate << std::endl;
  rate.print();
  return rate;
}

/*void init(DD_Vector *idata,DD_Scalar num){ 
  for(int i = 0;i < N;i++){
    idata[i] = num;
  }
}
*****/

/*DD_Scalar dotProduct(double *vec0,double *vec1){
  DD_Scalar sum = 0;
  for(int i=0;i<N;i++){
    sum += vec0[i] * vec1[i];
  }
  return sum;
}
*/

void bicg(DD_Vector *val,int *col_ind ,int *row_ptr,DD_Vector *vec){
  DD_Vector r,rr,p,pp,q,qq;
  r.malloc( N );
  rr.malloc( N );
  p.malloc( N );
  pp.malloc( N );
  q.malloc( N );
  qq.malloc( N );
  srand( (unsigned)time(NULL) );
  DD_Scalar num = ( double )( rand() % 100 + 2 ); /*初期値は適当*/
  DD_Vector vec_x;
  vec_x.malloc( N );
  vec_x.broadcast( num );

  /*  for(int i=0;i<N;i++){
    DD_Scalar tmp = 0;
    int count = 0;
    for(int j=row_ptr[i];j<row_ptr[i+1];j++){
      tmp += val[j] * vec_x[count];
      count++;
    }
    r[i] = vec[i] - tmp;/matrix_multiple/
  */
    rr[i] = ( DD_Scalar )( rand() % 100 + 2 );  /* choose r_*_0 s.t. (r*,r) != 0(変える?) */
    pp[i] = rr[i];
    p[i] = r[i];
  }  

  spMV
  int num2 = 1; 
  printf("%d: ",num2);
  num2++;
  DD_Scalar d_rate = distanceRate(r,vec);
  while( d_rate > pow(10,-12) ){
    printf("%d: ",num2);
    int count = 0; /* q_k = A * p_k 計算過程 */
    for(int i=0;i<N;i++){
      DD_Scalar tmp = 0;
      for(int j=0;j<N;j++){
	if( col_ind[count]  == j){
	  tmp += val[count] * p[j];
	  count++;
	}
      }
      q[i] = tmp;
    }
    /* qq_k = A_t * pp_k */
    for(int i=0;i<N;i++){
      DD_Scalar tmp = 0;
      count = 0;
      for(int j=0;j<3*N-2;j++){
	if( j >= row_ptr[count+1])
	    count++;
	if( col_ind[j] == i)
	  tmp += val[j] * pp[count];
      }
      qq[i] = tmp;
    }
    
    DD_Scalar tmp_r,tmp_p;
    tmp_r.dot(r,rr); 
    tmp_p.dot(pp,q);
    DD_Scalar scl_alf = tmp_r / tmp_p;

    DD_Scalar tmp_r2 = 0;
    DD_AVX_axpy(scl_alf,p,vec_x);
    DD_AVX_axpy(-scl_alf,q,r);
    DD_AVX_axpy(-scl_alf,qq,rr);
    tmp_r2.dot(r,rr);
    /*    for(int i=0;i<N;i++){
      /* vec_x[i] = vec_x[i] + scl_alf * p[i];*/
      /*r[i] = r[i] - scl_alf * q[i];*/
      /*rr[i] = rr[i] - scl_alf * qq[i];*/
    //   tmp_r2 += r[i] * rr[i];
      //}
    DD_Scalar scl_bet = tmp_r2/tmp_r;

    /*    for(int i=0;i<N;i++){
      p[i] = r[i] + scl_bet * p[i];
      pp[i] = rr[i] + scl_bet * pp[i];
      }*/
    DD_AVX_axpy(scl_bet,pp,rr);
    d_rate = distanceRate(r,vec);
    num2++;
  }
  /*
  printf("----------------------result----------------------\n");
  for(int i=0;i<N;i++){
    printf("%f\n",vec_x[i]);
  }
  */
}

int main(){

  DD_Vector b,vector,val;

  b.malloc(sizeof(double)*N);/*ベクトルb*/
  vector.malloc(sizeof(double)*N); /*初期のベクトル*/
  /*以下CRS形式の疎行列表現*/
  val.malloc(sizeof(double)*(3*N-2));
  /*非零要素を格納するための配列*/
  int col_ind = (int*)malloc(sizeof(int)*(3*N-2));
  /*非零要素の列番号を格納*/
  int row_ptr = (int*)malloc(sizeof(int)*(N+1));
  /*先頭の非零要素の位置を示す*/
  //init( vector,1 );
  //init( b,1 );
  vector.broadcast(DD_Scalar,1);
  b.broadcast(DD_Scalar,1);

  int i=0,j=0;
  for(i=0;i<3*N-3;i+=3){
    val[i] = 2;
    val[i+1] = 1;
    val[i+2] = GAMMA;/*gammaの値(変更して実験)*/
  }
  val[i] = 2;

  col_ind[0] = 0;
  col_ind[1] = 1;
  int param = 0;
  for(i=2;i<3*(N-2);i+=3){ 
    col_ind[i] = param;
    col_ind[i+1] = param + 1;
    col_ind[i+2] =  param + 2;
    param += 1;
  }
  col_ind[3*N-4] = param;
  col_ind[3*N-3] = param + 1;

  row_ptr[0] = 0;
  for(i=1;i<N;i++){
    row_ptr[i] = 3*i-1;
  }
  row_ptr[i] = 3*N-1; /*最後は非零要素数+1の値を格納する*/

  for(i=0;i<N;i++){
    DD_Scalar tmp = 0;
    int count = 0;
    while(count < row_ptr[i+1] - row_ptr[i]){
      tmp += val[count] * vector[i];
      count++;
    }
    b[i] = tmp;
  }
  b[N-1] = val[3*N-3] * vector[N-1] + val[3*N-4] * vector[N-2];

  free( vector );
  bicg( val,col_ind,row_ptr,b );
  free( val );
  free( col_ind );
  free( row_ptr );
  return 0;
}
