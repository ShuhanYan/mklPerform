/*******************************************************************************
* Copyright 2001-2016 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
!  Content:
!    vsAdd  Example Program Text
!******************************************************************************/
#include "mkl.h"
#include <sys/timeb.h>
#include <stdio.h>
#include "mkl_vml.h"
#include <time.h>
#include "_rms.h"
#include "mkl_types.h"
#include "mkl_cblas.h"
#include <stdlib.h>
long long getSystemTime() {
    struct timeb t;
    ftime(&t);
    return 1000 * t.time + t.millitm;
}

int sgemvPerform(int length,int iter){
      MKL_INT         m=length, n=length, lda=length, incx=length, incy=length;
      MKL_INT         rmax=length, cmax=length;
      float           alpha=1, beta=0;
      float          *a, *x, *y;
      CBLAS_LAYOUT    layout=CblasRowMajor;
      CBLAS_TRANSPOSE trans=CblasNoTrans;
      int len_x = 1+(n-1)*abs(incx),len_y = 1+(m-1)*abs(incy);
      int i,j;
      a = (float *)mkl_calloc(rmax*cmax, sizeof(float), 64);
      x = (float *)mkl_calloc(len_x, sizeof(float), 64);
      y = (float *)mkl_calloc(len_y, sizeof(float), 64);

      //generate data
      for(i = 0; i < rmax; i++){
         for(j = 0; j < cmax; j++){
             a[i*cmax+j] = rand();
         }
      }

      for (i=0;i<len_x;i++){
         x[i]=rand();
         y[i]=0;
      }
      long long start,end;
      srand((unsigned)time(NULL));
      //warm up
      cblas_sgemv(layout, trans, m, n, alpha, a, lda, x, incx, beta, y, incy);

      //test

      start = getSystemTime();
      for(i=0;i<iter;i++){
          cblas_sgemv(layout, trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
      }
      end = getSystemTime();
      printf("%i*%ix%i*%i: %f ms\n",length,length,length,length,(double)(end-start)/iter);

      mkl_free(a);
      mkl_free(x);
      mkl_free(y);

      return 0;
}

int sgemmPerform(int length,int iter){
      //parameter
      MKL_INT         m=length, n=length, k=length;
      MKL_INT         rmax=length, cmax=length;
      MKL_INT         lda, ldb, ldc;
      // CblasRowMajor
      lda=cmax;
      ldb=cmax;
      ldc=cmax;
      int i =0, j =0;
      float           alpha=1, beta=0;
      float          *a, *b, *c;
      CBLAS_LAYOUT    layout=CblasRowMajor;
      CBLAS_TRANSPOSE transA=CblasNoTrans, transB=CblasNoTrans;

      a = (float *)mkl_calloc(rmax*cmax, sizeof( float ), 64);
      b = (float *)mkl_calloc(rmax*cmax, sizeof( float ), 64);
      c = (float *)mkl_calloc(rmax*cmax, sizeof( float ), 64);
      if( a == NULL || b == NULL || c == NULL ) {
          printf( "\n Can't allocate memory for arrays\n");
          return 1;
      }

      long long start,end;
      srand((unsigned)time(NULL));

      //generate data
      for(i = 0; i < rmax; i++){
         for(j = 0; j < cmax; j++){
             a[i*cmax+j] = rand();
             b[i*cmax+j] = rand();
             c[i*cmax+j] = 0;
         }
      }

      //warm up
      cblas_sgemm(layout, transA, transB, m, n, k, alpha,
                  a, lda, b, ldb, beta, c, ldc);

      //test
      start = getSystemTime();
      for (i=0;i<iter;i++){
          cblas_sgemm(layout, transA, transB, m, n, k, alpha,
                  a, lda, b, ldb, beta, c, ldc);
      }
      end = getSystemTime();
      printf("%i*%ix%i*%i: %f ms\n",length,length,length,length,(double)(end-start)/iter);

      //free
      mkl_free(a);
      mkl_free(b);
      mkl_free(c);

      return 0;
}
 
int algoPerform(int length,int iter){
  long long start,end;
  MKL_INT i=0,vec_len=length;
  float fA[vec_len],fB[vec_len],fBha[vec_len];
  srand((unsigned)time(NULL));

  for(;i<vec_len;i++){
    fA[i]=rand();
    fB[i]=rand();
  }
  start = getSystemTime();
  vsAdd(vec_len,fA,fA,fBha);
  for(i=0;i<iter;i++){
    vsAdd(vec_len,fA,fA,fBha);
  }
  end = getSystemTime(); 

  printf("vsadd Performtest: %i*%i %f millis\n",length,length,(double)(end-start)/iter); 

  start = getSystemTime();
  vsSub(vec_len,fA,fA,fBha);
  for(i=0;i<iter;i++){
    vsSub(vec_len,fA,fA,fBha);
  }
  end = getSystemTime();
  printf("vssub Performtest: %i*%i %f millis\n",length,length,(double)(end-start)/iter);

  start = getSystemTime();
  vsMul(vec_len,fA,fA,fBha);
  for(i=0;i<iter;i++){
    vsMul(vec_len,fA,fA,fBha);
  }
  end = getSystemTime();
  printf("vsMul Performtest: %i*%i %f millis\n",length,length,(double)(end-start)/iter);

  start = getSystemTime();
  vsDiv(vec_len,fA,fA,fBha);
  for(i=0;i<iter;i++){
    vsDiv(vec_len,fA,fA,fBha);
  }
  end = getSystemTime();
  printf("vsdiv Performtest: %i*%i %f millis\n",length,length,(double)(end-start)/iter);

  return 0;
}

int main(){
 // algoPerform(4096,1000000);
 // algoPerform(512,10000000);
 // algoPerform(32,100000000);

      sgemmPerform(4096,10);
      sgemmPerform(512,200);
      sgemmPerform(32,100000);

      sgemvPerform(4096,300);
      sgemvPerform(512,100000);
      sgemvPerform(32,1000000);

}

