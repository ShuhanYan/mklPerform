#include "mkl.h"
#include <sys/timeb.h>
#include <stdio.h>
#include "mkl_vml.h"
#include "mkl_types.h"
#include "mkl_cblas.h"
#include <stdlib.h>
#include <sys/time.h>

long long getSystemTime() {
    struct timeb t;
    ftime(&t);
    return 1000 * t.time + t.millitm;
}

void init(float * data, const int length, const float min, const float max, const float interval){
  int i = 0, now = min;
  for(i = 0; i < length - 1; i++){
    if(now > max) now = min;
    data[i] = now;
    now = now + interval;
  }
}

int sgemvPerform(int length,int warmIter,int iter){
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

      struct timeval start,end;
      float timestamp = 0;

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


      //warm up
      for (i=0;i<warmIter;i++){
        cblas_sgemv(layout, trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
      }

      //test
      gettimeofday(&start,NULL);
      for(i=0;i<iter;i++){
          cblas_sgemv(layout, trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
      }
      gettimeofday(&end,NULL);
      timestamp = (end.tv_sec-start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
      printf("%i*%ix%i*%i: %f ms\n",length,length,length,length,(double)(timestamp)/iter);

      mkl_free(a);
      mkl_free(x);
      mkl_free(y);

      return 0;
}

int sgemmPerform(int length,int warmIter,int iter){
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

      struct timeval start,end;

      //generate data
      for(i = 0; i < rmax; i++){
         for(j = 0; j < cmax; j++){
             a[i*cmax+j] = rand();
             b[i*cmax+j] = rand();
             c[i*cmax+j] = 0;
         }
      }

      //warm up
      for (i=0;i<warmIter;i++){
          cblas_sgemm(layout, transA, transB, m, n, k, alpha,
                  a, lda, b, ldb, beta, c, ldc);
      }
      //test
      gettimeofday(&start,NULL);
      for (i=0;i<iter;i++){
          cblas_sgemm(layout, transA, transB, m, n, k, alpha,
                  a, lda, b, ldb, beta, c, ldc);
      }
      gettimeofday(&end,NULL);
      printf("%i*%ix%i*%i: %f ms\n",length,length,length,length,(double)(end-start)/iter);

      //free
      mkl_free(a);
      mkl_free(b);
      mkl_free(c);

      return 0;
}


int vectorMathPerform(int length,int warmIter,int iter, int range, int interval){
  int i=0,vec_len=length*length;
  int scalar = 5; 
  float *fA = (float *)mkl_calloc(vec_len, sizeof( float ), 64);//[vec_len],fB[vec_len],fBha[vec_len];
  float *fB = (float *)mkl_calloc(vec_len, sizeof( float ), 64);
  float *fBha = (float *)mkl_calloc(vec_len, sizeof( float ), 64);
  struct timeval start, end;
  float timestamp = 0;

  init(fA, vec_len, -range, range, interval);
  init(fB, vec_len, -range, range, interval);
  init(fBha, vec_len, -range, range, interval);

  //add
  for (i=0;i<warmIter;i++){
    vsAdd(vec_len,fA,fB,fBha);
  }

  gettimeofday(&start,NULL);
  for(i=0;i<iter;i++){
    vsAdd(vec_len,fA,fB,fBha);
  }
  gettimeofday(&end,NULL);
  timestamp =  (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000;
  
  printf("vsadd Performtest: %i*%i %f millis\n", length, length, (double) timestamp / iter); 

  //sub
  for (i=0;i<warmIter;i++){
    vsSub(vec_len,fA,fB,fBha);
  }

  gettimeofday(&start,NULL);
  for(i=0;i<iter;i++){
    vsSub(vec_len,fA,fB,fBha);
  }
  gettimeofday(&end,NULL);
  timestamp =  (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000;
  printf("vssub Performtest: %i*%i %f millis\n", length, length, (double) timestamp / iter);

  //mul
  for(i=0;i<warmIter;i++){
    vsMul(vec_len,fA,fB,fBha);
  }
  
  gettimeofday(&start,NULL);
  for(i=0;i<iter;i++){
    vsMul(vec_len,fA,fB,fBha);
  }
  gettimeofday(&end,NULL);
  timestamp =  (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000;
  printf("vsMul Performtest: %i*%i %f millis\n", length, length, (double) timestamp / iter);

  //div
  for(i=0;i<warmIter;i++){
    vsDiv(vec_len,fA,fB,fBha);
  }

  gettimeofday(&start,NULL);
  for(i=0;i<iter;i++){
    vsDiv(vec_len,fA,fB,fBha);
  }
  gettimeofday(&end,NULL);
  timestamp =  (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000;
  printf("vsdiv Performtest: %i*%i %f millis\n", length, length, (double) timestamp / iter);

  // pow
  for(i=0;i<warmIter;i++){
    vsPowx(vec_len,fA,scalar,fBha);
  }

  gettimeofday(&start,NULL);
  for(i=0;i<iter;i++){
    vsPowx(vec_len,fA,scalar,fBha);
  }
  gettimeofday(&end,NULL);
  timestamp =  (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000;
  printf("vspow Performtest: %i*%i %f millis\n", length, length, (double) timestamp / iter);

  //log
  for(i=0;i<warmIter;i++){
    vsLn(vec_len,fA,fBha);
  }

  gettimeofday(&start,NULL);
  for(i=0;i<iter;i++){
    vsLn(vec_len,fA,fBha);
  }
  gettimeofday(&end,NULL);
  timestamp =  (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000;
  printf("vslog Performtest: %i*%i %f millis\n", length, length, (double) timestamp / iter);

  //exp
  for(i=0;i<warmIter;i++){
    vsExp(vec_len,fA,fBha);
  }

  gettimeofday(&start,NULL);
  for(i=0;i<iter;i++){
    vsExp(vec_len,fA,fBha);
  }
  gettimeofday(&end,NULL);
  timestamp =  (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000;
  printf("vsexp Performtest: %i*%i %f millis\n", length, length, (double) timestamp / iter);

  //sqrt
  for(i=0;i<warmIter;i++){
    vsSqrt(vec_len,fA,fBha);
  }

  gettimeofday(&start,NULL);
  for(i=0;i<iter;i++){
    vsSqrt(vec_len,fA,fBha);
  }
  gettimeofday(&end,NULL);
  timestamp =  (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000;
  printf("vssqrt Performtest: %i*%i %f millis\n", length, length, (double) timestamp / iter);

  //log1p
  for(i=0;i<warmIter;i++){
    vsLog1p(vec_len,fA,fBha);
  }

  gettimeofday(&start,NULL);
  for(i=0;i<iter;i++){
    vsLog1p(vec_len,fA,fBha);
  }
  gettimeofday(&end,NULL);
  timestamp =  (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000;
  printf("vslog1p Performtest: %i*%i %f millis\n", length, length, (double) timestamp / iter);

  return 0;
}

int main(){
      vectorMathPerform(4096,10,100,1000,0.5);
      vectorMathPerform(512,10,10000,1000,0.5);
      vectorMathPerform(32,10,100000,500,1);

      sgemmPerform(4096,10,10);
      sgemmPerform(512,10,200);
      sgemmPerform(32,10,100000);

      sgemvPerform(4096,10,300);
      sgemvPerform(512,10,100000);
      sgemvPerform(32,10,1000000);

}

