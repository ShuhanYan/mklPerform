#include "mkl.h"
#include <sys/timeb.h>
#include <stdio.h>
#include "mkl_vml.h"
#include "mkl_types.h"
#include "mkl_cblas.h"
#include <stdlib.h>
long long getSystemTime() {
    struct timeb t;
    ftime(&t);
    return 1000 * t.time + t.millitm;
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

      long long start,end;
      srand((unsigned)time(NULL));

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
      for (i=0;i<warmIter;i++){
          cblas_sgemm(layout, transA, transB, m, n, k, alpha,
                  a, lda, b, ldb, beta, c, ldc);
      }
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
 
int vectorMathPerform(int length,int warmIter,int iter){
  int i=0,vec_len=length;
  int scalar = 5; 
  float fA[vec_len],fB[vec_len],fBha[vec_len];

  long long start,end;
  srand((unsigned)time(NULL));

  //generate random data
  for(i=0;i<vec_len;i++){
    fA[i]=rand();
    fB[i]=rand();
  }

  //add
  for (i=0;i<warmIter;i++){
    vsAdd(vec_len,fA,fB,fBha);
  }

  start = getSystemTime();
  for(i=0;i<iter;i++){
    vsAdd(vec_len,fA,fB,fBha);
  }
  end = getSystemTime(); 

  printf("vsadd Performtest: %i*%i %f millis\n",length,length,(double)(end-start)/iter); 

  //sub
  for (i=0;i<warmIter;i++){
    vsSub(vec_len,fA,fB,fBha);
  }

  start = getSystemTime();
  for(i=0;i<iter;i++){
    vsSub(vec_len,fA,fB,fBha);
  }
  end = getSystemTime();
  printf("vssub Performtest: %i*%i %f millis\n",length,length,(double)(end-start)/iter);

  //mul
  for(i=0;i<warmIter;i++){
    vsMul(vec_len,fA,fB,fBha);
  }
  
  start = getSystemTime();
  for(i=0;i<iter;i++){
    vsMul(vec_len,fA,fB,fBha);
  }
  end = getSystemTime();
  printf("vsMul Performtest: %i*%i %f millis\n",length,length,(double)(end-start)/iter);

  //div
  for(i=0;i<warmIter;i++){
    vsDiv(vec_len,fA,fB,fBha);
  }

  start = getSystemTime();
  for(i=0;i<iter;i++){
    vsDiv(vec_len,fA,fB,fBha);
  }
  end = getSystemTime();
  printf("vsdiv Performtest: %i*%i %f millis\n",length,length,(double)(end-start)/iter);

  // pow
  for(i=0;i<warmIter;i++){
    vsPowx(vec_len,fA,scalar,fBha);
  }

  start = getSystemTime();
  for(i=0;i<iter;i++){
    vsPowx(vec_len,fA,scalar,fBha);
  }
  end = getSystemTime();
  printf("vspow Performtest: %i*%i %f millis\n",length,length,(double)(end-start)/iter);

  //log
  for(i=0;i<warmIter;i++){
    vsLn(vec_len,fA,fBha);
  }

  start = getSystemTime();
  for(i=0;i<iter;i++){
    vsLn(vec_len,fA,fBha);
  }
  end = getSystemTime();
  printf("vslog10 Performtest: %i*%i %f millis\n",length,length,(double)(end-start)/iter);

  //exp
  for(i=0;i<warmIter;i++){
    vsExp(vec_len,fA,fBha);
  }

  start = getSystemTime();
  for(i=0;i<iter;i++){
    vsExp(vec_len,fA,fBha);
  }
  end = getSystemTime();
  printf("vsexp Performtest: %i*%i %f millis\n",length,length,(double)(end-start)/iter);

  //sqrt
  for(i=0;i<warmIter;i++){
    vsSqrt(vec_len,fA,fBha);
  }

  start = getSystemTime();
  for(i=0;i<iter;i++){
    vsSqrt(vec_len,fA,fBha);
  }
  end = getSystemTime();
  printf("vssqrt Performtest: %i*%i %f millis\n",length,length,(double)(end-start)/iter);

  //log1p
  for(i=0;i<warmIter;i++){
    vsLog1p(vec_len,fA,fBha);
  }

  start = getSystemTime();
  for(i=0;i<iter;i++){
    vsLog1p(vec_len,fA,fBha);
  }
  end = getSystemTime();
  printf("vslog1p Performtest: %i*%i %f millis\n",length,length,(double)(end-start)/iter);

  return 0;
}

int main(){
      vectorMathPerform(4096*4096,10,100000);
      vectorMathPerform(512,10,1000000);
      vectorMathPerform(32,10,10000000);

     // sgemmPerform(4096,10,10);
     // sgemmPerform(512,10,200);
     // sgemmPerform(32,10,100000);

     // sgemvPerform(4096,10,300);
     // sgemvPerform(512,10,100000);
     // sgemvPerform(32,10,1000000);

}

