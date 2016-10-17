mkdir _result
icc -I/opt/intel/mkl/include -I /opt/jdk/include/ -I /opt/jdk/include/linux/ source/mklPerform.c  -liomp5 -lpthread -lm -liomp5 -lc -lmkl_rt -fopenmp -Wall -fPIC -std=c99  -o _results/mklPerform.out
_results/mklPerform.out >_results/mklPerform.res

