gcc -O3 -Lsexp/lib/ -Isexp/include/ -o qvm qvm.c -lsexp -lm -std=c99 -g3 -fopenmp
gcc -O3 -Lsexp/lib/ -Isexp/include/ -o qvm_ qvm.c -lsexp -lm -std=c99 -g3 -fopenmp -DNOPERM
