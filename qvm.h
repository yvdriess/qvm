#ifndef QVM_H
#define QVM_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif 

#ifndef LEVEL1_DCACHE_LINESIZE
#define LEVEL1_DCACHE_LINESIZE 64
#endif

#include <stdint.h>
#include <complex.h>

typedef complex float amplitude;

typedef uint8_t qid_t;
typedef uint8_t tangle_size_t;
typedef size_t  pos_t;

#endif
