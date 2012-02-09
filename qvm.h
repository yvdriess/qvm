#ifndef QVM_H
#define QVM_H

#include "libquantum/complex.h"
#include "libquantum/error.h"
#include <quantum.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif 

extern void quantum_copy_qureg(quantum_reg *src, quantum_reg *dst);
extern void quantum_delete_qureg_hashpreserve(quantum_reg *reg);

typedef int qid_t;
typedef int tangle_size_t;
typedef int pos_t;

#endif
