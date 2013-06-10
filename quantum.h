#ifndef QUANTUM_H
#define QUANTUM_H

/* quantum state abstraction
 * this was created to replace the libquantum sparse quantum state representation
 * with a dense vector version */

#include <stdlib.h>
#include <complex.h>
#include <assert.h>
#include <stdint.h>

typedef float complex amplitude;

/* warning, quantum state size is the number of qubits
 * NOT the number of amplitudes, which is always 2^size */
typedef struct quantum_state {
  uint8_t     size;
  amplitude * vector;
} quantum_state_t;

size_t num_amplitudes( const quantum_state_t * restrict qstate ) {
  return 1 << qstate->size;
}

quantum_state_t * init_quantum_state( quantum_state_t * qstate, uint8_t size ) {
  assert( qstate != NULL && size >= 0 );
  qstate->size = size;
  qstate->vector = calloc( num_amplitudes(qstate), sizeof(amplitude) );
  return qstate;
}

void free_quantum_state( quantum_state_t * qstate ) {
  qstate->size = 0;
  free( qstate->vector );
  qstate->vector = NULL;
}

void quantum_copy_qstate( const quantum_state_t * restrict src,
			        quantum_state_t * restrict dst ) {
  assert(dst->size == src->size);
  memcpy( dst->vector, src->vector, sizeof(amplitude)*num_amplitudes(src) );
}

float quantum_prob( amplitude amp ) {
  return crealf(amp) * crealf(amp) + cimagf(amp) * cimagf(amp);
}

void quantum_normalize( quantum_state_t * restrict qstate ) {
  double limit = 1.0e-8;
  double norm=0;
  for( pos_t i=0; i< num_amplitudes(qstate) ; ++i ) {
    norm += quantum_prob(qstate->vector[i]);
  }
  norm = sqrt(norm);
  if( abs(1-norm) < limit )
    for( pos_t i=0; i< num_amplitudes(qstate) ; ++i ) {
      qstate->vector[i] /= norm;
    }
}


quantum_state_t* quantum_kronecker( const quantum_state_t * restrict qstate1, 
 				    const quantum_state_t * restrict qstate2,
				          quantum_state_t * restrict out ) {
  const size_t siz1 = num_amplitudes(qstate1);
  const size_t siz2 = num_amplitudes(qstate2);
  assert( qstate1->size + qstate2->size == out->size );
  assert( num_amplitudes(out) == siz1*siz2);
  for( pos_t i1=0; i1<siz1; ++i1 ) {
    const pos_t block_idx = i1 * siz2;
    memcpy( out->vector+block_idx, qstate2->vector, sizeof(amplitude)*siz2 );
    for( pos_t i2=0; i2<siz2; ++i2 )
      out->vector[block_idx+i2] *= qstate1->vector[i1];
  }
  return out;
}

char* binrep(unsigned int val, char *buff, int sz);

void quantum_print_amplitude( const quantum_state_t * qstate, size_t index ) {
  char buffer[65];
  assert( index < num_amplitudes(qstate) );
  amplitude z = qstate->vector[index];
  printf( " %f %+fi|%s> (%e)\n", 
	  creal(z), cimag(z), binrep(index, buffer, qstate->size), cabs(z) );
}

void quantum_print_qstate( const quantum_state_t * qstate ) {
  for( size_t i=0 ; i < num_amplitudes(qstate) ; ++i )
    quantum_print_amplitude( qstate, i );
}

void test_quantum_state() {
  quantum_state_t * qstate;
  init_quantum_state( qstate, 8 );
  assert(qstate->vector != NULL);
  assert(qstate->size == 256);
  srand(293485);
  for( int i=0; i < num_amplitudes(qstate); ++i )
    qstate->vector[i] = cexp(rand());
  free_quantum_state( qstate );
}

/* Create a string of binary digits based on the input value.
   Input:
       val:  value to convert.
       buff: buffer to write to must be >= sz+1 chars.
       sz:   number of digits to convert
   Returns address of string or NULL if not enough space provided.
*/
char *binrep (unsigned int val, char *buff, int sz) {
  *(buff+sz--) = '\0';
  /* For each bit (going backwards) store character. */
  while (sz >= 0) {
    *(buff+sz--) = ((val & 1) == 1) ? '1' : '0';
    /* Get next bit. */
    val >>= 1;
  }
  return buff;
}

#endif
