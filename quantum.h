#ifndef QUANTUM_H
#define QUANTUM_H

/* quantum state abstraction
 * this was created to replace the libquantum sparse quantum state representation
 * with a dense vector version */

#include <stdlib.h>
#include <complex.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>

typedef float complex amplitude;

/* warning, quantum state size is the number of qubits
 * NOT the number of amplitudes, which is always 2^size */
typedef struct quantum_state {
  uint8_t      size;
  uint8_t     cycled;
  amplitude * vector;
} quantum_state_t;

static quantum_state_t _the_empty_qstate_ = {0,0,NULL};

size_t num_amplitudes( const quantum_state_t * restrict qstate ) {
  return 1 << qstate->size;
}

quantum_state_t * init_quantum_state( quantum_state_t * qstate, uint8_t size ) {
  assert( qstate != NULL && size >= 0 );
  qstate->size = size;
  qstate->cycled = 0;
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

  
int permute(const int i, const int m, const int n) {
  return n * ( i % m ) + i / m;
}

void cycle(const quantum_state_t * restrict in, 
	   const int                        cycle,
	         quantum_state_t * restrict out) {
  assert( in );
  assert( out );
  assert( in->size == out->size );
  uint8_t old_cycle = in->cycled;
  uint8_t      size = in->size;
  uint8_t new_cycle = abs(old_cycle + cycle) % size;
  out->cycled = new_cycle;

  // dumb permutation to start with
  for( int i=0; i<num_amplitudes(in); ++i ) {
    out->vector[i] = in->vector[permute(i,  1<<new_cycle, 1<<(size - new_cycle) )];
  }
}

quantum_state_t* quantum_kronecker( const quantum_state_t * restrict qstate1, 
 				    const quantum_state_t * restrict qstate2,
				          quantum_state_t * restrict out ) {
  const size_t siz1 = num_amplitudes(qstate1);
  const size_t siz2 = num_amplitudes(qstate2);
  assert( qstate1->size + qstate2->size == out->size );
  assert( num_amplitudes(out) == siz1*siz2);
  #pragma omp parallel for
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
  const double limit = 1.0e-8;
  const size   = num_amplitudes(qstate);
  const stride = 1 << qstate->cycled;
  const size_t permuted_index = permute(index, size - stride, stride);
  const amplitude z = qstate->vector[permuted_index];
  const double prob = quantum_prob(z);
  char buffer[65];
  assert( index < num_amplitudes(qstate) );
  if( prob > limit )
    printf( "% f %+fi|%zd> (%e) (|%s>)\n", 
	    creal(z), cimag(z), index, 
	    quantum_prob(z), binrep(index, buffer, qstate->size) );
}

void quantum_print_qstate( const quantum_state_t * qstate ) {
  for( size_t i=0 ; i < num_amplitudes(qstate) ; ++i )
    quantum_print_amplitude( qstate, i );
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
