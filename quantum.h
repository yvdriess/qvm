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
  free( qstate->vector );
  qstate->size = 0;
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

  
size_t permute(const size_t i, const size_t m, const size_t n) {
  return n * ( i % m ) + i / m;
}

/* void test_permute(size_t m, size_t n) { */
/*   printf("testing permute m=%u, n=%u: { ",m,n); */
/*   for( int i=0; i<m*n; ++i ) */
/*     printf("%u, ", permute(i,m,n)); */
/*   printf("}\n"); */
/* } */


/* void cycle(const quantum_state_t * restrict in,  */
/* 	   const unsigned int               cycle, */
/* 	         quantum_state_t * restrict out) { */
/*   assert( in ); */
/*   assert( out ); */
/*   assert( in->size == out->size ); */
/*   const uint8_t old_cycle = in->cycled; */
/*   const uint8_t      size = in->size; */
/*   // to_cycle: how much times do I need to left cycle such that we achieve 'cycle' perm. */
/*   //    P_k P_l = P_(k+l)%n;   */
/*   // taking x = to_cycle, k = old_cycle, l = cycle and n = size: */
/*   //    P_k P_x = P_l  ==> P_x = P_(n-k) P_l ==> P_x = P_(n-k+l) */
/*   const uint8_t  to_cycle = cycles_to(qstate, cycle); */
/*   const size_t     stride = 1<<to_cycle; */
/*   const size_t  blocksize = 1<<(size - to_cycle); */
/*   out->cycled = cycle; */
  
/*   // dumb permutation to start with */
/*   for( int i=0; i<num_amplitudes(in); ++i ) { */
/*     size_t p_i = permute( i, stride, blocksize ); */
/*     out->vector[i] = in->vector[p_i]; */
/*   } */
/* } */

unsigned int 
cycles_to_back( const quantum_state_t qstate, 
		const unsigned int    new_cycle ) {
  return abs(qstate.size - qstate.cycled + new_cycle) % qstate.size;
}

// permute amplitudes, realizing qubit position cyclic shift such that 
//  qubit position 'target' comes last
// In most cases, this function will allocate and deallocate an amplitude vector
quantum_state_t 
shift_to_back(       quantum_state_t qstate, 
	       const pos_t           target ) {
  assert( qstate.size > 0 );
  assert( target < qstate.size );
  const uint8_t cycles = cycles_to_back( qstate, target );
  if( cycles == 0 || cycles == qstate.size )
    return qstate;
  else {
    quantum_state_t out;
    init_quantum_state( &out, qstate.size );
    const size_t    stride = 1<<cycles;
    const size_t blocksize = 1<<(qstate.size - cycles);
    // dumb permutation to start with
    for( int i=0; i<num_amplitudes(&qstate); ++i ) {
      size_t p_i = permute( i, stride, blocksize );
      out.vector[i] = qstate.vector[p_i];
    }
    free_quantum_state( &qstate );
    return out;
  }
}

quantum_state_t shift_back( const quantum_state_t qstate ) {
  return shift_to_back( qstate, 0 ); 
}
							   
quantum_state_t* 
quantum_kronecker( const quantum_state_t * restrict qstate1, 
		   const quantum_state_t * restrict qstate2,
		         quantum_state_t * restrict out ) {
  const size_t siz1 = num_amplitudes(qstate1);
  const size_t siz2 = num_amplitudes(qstate2);
  assert( qstate1->size + qstate2->size == out->size );
  assert( num_amplitudes(out) == siz1*siz2);
  quantum_state_t left  = shift_back( *qstate1 );
  quantum_state_t right = shift_back( *qstate2 );
  #pragma omp parallel for
  for( pos_t i1=0; i1<siz1; ++i1 ) {
    const pos_t block_idx = i1 * siz2;
    const amplitude left_element = left.vector[i1];
    // copy all 'right' elements as a block
    memcpy( out->vector+block_idx, right.vector, sizeof(amplitude)*siz2 );
    for( pos_t i2=0; i2<siz2; ++i2 )
      out->vector[block_idx+i2] *= left_element;
  }
  return out;
}

char* binrep(unsigned int val, char *buff, int sz);

void quantum_print_amplitude( const quantum_state_t * qstate, size_t index ) {
  const double          limit = 1.0e-8;
  const size_t           size = qstate->size;
  const size_t         stride = 1 << qstate->cycled;
  const size_t      blocksize = 1 << (size - qstate->cycled);
  // inverted permutation of permute(index, stride, blocksize): P_m,n P_n,m = I
  const size_t permuted_index = permute(index, blocksize, stride);
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
  printf("printing quantum state cycled by %u\n", qstate->cycled);
  for( size_t i=0 ; i < num_amplitudes(qstate) ; ++i )
    quantum_print_amplitude( qstate, i );
  //test_permute(1<<(qstate->size-qstate->cycled), 1<<qstate->cycled);
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
