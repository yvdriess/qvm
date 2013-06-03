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
};

size_t num_amplitudes( quantum_state * qstate ) {
  return 1 << qstate->size;
}

quantum_state * init_quantum_state( quantum_state * qstate, uint8_t size ) {
  qstate->vector = malloc( sizeof(amplitude) * num_amplitudes(qstate) );
  assert(qstate->vector);
  return qstate;
}

void free_quantum_state( quantum_state * qstate ) {
  free( qstate->vector );
  free( qstate );
  qstate = NULL;
}

void quantum_print_amplitude( quantum_state * qstate, size_t index ) {
  char buffer[65];
  itoa(index,buffer,2);
  assert( index < num_amplitudes(qstate) );
  amplitude z = qstate->vector[index];
  printf( " %f +%fi|%s> (%e)\n", 
	  creal(z), cimag(z), buffer, cabs(z) );
}

void quantum_print_qstate( quantum_state * qstate ) {
  for( size_t i=0; i < num_amplitudes(qstate); ++i )
    quantum_print_amplitude( qstate, i );
}

void test_quantum_state() {
  quantum_state * qstate;
  init_quantum_state( qstate, 8 );
  assert(qstate->vector != NULL);
  assert(qstate->size == 256);
  srand(293485);
  for( int i=0; i < num_amplitudes(qstate); ++i )
    qstate->vector[i] = cexp(rand());
  free_quantum_state( qstate );
}
