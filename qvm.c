#include <fcntl.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <stdint.h>
#include <complex.h>

#include <sexp.h>
#include <sexp_ops.h>
#include <sexp_vis.h>

//#include "libquantum/config.h"

#include "bitmask.h"
#include "qvm.h"
#include "quantum.h"

#define STRING_SIZE (size_t)UCHAR_MAX	
#define MAX_TANGLES (size_t)SHRT_MAX
#define MAX_QUBITS  (size_t)SHRT_MAX

#define car hd_sexp
#define cdr next_sexp
#define max(x,y) x < y ? x : y

#define _inv_sqrt2_ 0.707106781186548f
#define _diag_amp_ _inv_sqrt2_

int _verbose_ = 0;
int _interactive_ = 0;
int _alt_measure_ = 1;

// |+> state prototype
static quantum_state_t _proto_diag_qubit_;
// CZ|++> state prototype
static quantum_state_t _proto_dual_diag_qubit_;

/************
 ** TANGLE **
 ************/

typedef struct qid_list {
  qid_t qid;
  struct qid_list* rest;
} qid_list_t;

typedef struct tangle { 
  qid_t             size;
  qid_list_t      * qids;
  quantum_state_t   qstate;
} __attribute__ ((aligned(LEVEL1_DCACHE_LINESIZE))) tangle_t;

tangle_t* init_tangle() {
  tangle_t *restrict tangle = malloc( sizeof(tangle_t) );
  *tangle = (tangle_t){0,NULL,_the_empty_qstate_};
  return tangle;
}

void free_qid_list( qid_list_t* qids) {
  qid_list_t* rest = NULL;
  while( qids ) {
    rest = qids->rest;
    free( qids ); 
    qids = rest;
  }
}

void free_tangle( tangle_t* tangle ) {
  free_qid_list( tangle->qids );
  tangle->qids = NULL;
  tangle->size = 0;
  if( tangle->qstate.size > 0 ) {
    //free_quantum_state( &tangle->qstate );
    FREE_QUANTUM_STATE( tangle->qstate );
  }
  free( tangle );                                                        //FREE tangle
}

void print_qids( const qid_list_t* qids ) {
  printf("[");
  const qid_list_t* cons = qids;
  for( ;
       cons;
       cons=cons->rest ) {
    printf("%d", cons->qid);
    if( cons->rest )
      printf(", ");
  }
  printf("]");
}

void print_tangle( const tangle_t* tangle ) {
  assert( tangle );
  print_qids( tangle->qids );
  printf(" ,\n    {\n");
  if( tangle->size > 5 ) {
    printf("<a large quantum state>, ");
    if( _interactive_ ) {
      printf("really print? (y/N): ");
      if( getchar() == 'y' )
	quantum_print_qstate( tangle->qstate );
    }
    else
      printf("omitted due to size > 128");
  }
  else
    quantum_print_qstate( tangle->qstate );
  printf("}\n");
}

// destructively permutes the qid list, cyclically shifting 'cycle' elements to the left
qid_list_t* cycle_qid_list( qid_list_t * qids, unsigned int cycle ) {
  qid_list_t* old_head = qids;
  qid_list_t* current = qids;
  if( cycle == 0 ) // identity
    return qids; 
  while( --cycle > 0 ) {    // walk to the cutoff point
    assert( current );      // protect against inf. looping on bad input
    current = current->rest;
  }
  qid_list_t* new_tail = current;
  qid_list_t* new_head = current->rest;
  if( new_head == NULL )
    return qids; // cycle is as long as list => identity
  current = new_head;
  new_tail->rest = NULL; // snip
  while( current->rest ) { //seek end
    current = current->rest;
  }
  current->rest = old_head; // paste
  return new_head;
}

tangle_t * position_shift( tangle_t * restrict tangle, const qid_t qid ) {
  assert(tangle);
  assert(tangle->qstate.size > 0);
  quantum_state_t new_qstate = shift_to_back( tangle->qstate, qid );
  qid_list_t    * new_qids   = cycle_qid_list( tangle->qids, new_qstate.cycled );
  new_qstate.cycled = 0;
  tangle->qids   = new_qids;
  tangle->qstate = new_qstate;
  return tangle;
}

/***********
 ** QUBIT **
 ***********/
typedef struct qubit {
  tangle_t* tangle;
  qid_t qid;
  pos_t pos;
} qubit_t;

qubit_t _invalid_qubit_ = { NULL, -1, -1 };

// use this function to return a correct qstate position
//  libquantum uses an reverse order (least significant == 0)
pos_t get_target( const qubit_t qubit ) {
  return qubit.tangle->size - qubit.pos - 1;
}

bool invalid( const qubit_t qubit ) {
  return qubit.tangle == NULL;
}

quantum_state_t get_qstate( const qubit_t qubit ) {
  return qubit.tangle->qstate;
}

void set_qstate( qubit_t         qubit, 
		 quantum_state_t new_state ) {
  qubit.tangle->qstate = new_state;
}

/**********
 ** QMEM **
 **********/
typedef struct signal_map {
  // two bitfields
  //  entries : if qid has an entry, not needed for correct programs
  //  signals : value of the signal
  unsigned char entries[BITNSLOTS(MAX_QUBITS)];
  unsigned char signals[BITNSLOTS(MAX_QUBITS)];
} signal_map_t;

typedef struct qmem {
  size_t         size;
  signal_map_t   signal_map;
  tangle_t     * tangles[MAX_TANGLES];
} qmem_t;


void print_signal_map( const signal_map_t* signal_map ) {
  printf(" {\n");
  for( int qid=0 ; qid<MAX_QUBITS ; ++qid ) {
    if( BITTEST(signal_map->entries, qid) )
      printf("  %d -> %d,\n", qid,
	     BITTEST(signal_map->signals, qid) ? 1 : 0 );
  }  
  printf(" }\n");
}

bool get_signal( const qid_t qid, 
		 const signal_map_t* signal_map ) {
  if( BITTEST(signal_map->entries,qid) )
    return BITTEST(signal_map->signals, qid);
  else {
    printf( "ERROR: I was asked a signal map entry (qid:%d) that wasn't there,\n\
  check quantum program correctness.\n", qid);
    printf( "   signal map:\n     ");
    print_signal_map( signal_map );
    exit(EXIT_FAILURE);
  }

}
void set_signal( const qid_t qid, 
		 const bool signal, 
		 signal_map_t* signal_map ) {
  if( BITTEST(signal_map->entries, qid) ) {
    printf( "ERROR: I was asked to set an already existing signal (qid: %d),\n\
  check quantum program correctness.\n", qid);
    printf( "   signal map:\n     ");
    print_signal_map( signal_map );
    exit(EXIT_FAILURE);
  }
  BITSET(signal_map->entries, qid);
  if( signal )
    BITSET(signal_map->signals, qid);
}

qubit_t 
find_qubit_in_tangle( const qid_t qid, 
		      const tangle_t* restrict tangle )
{
  assert(tangle);
  qid_list_t * restrict qids = tangle->qids;
  if( tangle->size == 0 ) {
    printf("WARNING: looking for qid in empty tangle, this is not "
	   "supposed to happen (deallocate this tangle)\n");
    return _invalid_qubit_;
  }
  for( int i=0;
       qids;
       ++i, qids=qids->rest ) {
    if( qids->qid == qid ) 
      return (qubit_t){ (tangle_t*)tangle, qid, i };
  }
  return _invalid_qubit_;
}

qubit_t 
find_qubit(const qid_t qid, const qmem_t * restrict qmem) {
  tangle_t * restrict tangle;
  qubit_t qubit;
  for( int i=0, tally=0 ; tally < qmem->size ; ++i ) {
    tangle = qmem->tangles[i];
    if( tangle ) {
      qubit = find_qubit_in_tangle(qid, tangle);
      if( !invalid(qubit) )
	return qubit;      
      ++tally;
    }
  }
  return _invalid_qubit_;
}

qid_list_t* add_qid( const qid_t qid, qid_list_t * restrict qids ) {
  // assuming qid is NOT already in qids
  // ALLOC QUBIT LIST
  qid_list_t * restrict new_qids = (qid_list_t*) malloc(sizeof(qid_list_t));
  new_qids->qid = qid;
  new_qids->rest = qids;
  return new_qids;
}

void append_qids( qid_list_t * restrict new_qids, qid_list_t * restrict target_qids ) {
  assert( new_qids && target_qids );
  while( target_qids->rest ) {
    target_qids = target_qids->rest;
  }
  target_qids->rest = new_qids;  
}

/* void remove_qid( qid_t qid, qid_list_t* restrict qids ) { */
/*   qid_list_t* restrict next; */
/*   if( qids ) { */
/*     next = qids->rest; */
/*     if( next ) { */
/*       if( next->qid == qid ) { */
/* 	qids->rest = next->rest; */
/* 	free( next ); */
/*       } */
/*       else { */
/* 	remove_qid( qid, next ); */
/*       } */
/*     } */
/*   } */
/*   else */
/*     // recursion stops if qid was not found */
/*     printf("Warning: I was asked to remove qid %d from a tangle that did" */
/* 	   "not have it\n", qid); */
/* } */

void print_qmem( const qmem_t* qmem ) {
  assert(qmem);
  printf("qmem has %d tangles:\n  {", (int)qmem->size);
  for( int i=0, tally=0 ; tally < qmem->size ; ++i ) {
    assert(i<MAX_TANGLES);
    if( qmem->tangles[i] ) {
      if( tally>0 )
	printf(",\n   ");
      print_tangle(qmem->tangles[i]);
      ++tally;
    }
  }
  printf("}\n");
  //printf("signal map:");
  //print_signal_map( &(qmem->signal_map) );
}

qmem_t* init_qmem( qmem_t* restrict qmem ) {
  qmem->size = 0;
  //qmem->tangles = calloc(MAX_TANGLES,sizeof(tangle_t*)); //ALLOC tangles
  memset(&(qmem->tangles), 0, sizeof(tangle_t*) * MAX_TANGLES);
  //for( int i=0; i<MAX_TANGLES; ++i )
  //  qmem->tangles[i] = NULL;
  //  memset( &(qmem->signal_map),0,sizeof(signal_map_t) );
  qmem->signal_map = (signal_map_t){ {0}, {0} };
  
  _proto_diag_qubit_ = NEW_QUANTUM_STATE(1);
  _proto_dual_diag_qubit_ = NEW_QUANTUM_STATE(2);
  _proto_diag_qubit_.vector[0] = _diag_amp_;
  _proto_diag_qubit_.vector[1] = _diag_amp_;
  _proto_dual_diag_qubit_.vector[0] = _diag_amp_/2;
  _proto_dual_diag_qubit_.vector[1] = _diag_amp_/2;
  _proto_dual_diag_qubit_.vector[2] = _diag_amp_/2;
  _proto_dual_diag_qubit_.vector[3] = -_diag_amp_/2;

  // in the dense amplitide vector implementation, these prototypes are initialized as static
  // instantiate prototypes (libquantum qstates)
  //_proto_diag_qubit_ = quantum_new_qstate(0, 1);
  //_proto_dual_diag_qubit_ = quantum_new_qstate(0, 2);
  //  quantum_hadamard(0, &_proto_diag_qubit_);
  //quantum_hadamard(0, &_proto_dual_diag_qubit_);
  //quantum_hadamard(1, &_proto_dual_diag_qubit_);
  //quantum_gate2(0, 1, _cz_gate_, &_proto_dual_diag_qubit_);  

  // seed RNG
  //sranddev();
  srand(time(0));
  return qmem;
}

void free_qmem(qmem_t* qmem) {
  for( int i=0, tally=0 ; tally < qmem->size ; i++ ) {
    assert(i<MAX_TANGLES);
    if( qmem->tangles[i] ) {
      free_tangle(qmem->tangles[i]);
      qmem->tangles[i] = NULL;
      ++tally;
    }
  }
}

tangle_t* get_free_tangle(qmem_t* qmem) {
  tangle_t* restrict new_tangle = init_tangle();
  assert(new_tangle);
  // I simply loop here because tangles can get de-allocated (NULL-ed)
  for(int i=0; i<MAX_TANGLES; ++i) {
    if( qmem->tangles[i] == NULL ) {
      qmem->tangles[i] = new_tangle;
      return new_tangle;
    }
  }
  printf("Ran out of qmem memory, too many tangles! (> %lu)\n",MAX_TANGLES);
  exit(EXIT_FAILURE);
}

tangle_t*
add_dual_tangle( const qid_t qid1, 
		 const qid_t qid2, 
		 qmem_t* restrict qmem) {
  // allocate new tangle in qmem
  tangle_t*  restrict tangle = get_free_tangle(qmem);

  // init tangle
  tangle->qids = add_qid( qid2, tangle->qids );
  tangle->qids = add_qid( qid1, tangle->qids );
  tangle->size = 2;

  // update qmem info
  qmem->size += 1;

  // init quantum state
  tangle->qstate = NEW_QUANTUM_STATE(2);
  quantum_copy_qstate( &_proto_dual_diag_qubit_,
		       &tangle->qstate );
  return tangle;
}

tangle_t* 
add_tangle( const qid_t qid, 
	    qmem_t* restrict qmem ) {
  // allocate new tangle in qmem
  tangle_t*  restrict tangle = get_free_tangle(qmem);
  // init tangle
  tangle->qids = add_qid( qid, tangle->qids );
  tangle->size = 1;
  // update qmem info
  qmem->size += 1;
  // init quantum state
  tangle->qstate = NEW_QUANTUM_STATE(1);
  quantum_copy_qstate(&_proto_diag_qubit_,
 		      &tangle->qstate);
  return tangle;
}


/* Adds new qubit BEHIND existing state: |q> x |+>
 */
void
add_qubit( const qid_t qid, 
	   tangle_t* restrict tangle) {
  quantum_state_t  old_qstate = tangle->qstate;
  // tensor |+> to quantum state
  quantum_state_t  new_qstate = quantum_kronecker( old_qstate, _proto_diag_qubit_ );
  // appends new qid:  qids := [[qids...],qid]
  append_qids( add_qid(qid,NULL), tangle->qids );
  // out with the old
  FREE_QUANTUM_STATE( old_qstate );
  // in with the new
  tangle->qstate = new_qstate;
  tangle->size += 1;
}

void 
delete_tangle( tangle_t * restrict tangle,
	       qmem_t   * restrict qmem ) {
  assert( tangle );
  assert( tangle->qids == NULL );
  qmem->size -= 1;
  // null the tangle entry in qmem
  for(int i=0; i<MAX_TANGLES; ++i) {
    if( qmem->tangles[i] == tangle ) {
      qmem->tangles[i] = NULL;
      free_tangle( tangle );
      return;
    }
  }
  free_tangle( tangle );
  printf("ERROR: I was asked to delete an unknown tangle in qmem\n");
  exit(EXIT_FAILURE);
}

void 
delete_qubit( const qubit_t            qubit, 
	            qmem_t  * restrict qmem) {
  assert( !invalid(qubit) );
  tangle_t * tangle = qubit.tangle;
  // handle points to where current qid entry is stored
  qid_list_t** handle = &tangle->qids;
  qid_list_t*    qids = tangle->qids;
  
  // listref by qubit.pos, handle is the previous 'next' pointer 
  //   that will need to be changed
  for( int i=0 ; i < qubit.pos ; ++i ) {
    handle = &(qids->rest);
    qids = qids->rest;
  }
  assert(qids);
  assert(qids->qid == qubit.qid);

  // pointer plumbing:
  //  bypass current qid entry
  *handle = qids->rest;
  tangle->size -= 1;

  free(qids); // FREE QUBIT LIST element
 // when empty, dealloc tangle
  if( tangle->size == 0 ) {
    tangle->qids = NULL;
    delete_tangle( tangle, qmem );
  }
}

// WARNING: re-uses tangle_1, the resulting merged tangle will be put into tangle_1
void 
merge_tangles( tangle_t * restrict tangle_1, 
	       tangle_t * restrict tangle_2,
	       qmem_t   * restrict qmem) {
  assert( tangle_1 && tangle_2 );
  // append qids of tangle_2 to tangle_1, destructively
  append_qids( tangle_2->qids, tangle_1->qids);
  // tensor both qstates
  quantum_state_t qstate_1 = tangle_1->qstate;
  quantum_state_t qstate_2 = tangle_2->qstate;
  quantum_state_t new_qstate = quantum_kronecker( qstate_1, qstate_2 );
  // out with the old, but we are keeping tangle_1
  FREE_QUANTUM_STATE( qstate_1 ); // will be replaced, manually free the old
  tangle_2->qids = NULL; // avoids its qid_list nodes (now managed by tangle_1) from being freed
  delete_tangle( tangle_2, qmem ); //free the tangle (and qstate)
  // in with the new
  tangle_1->qstate = new_qstate;
  tangle_1->size   = new_qstate.size;
}

void ensure_list( sexp_t* exp ) {
  CSTRING* str = NULL;
  if( exp->ty != SEXP_LIST ) {
    print_sexp_cstr( &str, exp, STRING_SIZE );
    printf("ERROR: malformed expression, expecting a list expression and\
  got:\n  %s", toCharPtr( str ));
    exit(EXIT_FAILURE);
  }


}

void ensure_value( sexp_t* exp ) {
  CSTRING* str = NULL;

  if ( exp->ty != SEXP_VALUE ) {
    print_sexp_cstr( &str, exp, STRING_SIZE );
    printf("ERROR: malformed expression, expecting a value expression and\
  got:\n  %s", toCharPtr( str ));
    sdestroy( str );
    exit(EXIT_FAILURE);
    sdestroy( str );
  }
}

char get_opname( sexp_t* exp ) {
  return exp->val[0];
}

int get_qid( sexp_t* exp ) {
  return atoi( exp->val );
}

typedef struct angle_constant {
  const char* name;
  double value;
} angle_constant_t;

#define ANGLE_CONSTANT_MAX_CHARS 8
#define ANGLE_CONSTANTS_MAX 32
angle_constant_t _angle_constants_[ANGLE_CONSTANTS_MAX] = {
  {"PI",M_PI},
  {"PI/2",M_PI/2},
  {"PI/4",M_PI/4},
  {"PI/8",M_PI/8},
  {"-PI",-M_PI},
  {"-PI/2",-M_PI/2},
  {"-PI/4",-M_PI/4},
  {"-PI/8",-M_PI/8}
};
int _angle_constants_free_ = 8;

void add_new_constant(const char* name, double value) {
  
  if( _angle_constants_free_ > sizeof(_angle_constants_) ) {
    printf("ERROR: I can only remember %lu angle constants, I was asked"
	   " to add one more\n", sizeof(_angle_constants_));
    exit(EXIT_FAILURE);
  }
  angle_constant_t* restrict entry = 
    &_angle_constants_[_angle_constants_free_++];
  angle_constant_t new_entry = { name, value }; 
  *entry = new_entry;
}

double lookup_angle_constant(const char* str) {
  if( str )
    for( int i=0; i<_angle_constants_free_; ++i ) {
      if( strcmp(str, _angle_constants_[i].name) == 0 )
	return _angle_constants_[i].value;
    }
  return 0.0;
}


double parse_angle( const sexp_t* exp ) {
  /* syntax:
       <angle>  ::=  (- <angle>) | *angle_constant* | *float*  
  */ 
  //I can be more advanced and add some calc functionality,
  // but I'm not that insane atm. (some lib?)

  if( exp->ty == SEXP_LIST ) {
    // (- <angle>)
    const sexp_t* sign = exp->list;
    if( strcmp(sign->val, "-") == 0 )
      return -parse_angle(sign->next);
    else {
      printf("ERROR: expected (- ...) while parsing angle,"
	     "gotten:%s\n", sign->val);
      exit(EXIT_FAILURE);
    }
  }

  double angle = strtod( exp->val, NULL );
  if( angle == 0.0 && exp->val[0]!='0' ) { // atof failed
    // upper case 'str'
    char str[ANGLE_CONSTANT_MAX_CHARS];
    strcpy( str, exp->val );
    // ensure there is a termination string
    str[ANGLE_CONSTANT_MAX_CHARS-1] = 0;
    for(int i=0; str[i]; ++i) 
      str[i] = toupper(str[i]);

    // maybe it is specified in the environment?
    char* env_result = getenv(str);
    if( env_result ) {
      angle = atof( env_result );
      if( angle != 0.0 )
	return angle;
      else { // perhaps env contains one of our internal constants?
	for(int i=0; env_result[i]; ++i) 
	  env_result[i] = toupper(env_result[i]);
	return lookup_angle_constant(env_result);
      }
    }
    else {
      // maybe it's one of our constants
      angle = lookup_angle_constant(str);
      if( angle != 0.0 ) 
	return angle;
      else {
	// fallthrough: I really don't know what to do with this constant,
	// ask for input to the user
	printf(" angle \"%s\" is not a recognised constant, "
	       "please insert value: \n", str);
	int scanresult = scanf("%lf",&angle);
	if( scanresult ) {
	  printf(" added %s as %lf\n",str,angle);
	  add_new_constant(str, angle); 
        }
        else
          exit(EXIT_FAILURE);
      }
    }
  }
  return angle;
}

/************************
 ** QUANTUM OPERATIONS **
 ************************/

quantum_state_t
quantum_diag_measure( pos_t pos, double angle, const quantum_state_t qstate )
{
  //  double limit = (1.0 / ((MAX_UNSIGNED) 1 << qstate->width)) / 1000000;
  //  double prob=0, norm = 0;
  quantum_state_t out = NEW_QUANTUM_STATE( qstate.size - 1 );
  assert( pos < qstate.size );

  // take pos2 == 0b0010, 
  //  looping through the amplitude indices flips the bit at pos2
  //  after every pos2 iterations, this is the 'block size'
  //  this 'block' is repeated N/pos2 times
  // the way the diagonal measurement works, 
  //  every amplitude whose index at bit pos2 is 0, 
  //  is simply copied to the output vector
  //  amplitudes with index at pos2 of 1,
  //  are first multiplied by some factor
  //  both results are then summed

  // Method one: 
  //  1) loop through 'even' blocks, copy the blocks to the result
  //  2) loop through 'odd' blocks, multiply and add to result
 
  const amplitude factor = cexp(-angle*I);

#ifdef NOPERM
  // Method 1: nested loop per two blocks
  //   011101010110100110101
  //   | odd| ev.| odd| ev.|
  //   |         |    ^----^ : block_size
  //   |         |         |
  //             ^---------^ : block_stride
  const pos_t block_size   = 1 << pos;
  const pos_t block_stride = 1 << (pos + 1);
  assert( block_stride <= num_amplitudes(qstate) );
#pragma omp parallel for collapse(2)
  //         1000000000000 : block_size
  //        10000000000000 : block_stride
  // 0b1110101001111001001
  //         ^- target qubit (position 'pos'), pos=0 is lsb
  for( pos_t even_block=0;
       even_block < num_amplitudes(qstate) ; 
       even_block += block_stride) {
    for( pos_t i=0 ; i<block_size ; ++i ) {
      const pos_t odd_block = even_block + block_size;
      const pos_t out_block = even_block >> 1;
      const amplitude even = qstate.vector[even_block+i];
      const amplitude odd  = qstate.vector[odd_block+i];
      out.vector[out_block+i] = even - odd * factor;
    }
  }
  FREE_QUANTUM_STATE( qstate );
#else
  // Method 2: permutation step to collect contiguously accessed blocks
  // two permutations: a gather operation followed by block permutation (block permutation is achieved by prefetching/blocked access)
  quantum_state_t in = shift_to_back( qstate, pos );
  // Perform M on the 'last' position
#pragma omp parallel for
  for( pos_t i=0; i<num_amplitudes(out); i+=1 ) {
    const amplitude even = in.vector[i*2];
    const amplitude odd  = in.vector[i*2+1];
    out.vector[i] = even - odd * factor;
  }
  // measured out a qubit, adjust cycle (one less to cycle)
  if( out.cycled )
    out.cycled--;
#endif
  return out;
}

void qop_cz( const qubit_t qubit_1, const qubit_t qubit_2 ) {
  quantum_state_t qstate = get_qstate( qubit_1 );
  const pos_t tar1 = get_target(qubit_1);
  const pos_t tar2 = get_target(qubit_2);
  assert( !(invalid(qubit_1) || invalid(qubit_2)) );
  assert( qubit_1.tangle == qubit_2.tangle );
  if( qubit_1.pos + qubit_2.pos == 1 ) // parallel position
#pragma omp parallel for
    for( int i=0 ; i<num_amplitudes(qstate) ; i+=4 )
      qstate.vector[i] *= -1;
  else {
    const size_t bitmask = (1 << tar1) | (1 << tar2);
#pragma omp parallel for
    for( int i=0 ; i<num_amplitudes(qstate) ; i++ )
      /* Flip the target bit of a basis state if the control bit is set */     
      if( (i & bitmask) == bitmask )
	qstate.vector[i] *= -1;
  }
}

int qop_m( pos_t pos, double angle, const qubit_t qubit ) {
  assert( !invalid(qubit) );
  quantum_state_t qstate = get_qstate(qubit);
  quantum_state_t out = quantum_diag_measure( pos, angle, qstate );
  set_qstate( qubit, out );
  return 1;
}

void qop_x( const qubit_t qubit ) {
  assert( !invalid(qubit) );
  const pos_t target = get_target(qubit);
  quantum_state_t qstate = get_qstate(qubit);

#ifdef NOPERM
  const pos_t block_size   = 1 << target;
  const pos_t block_stride = 1 << (target + 1);
  // Method 1: nested loop per two blocks
  #pragma omp parallel for
  for( pos_t even_block=0;
       even_block < num_amplitudes(qstate) ;
       even_block += block_stride ) {
    amplitude tmp[block_size];
    const pos_t odd_block = even_block + block_size;
    // swap, using tmp as scratchpad
    //  TODO : might be interesting to see if three for loops might give gcc the chance to better optimize
    memcpy( tmp,                      qstate.vector+even_block, sizeof(amplitude)*block_size );
    memcpy( qstate.vector+even_block, qstate.vector+odd_block,  sizeof(amplitude)*block_size );
    memcpy( qstate.vector+odd_block,  tmp,                      sizeof(amplitude)*block_size );
  }
#else
  quantum_state_t out = shift_to_back( qstate, target );
  // Perform X on the 'last' position, basically a swap
#pragma omp parallel for
  for( pos_t i=0; i<num_amplitudes(qstate); i+=2 ) {
    amplitude tmp = out.vector[i];
    out.vector[i] = out.vector[i+1];
    out.vector[i+1] = tmp;
  }
  set_qstate( qubit, out );
#endif
}  

void qop_z( const qubit_t qubit ) {
  assert( !invalid(qubit) );
  const pos_t target = get_target(qubit);
  quantum_state_t qstate = get_qstate(qubit);
#ifdef NOPERM
  const pos_t block_size   = 1 << target;
  const pos_t block_stride = 1 << (target + 1);
  // Method 1: nested loop per two blocks
#pragma omp parallel for collapse(2)
  for( pos_t even_block=0;
       even_block < num_amplitudes(qstate) ; 
       even_block += block_stride ) {
    for( pos_t i=0 ; i<block_size; ++i ) {
      const pos_t odd_block = even_block + block_size;
      qstate.vector[odd_block+i] *= -1;
    }
  }
#else
  quantum_state_t out = shift_to_back( qstate, target );
  // perform Z on the 'last' position
  #pragma omp parallel for
  for( pos_t i=1; i < num_amplitudes(qstate); i+=2 ) {
    out.vector[i] *= -1;
  }
  set_qstate( qubit, out );
#endif
}

/***************
 ** EVALUATOR **
 ***************/
void eval_E(sexp_t* exp, qmem_t* qmem) {
  int qid1, qid2;
  qubit_t qubit_1;
  qubit_t qubit_2;

  assert( qmem );

  // move to the first argument
  exp = cdr(exp);
  if( !exp ) {
    printf("Entangle did not have any qubit arguments");
    exit(EXIT_FAILURE);
  }
  qid1 = get_qid( exp );
  
  // move to the second argument
  exp = cdr(exp);
  if( !exp ) {
    printf("Entangle did not have a second argument");
    exit(EXIT_FAILURE);
  }
  qid2 = get_qid( exp );

  // get tangle for qid1
  qubit_1 = find_qubit( qid1, qmem );
  qubit_2 = find_qubit( qid2, qmem );

  if( invalid(qubit_1) )
    if( invalid(qubit_2) ) {
      // if both unknown, create new tangle with two |+> states
      add_dual_tangle(qid1, qid2, qmem);
      return; // already in correct state by construction
    }
    else {
#ifndef NOPERM
      // make sure that qubit_2 is the in last position, 
      //  perform the necessary permutations if necessary
      qubit_2.tangle = position_shift( qubit_2.tangle, qubit_2.pos );
#endif
      // add qid1 to the 'end' of qid2's tangle
      add_qubit( qid1, qubit_2.tangle );
    }
  else
    if( invalid(qubit_2) ) {
#ifndef NOPERM
      // make sure that qubit_1 is the in last position, 
      //  perform the necessary permutations if necessary
      qubit_1.tangle = position_shift( qubit_1.tangle, qubit_1.pos );
#endif
      // add qid2 to qid1's tangle
      add_qubit( qid2, qubit_1.tangle );
    }
    else
      if( qubit_1.tangle == qubit_2.tangle ) {
	// if not, qubit entries are already valid
	qop_cz( qubit_1, qubit_2 );
	return;
      }
      else {
/* #ifndef NOPERM */
/*       // make sure that qubit_1 is the in last position,  */
/*       // make sure that qubit_2 is in the FIRST position */
/* 	const post_t new_q2_pos = (qubit_2.pos - 1 + qubit_2.tangle.size) % qubit_2.tangle.size; */
/* 	qubit_1.tangle = position_shift( qubit_1.tangle, qubit_1.pos ); */
/* 	qubit_2.tangle = position_shift( qubit_2.tangle, new_q2_pos  ); */
/* #endif */
	// both tangles are non-NULL, merge both
	merge_tangles(qubit_1.tangle, qubit_2.tangle, qmem);
      }
  // get valid qubit entries
  qubit_1 = find_qubit( qid1, qmem );
  qubit_2 = find_qubit( qid2, qmem );
  assert( qubit_1.tangle == qubit_2.tangle );
  qop_cz( qubit_1, qubit_2 );
}

/* Parses and checks the value of the given signal(s) */
/*   Syntax:  <identifier> | 0 | 1 | (q <qubit>) | (+ {<signal>}+ ) */
bool satisfy_signals( const sexp_t* restrict exp, 
		      const qmem_t* restrict qmem) {
  const sexp_t* args;
  const sexp_t* first_arg;
  bool signal;
  CSTRING* str = NULL;

  if( exp->ty == SEXP_LIST ) {
    args = exp->list;
    if( args->ty == SEXP_VALUE ) {
      if( strcmp(args->val, "q")==0 || 
	  strcmp(args->val, "Q")==0 ||
	  strcmp(args->val, "s")==0 ||
	  strcmp(args->val, "S")==0) {
	first_arg = args->next;
	return get_signal( atoi(first_arg->val), &qmem->signal_map );
      }
      else
	if( strcmp(args->val, "+")==0 ) {
	  first_arg = args->next;
	  signal = satisfy_signals(first_arg, qmem);
	  for(sexp_t* arg=first_arg->next; arg; arg=arg->next) {
	    signal ^= satisfy_signals(arg, qmem);
	  }
	  return signal;
	} 
    }// otherwise, fall through to parse_error
  }
  else
    if( exp->ty == SEXP_VALUE ) {
      if( strcmp(exp->val, "0") == 0 )
	return false;
      if( strcmp(exp->val, "1") == 0 )
	return true;
    } // otherwise, fall through to parse_error

  print_sexp_cstr( &str, exp, STRING_SIZE );
  printf("ERROR: I got confused parsing signal: %s\n", toCharPtr( str ));
  printf("  signal syntax:  <identifier> | 0 | 1 | (q <qubit>) |"
	 " (+ {<signal>}+ )\n");
  sdestroy(str);
  exit(EXIT_FAILURE);
}

void eval_M(sexp_t* exp, qmem_t* qmem) {
  int qid;
  double angle = 0.0;
  int signal;
  assert( qmem );

  // move to the first argument
  exp = cdr(exp);
  if( !exp ) {
    printf("Measurement did not have any target qubit argument\n");
    exit(EXIT_FAILURE);
  }
  qid = get_qid( exp );
  
  // move to the second argument
  exp = cdr(exp);
  if( exp ) { // default is 0
    angle = parse_angle( exp );
    // change angles by s- and t-signals when available
    exp = cdr(exp);
    if( exp ) { //s-signal, flips sign
      if( _verbose_ )
	printf("before angle correction, angle: %f\n", angle);
      if( satisfy_signals(exp, qmem) )
	angle = -angle;
      exp = cdr(exp);
      if( exp )  //t-signal, adds PI to angle
	if( satisfy_signals(exp, qmem) )
	  angle += M_PI;
    }
  }
  //  printf("  Measuring qubits %d\n",qid);
  qubit_t qubit = find_qubit( qid, qmem );
  if( invalid(qubit) ) {
    // create new qubit
    tangle_t* tangle = add_tangle( qid, qmem );
    qubit = find_qubit_in_tangle( qid, tangle );
  }
  if( _verbose_ )
    printf("  measuring qubit %d on angle %2.4f\n", qid, angle);
  signal = qop_m( get_target(qubit), angle, qubit );
  delete_qubit( qubit, qmem );
  /* printf("   result is %d\n",signal); */
  set_signal( qid, signal, &qmem->signal_map );
}

void eval_X(sexp_t* exp, qmem_t* qmem) {
  qid_t qid;
  qubit_t qubit;
  tangle_t* restrict tangle;
  assert( qmem );

  // move to the first argument
  exp = cdr(exp);
  if( !exp ) {
    printf("X-correction did not have any target qubit argument\n");
    exit(EXIT_FAILURE);
  }
  qid = get_qid( exp );

  if( cdr(exp) ) { 
    if( _verbose_ )
      printf(" (signal was: %d)\n", satisfy_signals( cdr(exp), qmem ));
    // there is a signal argument, bail out early if not satisfied
    if( satisfy_signals( cdr(exp), qmem ) == 0)
      return;
  }

  qubit = find_qubit( qid, qmem );
  if( invalid(qubit) ) {
    // create new qubit
    tangle = add_tangle( qid, qmem );
    qubit = find_qubit_in_tangle( qid, tangle );
  }
  qop_x( qubit );
}

void eval_Z(sexp_t* exp, qmem_t* qmem) {
  qid_t qid;
  qubit_t qubit;
  tangle_t* restrict tangle;
  assert( qmem );

  // move to the first argument
  exp = cdr(exp);
  if( !exp ) {
    printf("Z-correction did not have any target qubit argument\n");
    exit(EXIT_FAILURE);
  }
  qid = get_qid( exp );

  if( cdr(exp) ) { 
    // there is a signal argument, bail out early if not satisfied
    if( _verbose_ )
      printf(" (signal was: %d)\n", satisfy_signals( cdr(exp), qmem ));
    if( satisfy_signals( cdr(exp), qmem ) == 0 )
      return;
  }

  qubit = find_qubit( qid, qmem );
  if( invalid(qubit) ) {
    // create new qubit
    tangle = add_tangle( qid, qmem );
    qubit = find_qubit_in_tangle( qid, tangle );
  }
  qop_z( qubit );
}
 
// expects a list, evals the first argument and calls itself tail-recursively
void eval( sexp_t* restrict exp, qmem_t* restrict qmem ) {
  CSTRING* str = snew(0);
  sexp_t* command;
  sexp_t* rest;
  char opname;

  assert( qmem );

  if( exp == NULL )
    return;
  
  //ensure_list( exp );
  if( exp->ty == SEXP_LIST ) {
    command = car(exp);
    rest = cdr(exp);
    if( _verbose_ ) {
      print_sexp_cstr( &str, exp, STRING_SIZE );
      printf("evaluating %s\n", toCharPtr(str));
    }
  }
  else {
    assert( exp->ty == SEXP_VALUE );
    command = exp;
    rest = NULL;
    sexp_t tmp_list = (sexp_t){SEXP_LIST, NULL, 0, 0, exp, NULL, 0,
			       NULL, 0};
    //    print_sexp_cstr( &str, new_sexp_list(exp), STRING_SIZE );
    if( _verbose_ ) {
      print_sexp_cstr( &str, &tmp_list, STRING_SIZE );
      printf("evaluating %s\n", toCharPtr(str));
    }
  }
  
  opname = get_opname( command );

  sdestroy( str );
    
  switch ( opname ) {
  case 'E': 
    eval_E( command, qmem ); 
    if( _verbose_ )
      print_qmem(qmem);
    eval( rest, qmem ); 
    break;
  case 'M': 
    eval_M( command, qmem ); 
    if( _verbose_ )
      print_qmem(qmem);
    eval( rest, qmem ); 
    break;
  case 'X': 
    eval_X( command, qmem ); 
    if( _verbose_ )
      print_qmem(qmem);
    eval( rest, qmem ); 
    break;
  case 'Z': 
    eval_Z( command, qmem );
    if( _verbose_ )
      print_qmem(qmem);
    eval( rest, qmem ); 
    break;
  default: 
    printf("unknown command: %c\n", opname);
  }
}

amplitude parse_complex( const char* str ) {
  char* next_str = NULL;
  char* last_str = NULL;
  double real = strtod(str, &next_str);
  double imaginary = strtod(next_str, &last_str);
  if( last_str && (last_str[0] == 'i') )
    return real + imaginary * I;
  else 
    return real + 0 * I;
}


void parse_tangle( const sexp_t* exp, qmem_t* restrict qmem ) {
  qubit_t qubit;
  tangle_t* tangle = NULL;
  const sexp_t* qids_exp = exp->list;
  const sexp_t* amps_exp = exp->list->next;

  sexp_to_dotfile(exp,"input.dot");

  assert( qids_exp->ty == SEXP_LIST &&
	  amps_exp->ty == SEXP_LIST );

  const int num_amps = sexp_list_length(amps_exp);

  //at least one qid  
  sexp_t* qids = qids_exp->list;
  assert( qids && qids->val );
  for( sexp_t* qid_exp = qids; qid_exp; qid_exp = qid_exp->next) {
    qubit = find_qubit(get_qid(qid_exp), qmem);
    if( !invalid(qubit) ) {
      fprintf( stderr, 
	       "ERROR: trying to add already existing qubit "
	       "during input file initialization (qid:%d)\n",
	       qubit.qid );
      exit(EXIT_FAILURE);
    }
  }
  
  tangle = get_free_tangle(qmem);
 
  qmem->size += 1;
  tangle->size = sexp_list_length(qids_exp);
  tangle->qstate = NEW_QUANTUM_STATE( tangle->size );
  tangle->qids = add_qid( get_qid(qids), tangle->qids );

  while(qids->next) {
    qids=qids->next;
    append_qids( add_qid(get_qid(qids), NULL), tangle->qids );
  }

  sexp_t* amp = amps_exp->list;
  for( int i=0; i<num_amps ;  ++i ) {
    const pos_t index = atoi(amp->list->val);
    const amplitude a = parse_complex(amp->list->next->val);
    // init_quantum_state has (should) set all other amplitudes to 0
    tangle->qstate.vector[index] = a;
    amp=amp->next;
  } 
}

const tangle_t* fetch_first_tangle( const qmem_t* restrict qmem ) {

  for(int i=0; i<MAX_TANGLES; ++i) {
    const tangle_t* tangle = qmem->tangles[i];
    if( tangle )
      return tangle;
  }
  return NULL;
}


/* prints ONLY THE FIRST TANGLE in sexpr form to file, same format as input file, but
   also produces 0's */
void 
produce_output_file( const char* output_file, 
		     const qmem_t* restrict qmem ) {
  CSTRING* out = snew(STRING_SIZE); 
  char str[STRING_SIZE];
  //  int count=0;
  const tangle_t* restrict tangle = fetch_first_tangle(qmem);
  const quantum_state_t qstate = tangle->qstate;
  assert( tangle );
  assert( output_file );
  saddch(out,'(');
  // print qids
  saddch(out, '(');
  for( const qid_list_t* cons = tangle->qids;
       cons;
       cons=cons->rest ) {
    sprintf(str,"%d", cons->qid);
    sadd(out, str);
    if( cons->rest )
      saddch(out, ' ');
  }
  // end qids
  sadd(out, ")\n ");

  // print (basis amplitude)
  saddch(out, '(');
  for( unsigned int i=0; i<num_amplitudes(qstate); ++i ) {
    sprintf(str,"(%ui ", i);
    sadd(out, str);
    sprintf( str, "% .12g%+.12gi)", 
	     creal(qstate.vector[i]),
	     cimag(qstate.vector[i]) );
    sadd(out, str);
  }
  sadd(out, "\n  ");

  // end amplitudes
  saddch(out, ')');
  saddch(out, ')');
  saddch(out, '\n');
  FILE* file = fopen(output_file, "w");
  fputs(toCharPtr(out), file);
  fclose(file);
  sdestroy(out);
}

void initialize_input_state( const char* input_file, qmem_t* qmem ) {
  if( input_file == NULL )
    return;
  int fd = open( input_file, O_RDONLY );
  sexp_iowrap_t* input_port = init_iowrap( fd );
  sexp_t* exp = read_one_sexp(input_port);
  parse_tangle( exp, qmem );
  destroy_sexp( exp );
  destroy_iowrap( input_port );
  close(fd);  
}


int main(int argc, char* argv[]) {
  sexp_iowrap_t* input_port;
  sexp_t* mc_program;
  qmem_t qmem;
  CSTRING* str = snew( 0 );

  int silent = 0;
  char* output_file = NULL;
  int program_fd = 0;
  int c;
     
  opterr = 0;
  init_qmem( &qmem );
    
  while ((c = getopt (argc, argv, "isvmf:o::")) != -1)
    switch (c)
      {
      case 'i':
	_interactive_ = 1;
	break;
      case 's':
	silent = 1;
	break;
      case 'v':
	_verbose_ = 1;
	break;
      case 'm':
	_alt_measure_ = 1;
	break;
      case 'f':
	initialize_input_state(optarg, &qmem);
	break;
      case 'o':
	output_file = optarg;
	break;
      case '?':
	if (optopt == 'f')
	  fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	else if (optopt == 'o') {
	  output_file = "out";
	  break;
	}
	else if (isprint (optopt))
	  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	else
	  fprintf (stderr, 
		   "Unknown option character `\\x%x'.\n",
		   optopt);
	return 1;
      default:
	abort ();
      }
     
  //  

  if (_verbose_) {
    printf("Initial QMEM:\n ");
    print_qmem( &qmem );
  }
  if( _interactive_ ) {
    printf("Starting QVM in interactive mode.\n qvm> ");
    input_port = init_iowrap( 0 );  // we are going to read from stdin
    mc_program = read_one_sexp( input_port );
    while( mc_program ) {
      eval( mc_program->list, &qmem );
      print_qmem( &qmem );
      printf("\n qvm> ");
      destroy_sexp( mc_program );
      mc_program = read_one_sexp( input_port );
    }
  }
  else {
    // read input program
    program_fd = 
      optind < argc ?                // did the user pass a non-option argument?
      open(argv[optind], O_RDONLY) : // open the file
      0;                             // otherwise, use stdin
    input_port = init_iowrap( program_fd );
    if( input_port == NULL ) {
      printf(" ERROR while wrapping i/o, sexp_errno code is: %d\n", sexp_errno);
      abort();
    }
    mc_program = read_one_sexp( input_port );
    
    if (!silent) {
      print_sexp_cstr( &str, mc_program, STRING_SIZE );
      printf("I have read: \n%s\n", toCharPtr(str) );
    }
    // emit dot file
    /* sexp_to_dotfile( mc_program->list, "mc_program.dot" ); */
    
    eval( mc_program->list, &qmem );
    if( program_fd )
      close( program_fd );
  }

  //normalize at the end, not during measurement

  int tally=0;
  for( int t=0; tally<qmem.size; ++t ) {
    tangle_t* tangle = qmem.tangles[t];
    if( tangle ) {
      quantum_normalize( tangle->qstate );
      ++tally;
    }
  }
  
  if (!silent) {
    printf("Resulting quantum memory is:\n");
    print_qmem( &qmem );
  }

  if( output_file ) {
    produce_output_file(output_file, &qmem);
  }
  
  destroy_iowrap( input_port );
  sdestroy( str );
  destroy_sexp( mc_program );
  sexp_cleanup();
  free_qmem( &qmem );
  return 0;
}
