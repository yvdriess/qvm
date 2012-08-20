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

#include <sexp.h>
#include <sexp_ops.h>
#include <sexp_vis.h>
#include <quantum.h>


//#include "libquantum/config.h"

#include "bitmask.h"
#include "qvm.h"

#define STRING_SIZE (size_t)UCHAR_MAX	
#define MAX_TANGLES (size_t)SHRT_MAX
#define MAX_QUBITS  (size_t)SHRT_MAX

#define car hd_sexp
#define cdr next_sexp
#define max(x,y) x < y ? x : y


int _verbose_ = 0;
int _alt_measure_ = 0;
quantum_reg _proto_diag_qubit_;
quantum_reg _proto_dual_diag_qubit_;
  quantum_matrix _cz_gate_ = 
    { 4,4, (COMPLEX_FLOAT[16]){1,0,0,0,
			       0,1,0,0,
			       0,0,1,0,
			       0,0,0,-1} };

/************
 ** TANGLE **
 ************/
typedef struct qid_list {
  qid_t qid;
  struct qid_list* rest;
} qid_list_t;

typedef struct tangle { 
  tangle_size_t size;
  qid_list_t* qids;
  quantum_reg qureg;
 } tangle_t;  

tangle_t* init_tangle() {
  tangle_t* tangle = (tangle_t*) malloc(sizeof(tangle_t));   //ALLOC tangle
  tangle->size = 0;
  tangle->qids = NULL;
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
  tangle->size = 0;
  tangle->qids = NULL;
  quantum_delete_qureg( &tangle->qureg );
  free( tangle ); //FREE tangle
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

void print_tangle( const tangle_t* restrict tangle ) {
  assert( tangle );
  print_qids( tangle->qids );
  printf(" ,\n    {\n");
  if( tangle->qureg.size > 32 ) {
    printf("<a large quantum state>, really print? (y/N): ");
    /* if( getchar() == 'y' ) */
    /*   quantum_print_qureg( tangle->qureg ); */
  }
  else
    quantum_print_qureg( tangle->qureg );
  printf("}");
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

// use this function to return a correct qureg position
//  libquantum uses an reverse order (least significant == 0)
int get_target( const qubit_t qubit ) {
  return qubit.tangle->size - qubit.pos - 1;
}

bool invalid( const qubit_t qubit ) {
  return qubit.tangle == NULL;
}
quantum_reg* get_qureg( const qubit_t qubit ) {
  return &qubit.tangle->qureg;
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
  size_t size;
  signal_map_t signal_map;
  tangle_t* tangles[MAX_TANGLES];
} qmem_t;


void print_signal_map( const signal_map_t* restrict signal_map ) {
  printf(" {\n");
  for( int qid=0 ; qid<MAX_QUBITS ; ++qid ) {
    if( BITTEST(signal_map->entries, qid) )
      printf("  %d -> %d,\n", qid,
	     BITTEST(signal_map->signals, qid) ? 1 : 0 );
  }  
  printf(" }\n");
}

bool get_signal( const qid_t qid, 
		 const signal_map_t* restrict signal_map ) {
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
		 signal_map_t* restrict signal_map ) {
  if( BITTEST(signal_map->entries, qid) ) {
    printf( "ERROR: I was asked to set an already existing signal,\n\
  check quantum program correctness.\n");
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
  qid_list_t* qids = tangle->qids;
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
find_qubit(const qid_t qid, const qmem_t* restrict qmem) {
  tangle_t* tangle;
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

qid_list_t* add_qid( const qid_t qid, qid_list_t* restrict qids ) {
  // assuming qid is NOT already in qids
  // ALLOC QUBIT LIST
  qid_list_t* restrict new_qids = (qid_list_t*) malloc(sizeof(qid_list_t));
  new_qids->qid = qid;
  new_qids->rest = qids;
  return new_qids;
}

void append_qids( qid_list_t* new_qids, qid_list_t* target_qids ) {
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


void print_qmem( const qmem_t* restrict qmem ) {
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
  printf("signal map:");
  print_signal_map( &qmem->signal_map );
}

qmem_t* init_qmem() {
  qmem_t* restrict qmem = malloc(sizeof(qmem_t)); //ALLOC qmem

  qmem->size = 0;
  //qmem->tangles = calloc(MAX_TANGLES,sizeof(tangle_t*)); //ALLOC tangles
  //memset(&qmem->tangles, 0, sizeof(tangle_t*) * MAX_TANGLES);
  for( int i=0; i<MAX_TANGLES; ++i )
    qmem->tangles[i] = NULL;
  qmem->signal_map = (signal_map_t){{0},{0}};
  
  // instantiate prototypes (libquantum quregs)
  _proto_diag_qubit_ = quantum_new_qureg(0, 1);
  _proto_dual_diag_qubit_ = quantum_new_qureg(0, 2);
  quantum_hadamard(0, &_proto_diag_qubit_);
  quantum_hadamard(0, &_proto_dual_diag_qubit_);
  quantum_hadamard(1, &_proto_dual_diag_qubit_);
  quantum_gate2(0, 1, _cz_gate_, &_proto_dual_diag_qubit_);  

  // seed RNG
  //sranddev();
  srand(time(0));
  return qmem;
}

void free_qmem(qmem_t* qmem) {
  quantum_delete_qureg( &_proto_diag_qubit_ );
  quantum_delete_qureg( &_proto_dual_diag_qubit_ );

  for( int i=0, tally=0 ; tally < qmem->size ; i++ ) {
    assert(i<MAX_TANGLES);
    if( qmem->tangles[i] ) {
      free_tangle(qmem->tangles[i]);
      qmem->tangles[i] = NULL;
      ++tally;
    }
  }
  //free(qmem->tangles); //FREE tangles
  free(qmem); //FREE qmem
}

tangle_t* get_free_tangle(qmem_t* qmem) {
  tangle_t* restrict new_tangle = init_tangle();
  assert(new_tangle);
  // I loop here because tangles can get de-allocated (NULL-ed)
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
  quantum_copy_qureg(&_proto_dual_diag_qubit_,
		     &tangle->qureg);
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
  quantum_copy_qureg(&_proto_diag_qubit_,
		     &tangle->qureg);
  return tangle;
}


/* Adds new qubit BEHIND existing state: |q> x |+>
 */
void
add_qubit( const qid_t qid, 
	   tangle_t* restrict tangle) {
  assert(tangle);
  // appends new qid:  qids := [[qids...],qid]
  append_qids( add_qid(qid,NULL), tangle->qids );
  tangle->size += 1;
  // tensor |+> to tangle
  const quantum_reg new_qureg = 
    quantum_kronecker(&tangle->qureg,&_proto_diag_qubit_);
  // out with the old
  quantum_delete_qureg( &tangle->qureg );
  // in with the new
  tangle->qureg = new_qureg;
}

void 
delete_tangle( tangle_t* tangle,
	       qmem_t* restrict qmem ) {
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
delete_qubit(const qubit_t qubit, 
	     qmem_t* restrict qmem) {
  assert( !invalid(qubit) );
  tangle_t* tangle = qubit.tangle;
  // handle points to where current qid entry is stored
  qid_list_t** handle = &tangle->qids;
  qid_list_t* qids = tangle->qids;
  
  // listref by qubit.pos
  for( int i=0 ; i < qubit.pos ; ++i ) {
    handle = &qids->rest;
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

void 
merge_tangles(tangle_t* restrict tangle_1, 
	      tangle_t* restrict tangle_2, 
	      qmem_t* restrict qmem) {
  assert( tangle_1 && tangle_2 );
  tangle_1->size = tangle_1->size + tangle_2->size;
  // append qids of tangle_2 to tangle_1, destructively
  append_qids( tangle_2->qids, tangle_1->qids);
  // tensor both quregs
  const quantum_reg new_qureg = 
    quantum_kronecker( &tangle_1->qureg, &tangle_2->qureg );
  // out with the old
  quantum_delete_qureg( &tangle_1->qureg );
  quantum_delete_qureg( &tangle_2->qureg );
  // in with the new
  tangle_1->qureg = new_qureg;

  tangle_2->qids = NULL; // avoids the qid_list from being collected
  delete_tangle( tangle_2, qmem ); //free the tangle
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
  entry->name = name;
  entry->value = value;
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

static inline unsigned int
quantum_hash64(MAX_UNSIGNED key, int width)
{
  unsigned int k32;

  k32 = (key & 0xFFFFFFFF) ^ (key >> 32);

  k32 *= 0x9e370001UL;
  k32 = k32 >> (32-width);

  return k32;
}

static inline int
quantum_get_state(MAX_UNSIGNED a, quantum_reg reg)
{
  int i;

  if(!reg.hashw)
    return a;

  i = quantum_hash64(a, reg.hashw);

  while(reg.hash[i])
    {
      if(reg.node[reg.hash[i]-1].state == a)
	return reg.hash[i]-1;
      i++;
      if(i == (1 << reg.hashw))
	i = 0;
    }
  
  return -1;   
}
/* Add an element to the hash table */

static inline void
quantum_add_hash(MAX_UNSIGNED a, int pos, quantum_reg *reg)
{
  int i, mark = 0;

  i = quantum_hash64(a, reg->hashw);

  while(reg->hash[i])
    {
      i++;
      if(i == (1 << reg->hashw))
	{
	  if(!mark)
	    {
	      i = 0;
	      mark = 1;
	    }
	  else
	    quantum_error(QUANTUM_EHASHFULL);
	}
    }

  reg->hash[i] = pos+1;

}

/* Reconstruct hash table */

static inline void
quantum_reconstruct_hash(quantum_reg *reg)
{
  int i;

  /* Check whether register is sorted */

  if(!reg->hashw)
    return;
  
  for(i=0; i<(1 << reg->hashw); i++)
    reg->hash[i] = 0;
  for(i=0; i<reg->size; i++)
    quantum_add_hash(reg->node[i].state, i, reg);
}

int
quantum_diag_measure(int pos, double angle, quantum_reg* restrict reg)
{
  //int result=0;
  //int value=0;
  quantum_reg out;
  MAX_UNSIGNED pos2 = (MAX_UNSIGNED) 1 << pos;
  double limit = (1.0 / ((MAX_UNSIGNED) 1 << reg->width)) / 1000000;
  double prob=0, norm = 0;
  COMPLEX_FLOAT amp = 0;

  // TODO: currently just measures to <+_alpha|
  out.width = reg->width-1;
  out.size = reg->size;
  out.node = calloc(reg->size, sizeof(quantum_reg_node));
  //quantum_memman(size * sizeof(quantum_reg_node));
  out.hashw = reg->hashw;
  out.hash = reg->hash;
  
  for( int i=0 ; i<reg->size ; ++i ) {
    //    quantum_prob_inline( reg->node[i].ampl
  }

  if(reg->hashw)
    quantum_reconstruct_hash(reg);

  /* METHOD 1:
      loop through all collapsed basis and lookup the two contributing
      amplitudes.
      should have really rubbish cache usage
   */

  typedef unsigned int basis;
  basis upper_mask = ((basis)(-1/pos2))*pos2;
  basis lower_mask = -1 % pos2;
  assert( upper_mask + lower_mask == -1 );
  
  basis lpart,rpart;
  int free = 0;
  for(basis state=0; state<(1 << out.width); ++state ) {
    lpart = upper_mask & state<<1;
    rpart = lower_mask & state;
    basis k = lpart+rpart;
    int i = quantum_get_state(k, *reg);
    int j = quantum_get_state(k^pos2, *reg);
    int k_is_odd = k & pos2;
    if( i >= 0 )
      amp += k_is_odd ? -(reg->node[i].amplitude * quantum_cexp(-angle))
	              : reg->node[i].amplitude;
    if( j >= 0 )
      amp += k_is_odd ? reg->node[j].amplitude 
	              : -(reg->node[j].amplitude * quantum_cexp(-angle));
    if( i >= 0 || j >= 0 ) {
      prob = quantum_prob_inline( amp );
      if( prob > limit ) {
	assert(free<out.size);
	norm += prob;
	out.node[free].amplitude = amp;
	out.node[free].state = state;
	++free;
      }
      amp = 0;
    }
  }
  out.size = free;
  if( out.size != reg->size ) {
    out.node = realloc(out.node, (out.size)*sizeof(quantum_reg_node));
    if(out.node == NULL) 
      quantum_error(QUANTUM_ENOMEM);
  }
  // normalize, turned off
  /* norm = sqrt(norm); */
  /* if( abs(1-norm) > limit ) */
  /*   for( int i=0; i<out.size; ++i ) */
  /*     out.node[i].amplitude /= norm; */
  
  quantum_delete_qureg_hashpreserve(reg);
  *reg = out;
  return 1;


  /* METHOD 2: suggestion
      loop through all amplitudes
      
  */

}

void qop_cz( const qubit_t qubit_1, const qubit_t qubit_2 ) {
  const int tar1 = get_target(qubit_1);
  const int tar2 = get_target(qubit_2);
  assert( !(invalid(qubit_1) || invalid(qubit_2)) );
  assert( qubit_1.tangle == qubit_2.tangle );

  /* printf("Performing CZ on qubits %d and %d on tangle ",  */
  /* 	 qubit_1.qid, qubit_2.qid); */
  /* print_qids( qubit_1.tangle->qids ); */
  /* printf("\n"); */
  /* printf("  calling cz with targets %d and %d\n, ", tar1, tar2); */

  // manual cz because a) libquantum's gate2 appears to be bugggy and
  //  can be implemented optimally relatively easily, similar to cnot
  quantum_reg* reg = get_qureg( qubit_1 );
  MAX_UNSIGNED bitmask = 
    ((MAX_UNSIGNED) 1 << tar1) | ((MAX_UNSIGNED) 1 << tar2);
  for(int i=0; i<reg->size; i++)
    {
      /* Flip the target bit of a basis state if the control bit is set */     
      if((reg->node[i].state & bitmask) == bitmask)
	reg->node[i].amplitude *= (COMPLEX_FLOAT)-1;
    }
  //  quantum_gate2(tar1, tar2, _cz_gate_, get_qureg(qubit_1)); 
}

void qop_x( const qubit_t qubit ) {
  assert( !invalid(qubit) );
  quantum_sigma_x( get_target(qubit), get_qureg(qubit) );
}

void qop_z( const qubit_t qubit ) {
  assert( !invalid(qubit) );
  quantum_sigma_z( get_target(qubit), get_qureg(qubit) );
}

/* Apply a phase kick by the angle GAMMA */
void
quantum_inv_phase_kick(int target, double gamma, quantum_reg *reg)
{
  int i;
  COMPLEX_FLOAT z;

  z = quantum_conj(quantum_cexp(gamma));
  double* p = (double*)&z;
  printf("before phase kick (z=%f,%fi):\n",p[0],p[1]);
  //  quantum_print_qureg( *reg );
  
  
  for(i=0; i<reg->size; i++)
    {
      if(reg->node[i].state & ((MAX_UNSIGNED) 1 << target))
	reg->node[i].amplitude *= z;
    }

  printf("\nafter phase kick:\n");

  //  quantum_decohere(reg);
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
    else
      // add qid1 to qid2's tangle
      add_qubit( qid1, qubit_2.tangle );
  else
    if( invalid(qubit_2) )
      // add qid2 to qid1's tangle
      add_qubit( qid2, qubit_1.tangle );
    else
      if( qubit_1.tangle == qubit_2.tangle ) {
	// if not, qubit entries are already valid
	qop_cz( qubit_1, qubit_2 );
	return;
      }
      else
	// both tangles are non-NULL, merge both
	merge_tangles(qubit_1.tangle, qubit_2.tangle, qmem);
  // get valid qubit entries
  qubit_1 = find_qubit( qid1, qmem );
  qubit_2 = find_qubit( qid2, qmem );
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
  tangle_t* tangle;
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
    tangle = add_tangle( qid, qmem );
    qubit = find_qubit_in_tangle( qid, tangle );
  }
  // libquantum can only measure in ortho basis,
  //  but <+|q = <0|Hq makes it diagonal
  //  and <+_a| = <+|P_-a
  if( _verbose_ )
    printf("  measuring qubit %d on angle %2.4f\n", qid, angle);
  /* printf("   before + correction:\n"); */
  /* quantum_print_qureg( qubit.tangle->qureg ); */
  
  //  quantum_inv_phase_kick( get_target(qubit), angle, get_qureg(qubit) );

  if( _alt_measure_ )
    signal = quantum_diag_measure( get_target(qubit), 
				   angle,
				   get_qureg(qubit) );
  else {
    quantum_phase_kick( get_target(qubit), -angle, get_qureg( qubit ) );
  
    //printf("   after kick: \n");
    //  quantum_print_qureg( qubit.tangle->qureg );

    quantum_hadamard( get_target(qubit), get_qureg( qubit ) );
  
    //printf("   measuring : \n"     );
    // quantum_print_qureg( qubit.tangle->qureg );
    signal = quantum_bmeasure( get_target(qubit), get_qureg( qubit ) );

    //signal = quantum_diag_measure( get_target(qubit), angle, get_qureg(qubit) );
  }
    /* quantum_hadamard( get_target(qubit), get_qureg( qubit ) ); */
  /* quantum_phase_kick( get_target(qubit), angle, get_qureg( qubit ) ); */

  

  /* printf("   result is %d\n",signal); */
  set_signal( qid, signal, &qmem->signal_map );

  // remove measured qubit from memory
  delete_qubit( qubit, qmem );
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

COMPLEX_FLOAT parse_complex( const char* str ) {
  char* next_str = NULL;
  char* last_str = NULL;
  double real = strtod(str, &next_str);
  double imaginary = strtod(next_str, &last_str);
  if( last_str && (last_str[0] == 'i') )
    return real + imaginary * IMAGINARY;
  else 
    return real + 0 * IMAGINARY;
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
  tangle->qureg = quantum_new_qureg_size( num_amps, tangle->size  );
  tangle->qids = add_qid( get_qid(qids), tangle->qids );
  while(qids->next) {
    qids=qids->next;
    append_qids( add_qid(get_qid(qids), NULL), tangle->qids );
  }

  quantum_reg* reg = &tangle->qureg;
  sexp_t* amp = amps_exp->list;
  for( int i=0; i<num_amps ;  ++i ) {
    reg->node[i].state = atoi(amp->list->val);
    reg->node[i].amplitude = parse_complex(amp->list->next->val);
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
  quantum_reg reg;
  const tangle_t* tangle = fetch_first_tangle(qmem);
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
  reg = tangle->qureg;
  for( int i=0; i<reg.size; ++i ) {
    sprintf(str,"(%lli ", reg.node[i].state);
    sadd(out, str);
    sprintf( str, "% .12g%+.12gi)", 
	     quantum_real(reg.node[i].amplitude),
	     quantum_imag(reg.node[i].amplitude) );
    sadd(out, str);
    if( i+1<reg.size )
      sadd(out, "\n  ");
  }
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

void quantum_normalize( quantum_reg reg ) {
  double limit = 1.0e-8;
  COMPLEX_FLOAT norm=0;
  for( int i=0; i<reg.size; ++i ) {
    norm += quantum_prob_inline(reg.node[i].amplitude);
  }
  if( abs(1-norm) < limit )
    for( int i=0; i<reg.size; ++i ) {
      reg.node[i].amplitude /= norm;
    }
}

int main(int argc, char* argv[]) {
  sexp_iowrap_t* input_port;
  sexp_t* mc_program;
  qmem_t* restrict qmem = init_qmem();
  CSTRING* str = snew( 0 );

  int interactive = 0;
  int silent = 0;
  char* output_file = NULL;
  int program_fd;
  int c;
     
  opterr = 0;
    
  while ((c = getopt (argc, argv, "isvmf:o::")) != -1)
    switch (c)
      {
      case 'i':
	interactive = 1;
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
	initialize_input_state(optarg, qmem);
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
    print_qmem( qmem );
  }
  if( interactive ) {
    printf("Starting QVM in interactive mode.\n qvm> ");
    input_port = init_iowrap( 0 );  // we are going to read from stdin
    mc_program = read_one_sexp( input_port );
    while( mc_program ) {
      eval( mc_program->list, qmem );
      print_qmem( qmem );
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
    mc_program = read_one_sexp( input_port );
    if( program_fd )
      close( program_fd );
    
    if (!silent) {
      print_sexp_cstr( &str, mc_program, STRING_SIZE );
      printf("I have read: \n%s\n", toCharPtr(str) );
    }
    // emit dot file
    /* sexp_to_dotfile( mc_program->list, "mc_program.dot" ); */
    
    eval( mc_program->list, qmem );
  }

  //normalize at the end, not during measurement

  int tally=0;
  tangle_t* tangle=NULL;
  for( int t=0; tally<qmem->size; ++t ) {
    tangle = qmem->tangles[t];
    if( tangle ) {
      quantum_normalize( tangle->qureg );
      ++tally;
    }
  }
  
  if (!silent) {
    printf("Resulting quantum memory is:\n");
    print_qmem( qmem );
  }

  if( output_file ) {
    produce_output_file(output_file, qmem);
  }
  
  destroy_iowrap( input_port );
  sdestroy( str );
  destroy_sexp( mc_program );
  sexp_cleanup();
  free_qmem( qmem );
  return 0;

}

