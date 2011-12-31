#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sexp.h>
#include <sexp_ops.h>
#include <sexp_vis.h>
#include <quantum.h>

#define STRING_SIZE sizeof(char)
#define MAX_TANGLES ((size_t)32)
#define MAX_QUBITS ((size_t)32)
#define SIGNAL_MAP_SIZE (sizeof(qid_t)*8)

#define car hd_sexp
#define cdr next_sexp
#define max(x,y) x < y ? x : y

extern void quantum_copy_qureg(quantum_reg *src, quantum_reg *dst);

typedef unsigned char qid_t;
typedef unsigned char tangle_size_t;
typedef unsigned char pos_t;

typedef struct qid_list {
  qid_t qid;
  struct qid_list* rest;
} qid_list_t;

quantum_reg _proto_diag_qubit_;
quantum_reg _proto_dual_diag_qubit_;
  quantum_matrix _cz_gate_ = 
    { 4,4, (COMPLEX_FLOAT[16]){1,0,0,0,
			       0,1,0,0,
			       0,0,1,0,
			       0,0,0,-1} };


typedef struct tangle { 
  tangle_size_t size;
  qid_list_t* qids;
  quantum_reg qureg;
 } tangle_t;  

typedef struct qubit {
  tangle_t* tangle;
  qid_t qid;
  pos_t pos;
} qubit_t;

typedef struct signal_map {
  // two bitfields
  //  entries : if qid has an entry, not needed for correct programs
  //  signals : value of the signal
  unsigned int entries : SIGNAL_MAP_SIZE;
  unsigned int signals : SIGNAL_MAP_SIZE;
} signal_map_t;

typedef struct qmem {
  size_t size;
  signal_map_t signal_map;
  tangle_t** tangles;
} qmem_t;

qubit_t _invalid_qubit_ = { NULL, -1, -1 };

bool invalid( const qubit_t qubit ) {
  return qubit.tangle == NULL;
}
quantum_reg* get_qureg( const qubit_t qubit ) {
  return &qubit.tangle->qureg;
}

tangle_t* init_tangle() {
  tangle_t* tangle = (tangle_t*) malloc(sizeof(tangle_t));   //ALLOC tangle
  tangle->size = 0;
  tangle->qids = NULL;
  return tangle;
}

bool get_signal( const qid_t qid, 
		 const signal_map_t* restrict signal_map ) {
  if( (signal_map->entries) & (1 << qid) )
    return (signal_map->signals) & (1 << qid);
  else {
    printf( "ERROR: I was asked a signal map entry (qid:%d) that wasn't there,\n\
  check quantum program correctness.\n",qid);
    exit(EXIT_FAILURE);
  }

}
void set_signal( const qid_t qid, 
		 const bool signal, 
		 signal_map_t* restrict signal_map ) {
  if( (signal_map->entries) & (1 << qid) ) {
    printf( "ERROR: I was asked to set an already existing signal,\n\
  check quantum program correctness.\n");
    exit(EXIT_FAILURE);
  }
  signal_map->entries |= (1 << qid);
  signal_map->signals |= (1 << qid);
}

qubit_t 
find_qubit_in_tangle( qid_t qid, 
		      const tangle_t* tangle )
{
  qid_list_t* qids;
   if( tangle ) {
     qids = tangle->qids;
     for( int i=0;
	  qids;
	  ++i, qids=qids->rest ) {
      if( qids->qid == qid ) 
	return (qubit_t){ (tangle_t*)tangle, qid, i };
    }
  }
  return _invalid_qubit_;
}

qubit_t 
find_qubit(qid_t qid, qmem_t* qmem) {
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

qid_list_t* add_qid( qid_t qid, qid_list_t* qids ) {
  // assuming qid is NOT already in qids
  // ALLOC QUBIT LIST
  qid_list_t* new_qids = (qid_list_t*) malloc(sizeof(qid_list_t*));
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
  free( tangle ); //FREE tangle
  tangle = NULL;
}

void print_tangle( const tangle_t* restrict tangle ) {
  assert( tangle );
  printf("[");
  for( qid_list_t* qids = tangle->qids;
       qids;
       qids=qids->rest ) {
    printf("%d", qids->qid);
    if( qids->rest )
      printf(", ");
  }
  printf("] {");
  quantum_print_qureg( tangle->qureg );
  printf("}");
}

void print_signal_map( const signal_map_t* restrict signal_map ) {
  printf(" {\n");
  for( int qid=0 ; qid<SIGNAL_MAP_SIZE ; ++qid ) {
    if( signal_map->entries & (1<<qid) )
      printf("  %d -> %d,\n", qid,
	     signal_map->entries & (1<<qid) ? 1 : 0 );
  }  
  printf(" }\n");
}

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
  _proto_diag_qubit_ = quantum_new_qureg(0, 1);
  _proto_dual_diag_qubit_ = quantum_new_qureg(0, 2);
  quantum_hadamard(0, &_proto_diag_qubit_);
  quantum_hadamard(0, &_proto_dual_diag_qubit_);
  quantum_hadamard(1, &_proto_dual_diag_qubit_);
  quantum_gate2(0, 1, _cz_gate_, &_proto_dual_diag_qubit_);  

  quantum_print_qureg(_proto_diag_qubit_);
  quantum_print_qureg(_proto_dual_diag_qubit_);

  qmem_t* qmem = malloc(sizeof(qmem_t)); //ALLOC qmem
  qmem->size = 0;
  qmem->tangles = calloc(MAX_TANGLES,sizeof(tangle_t*)); //ALLOC tangles

  qmem->signal_map = (signal_map_t){0,0};

  srand(time(0));

  return qmem;
}

void free_qmem(qmem_t* qmem) {
  for( int i=0, tally=0 ; tally < qmem->size ; i++ ) {
    assert(i<MAX_TANGLES);
    if( qmem->tangles[i] ) {
      free_tangle((tangle_t*)qmem->tangles[i]);
      ++tally;
    }
  }
  free(qmem->tangles); //FREE tangles
  free(qmem); //FREE qmem
}

tangle_t* get_free_tangle(qmem_t* qmem) {
  // I loop here because tangles can get de-allocated (NULL-ed)
  printf("    searching for free tangle\n");
  for(int i=0; i<MAX_TANGLES; ++i) {
    if( qmem->tangles[i] == NULL ) {
      printf("      found one at %d!\n",i);
      qmem->tangles[i] = init_tangle();
      return qmem->tangles[i];
    }
    printf("      ...skipping %d\n",i);
  }
  printf("Ran out of qmem memory, too many tangles! (> %lu)\n",MAX_TANGLES);
  exit(EXIT_FAILURE);
}

void
add_tangle(const qid_t qid1, 
	   const qid_t qid2, 
	   qmem_t* restrict qmem) {
  // allocate new tangle in qmem
  tangle_t*  restrict tangle = get_free_tangle(qmem);

  // init tangle
  tangle->qids = add_qid( qid1, tangle->qids );
  tangle->qids = add_qid( qid2, tangle->qids );
  tangle->size = 2;

  // update qmem info
  qmem->size += 1;

  // init quantum state
  quantum_copy_qureg(&_proto_dual_diag_qubit_,
		     &tangle->qureg);
}

void add_qubit(const qid_t qid, tangle_t* tangle) {
  quantum_reg* qureg = &tangle->qureg;
  tangle->qids = add_qid( qid, tangle->qids );
  tangle->size += 1;
  // tensor |+> to tangle
  tangle->qureg = quantum_kronecker(&_proto_diag_qubit_,qureg);
}

void delete_qubit(qubit_t qubit) {
  assert( !invalid(qubit) );
  tangle_t* restrict tangle = qubit.tangle;
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

  free(qids); // FREE QUBIT LIST
  // OPTIONAL: dealloc tangle (requires qmem pointer)
}

void 
merge_tangles(tangle_t* tangle_1, 
	      tangle_t* tangle_2, 
	      qmem_t* qmem) {
  assert( tangle_1 && tangle_2 );
  tangle_1->size = tangle_1->size + tangle_2->size;
  // append tangle_2 to tangle_1, destructively
  append_qids( tangle_2->qids, tangle_1->qids);
  free( tangle_2 ); //free the tangle, but not the qid_list
}

void ensure_list( sexp_t* exp ) {
  CSTRING* str = NULL;
 
  if( exp->ty != SEXP_LIST ) {
    str = snew( STRING_SIZE );
    print_sexp_cstr( &str, exp, STRING_SIZE );
    printf("ERROR: malformed expression, expecting a list expression and\
  got:\n  %s", toCharPtr( str ));
    exit(EXIT_FAILURE);
  }
}

void ensure_value( sexp_t* exp ) {
  CSTRING* str = NULL;

  if ( exp->ty != SEXP_VALUE ) {
    str = snew( STRING_SIZE );
    print_sexp_cstr( &str, exp, STRING_SIZE );
    printf("ERROR: malformed expression, expecting a value expression and\
  got:\n  %s", toCharPtr( str ));
    exit(EXIT_FAILURE);
  }
}

char get_opname( sexp_t* exp ) {
  return exp->val[0];
}

int get_qid( sexp_t* exp ) {
  return atoi( exp->val );
}

double get_angle( sexp_t* exp ) {
  return atof( exp->val );
}

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
    if( invalid(qubit_2) )
      // if both unknown, create new tangle with two |+> states
      add_tangle(qid1, qid2, qmem);
    else
      // add qid1 to qid2's tangle
      add_qubit( qid1, qubit_2.tangle );
  else
    if( invalid(qubit_2) )
      // add qid2 to qid1's tangle
      add_qubit( qid2, qubit_1.tangle );
    else
      if( qubit_1.tangle == qubit_2.tangle )
	// qubit entries are already valid
	goto perform_cz; 
      else
	// both tangles are non-NULL, merge both
	merge_tangles(qubit_1.tangle, qubit_2.tangle, qmem);

  // get valid qubit entries
  qubit_1 = find_qubit( qid1, qmem );
  qubit_2 = find_qubit( qid2, qmem );
    
  // perform CZ on resulting tangle for both qubits
 perform_cz:
  quantum_gate2(qubit_1.pos,
		qubit_2.pos,
		_cz_gate_, 
		get_qureg(qubit_1));
}

void eval_M(sexp_t* exp, qmem_t* qmem) {
  int qid;
  double angle = 0.0;
  tangle_t* tangle;
  const tangle_t* new_tangle;
  bool signal;
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
    angle = get_angle( exp );
  } 

  // TODO support for signals

  //  printf("  Measuring qubits %d\n",qid);

  qubit_t qubit = find_qubit( qid, qmem );
  if( invalid(qubit) ) {
    printf("ERROR: Measurement got unknown target qubit %d\n",qid);
    print_qmem( qmem );
    exit(EXIT_FAILURE);
  }
  // libquantum can only measure in ortho basis,
  //  but <+|q = <0|Hq makes it diagonal
  // TODO: take angle into account too!
  quantum_sigma_x( qubit.pos, get_qureg( qubit ) );
  signal = quantum_bmeasure( qubit.pos, get_qureg( qubit ) );
  set_signal( qid, signal, &qmem->signal_map );

  // remove measured qubit from memory
  delete_qubit( qubit );
}

bool satisfy_signals( const sexp_t* restrict exp, 
		      const qmem_t* qmem) {
  qid_t qid;
  bool signal;
  // TODO: support complete signal syntax
  //  currently only supports single signal dependency
  if( (exp->ty != SEXP_LIST) || strcmp(exp->list->val, "q") ) {
    
    printf("PARSE ERROR: signal did not start with (q ...) in %s\n",
	   exp->list->val);
    exit(EXIT_FAILURE);
  }
  else {
    qid = atoi( exp->list->next->val );
    signal = get_signal( qid, &qmem->signal_map );
    printf("qid %d's signal was: %d\n",qid,signal);
    return signal;
  }
  
  
}

void eval_X(sexp_t* exp, qmem_t* qmem) {
  qid_t qid;
  qubit_t qubit;

  assert( qmem );

  // move to the first argument
  exp = cdr(exp);
  if( !exp ) {
    printf("X-correction did not have any target qubit argument\n");
    exit(EXIT_FAILURE);
  }
  qid = get_qid( exp );

  if( cdr(exp) ) { 
    // there is a signal argument, bail out early if not satisfied
    if( !satisfy_signals( cdr(exp), qmem ) )
      return;
  }

  qubit = find_qubit( qid, qmem );
  if( invalid(qubit) ) {
    printf("ERROR: X-correction got an unknown target qubit %d\n",qid);
    exit(EXIT_FAILURE);
  }
  quantum_sigma_x( qubit.pos, get_qureg(qubit) );
}

void eval_Z(sexp_t* exp, qmem_t* qmem) {
  qid_t qid;
  qubit_t qubit;
  
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
    if( !satisfy_signals( cdr(exp), qmem ) )
      return;
  }

  qubit = find_qubit( qid, qmem );
  if( invalid(qubit) ) {
    printf("ERROR: Z-correction got an unknown target qubit %d\n",qid);
    exit(EXIT_FAILURE);
  }
  quantum_sigma_z( qubit.pos, get_qureg(qubit) );
}
 
// expects a list, evals the first argument and calls itself tail-recursively
void eval( sexp_t* restrict exp, qmem_t* restrict qmem ) {
  sexp_t* command;
  sexp_t* rest;
  char opname;

  assert( qmem );

  if( exp == NULL )
    return;
  
  ensure_list( exp );
  command = car(exp);
  rest = cdr(exp);
  opname = get_opname( command );

  switch ( opname ) {
  case 'E': 
    printf("evaluating E ...\n"); 
    eval_E( command, qmem ); 
    eval( rest, qmem ); 
    break;
  case 'M': 
    printf("evaluating M ...\n"); 
    eval_M( command, qmem ); 
    eval( rest, qmem ); 
    break;
  case 'X': 
    printf("evaluating X ...\n"); 
    eval_X( command, qmem ); 
    eval( rest, qmem ); 
    break;
  case 'Z': 
    printf("evaluating Z ...\n"); 
    eval_Z( command, qmem );
    eval( rest, qmem ); 
    break;
  default: 
    printf("unknown command: %c\n", opname);
  }
}


int main(int argc, char* argv[]) {
  sexp_iowrap_t* input_port;
  sexp_t* mc_program;
  qmem_t* qmem;
  //  CSTRING* str = NULL;

  input_port = init_iowrap( 0 );  // we are going to read from stdin

  if (!input_port) { 
    printf("IO ERROR\n");
    return -1; 
  }

  // read input program
  mc_program = read_one_sexp( input_port );

  // test print
  /* str = snew( STRING_SIZE ); */
  /* print_sexp_cstr( &str, mc_program, STRING_SIZE ); */
  /* printf("I have read: \n%s\n", toCharPtr(str) ); */

  // emit dot file
  /* sexp_to_dotfile( mc_program->list, "mc_program.dot" ); */

  qmem = init_qmem();
  
  eval( mc_program->list, qmem );

  printf("Resulting quantum memory is:\n");
  print_qmem( qmem );

  free_qmem( qmem );
  printf("ok\n");
  return 0;

}

