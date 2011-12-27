#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>

#include <sexp.h>
#include <sexp_ops.h>
#include <sexp_vis.h>

#define STRING_SIZE sizeof(char)
#define MIN_TANGLE_SIZE 8
#define MAX_TANGLES sizeof(qid_t)

#define car hd_sexp
#define cdr next_sexp
#define max(x,y) x < y ? x : y

typedef unsigned char qid_t;
typedef unsigned char tangle_size_t;

typedef struct qid_list {
  qid_t qid;
  struct qid_list* rest;
} qid_list_t;

typedef struct tangle { 
  tangle_size_t size;
  qid_list_t* qids;
  //TODO: amp vec here
 } tangle_t;

void init_tangle( tangle_t* tangle) {
  assert( tangle == NULL );

  tangle = (tangle_t*) malloc(sizeof(tangle_t));   //ALLOC tangle
  tangle->size = 0;
  tangle->qids = NULL;
}

bool find_qid( qid_t qid, qid_list_t* restrict qids ) {
  assert( qids );
  if( qids->qid == qid )
    return true;
  else
    if( qids->rest )
      return find_qid( qid, qids->rest );
    else 
      return false;
}

void add_qid( qid_t qid, qid_list_t* qids ) {
  if( qids ) // follow list until end
    add_qid( qid, qids->rest );
  else {
    qids = (qid_list_t*) malloc(sizeof(qid_list_t*));
    qids->qid = qid;
    qids->rest = NULL;
  }
}

void append_qids( qid_list_t* new_qids, qid_list_t* target_qids ) {
  assert( new_qids && target_qids );
  while( target_qids->rest ) {
    target_qids = target_qids->rest;
  }
  target_qids->rest = new_qids;
  
}

void remove_qid( qid_t qid, qid_list_t* qids ) {

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
  free( tangle ); //FREE tangle
  tangle = NULL;
}

typedef struct qmem {
  size_t size;
  tangle_t** tangles;
} qmem_t;

void init_qmem(qmem_t* qmem) {
  qmem = malloc(sizeof(qmem_t)); //ALLOC qmem
  qmem->size = 0;
  qmem->tangles = calloc(MAX_TANGLES,sizeof(tangle_t*)); //ALLOC tangles
}

void free_qmem(qmem_t* qmem) {
  for(int i=0; i < qmem->size; ++i)
    free_tangle((tangle_t*)qmem->tangles[i]);  
  free(qmem->tangles); //FREE tangles
  free(qmem); //FREE qmem
  qmem = NULL;
}

tangle_t* get_tangle(const qid_t qid, const qmem_t* qmem) {
  tangle_t* tangle = NULL;
  for(int i=0; i<qmem->size; ++i) {
    tangle = qmem->tangles[i];
    if( find_qid(qid, tangle->qids) )
      return tangle;      
  }
  return NULL;
}

tangle_t* get_free_tangle(qmem_t* qmem) {
  // I loop here because tangles can get de-allocated (NULL-ed)
  for(int i=0; i<qmem->size; ++i) {
    if( qmem->tangles[i] == NULL )
      return qmem->tangles[i];
  }
  printf("Ran out of qmem memory, too many tangles!");
  exit(EXIT_FAILURE);
}

tangle_t* add_tangle(const qid_t qid1, 
		     const qid_t qid2, 
		     qmem_t* qmem) {
  tangle_t* tangle = get_free_tangle(qmem);
  init_tangle(tangle);
  tangle->size = 2;
  add_qid( qid1, tangle->qids );
  add_qid( qid2, tangle->qids );
  // init quantum state
  return tangle;
}

tangle_t* add_qubit(const qid_t qid, tangle_t* tangle) {
  add_qid( qid, tangle->qids );
  tangle->size += 1;
  // tensor |+> to tangle
  return tangle;
}

tangle_t* merge_tangles(tangle_t* tangle_1, 
			tangle_t* tangle_2, 
			qmem_t* qmem) {
  assert( tangle_1 && tangle_2 );
  tangle_1->size = tangle_1->size + tangle_2->size;
  // append tangle_2 to tangle_1, destructively
  append_qids( tangle_2->qids, tangle_1->qids);
  free( tangle_2 ); //free the tangle, but not the qid_list
  return tangle_1;
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
  return strtol( exp->val, NULL, 10 );
}

void eval_E(sexp_t* exp, qmem_t* qmem) {
  int qid1, qid2;
  tangle_t* tangle_1;
  tangle_t* tangle_2;
  const tangle_t* new_tangle;
 
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

  printf("  Entangling qubits %d and %d\n",qid1,qid2);
  printf("    grabbing tangle for qubit %d\n", qid1);

  // get tangle for qid1
  tangle_1 = get_tangle( qid1, qmem );

  printf("    grabbing tangle for qubit %d\n", qid1);
  // get tangle for qid2
  tangle_2 = get_tangle( qid2, qmem );

  if (tangle_1 && tangle_2) {
    // if both NULL, create new tangle with two |+> states
    new_tangle = add_tangle(qid1, qid2, qmem);
  }
  else {
    if ((tangle_1==NULL) || (tangle_2==NULL)) {
      // if either NULL
      if (tangle_1) {
	// add qid2 to tangle_1
	new_tangle = add_qubit( qid2, tangle_1 );
      }
      else {
	// add qid1 to tangle_2
	new_tangle = add_qubit( qid1, tangle_2 );
      }
    }
    else{
      // both tangles are non-NULL, merge both
      new_tangle = merge_tangles(tangle_1, tangle_2, qmem);
    }
  }
  // perform CZ on resulting tangle for both qubits
}

// expects a list, evals the first argument and calls itself tail-recursively
void eval( sexp_t* exp, qmem_t* qmem ) {
  sexp_t* command;
  sexp_t* rest;
  char opname;

  if( exp == NULL )
    return;
  
  ensure_list( exp );
  command = car(exp);
  rest = cdr(exp);
  opname = get_opname( command );

  switch ( opname ) {
  case 'E': printf("evaluating E ...\n"); eval_E( command, qmem ); eval( rest, qmem ); break;
  case 'M': printf("evaluating M ...\n"); eval( rest, qmem ); break;
  case 'X': printf("evaluating X ...\n"); eval( rest, qmem ); break;
  case 'Z': printf("evaluating Z ...\n"); eval( rest, qmem ); break;
  default: printf("unknown command: %c\n", opname);
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

  init_qmem( qmem );
  
  eval( mc_program->list, qmem );

  free_qmem( qmem );
  printf("ok\n");
  return 0;

}

