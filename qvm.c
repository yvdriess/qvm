#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include <ccan/jmap/jmap.h>

#include <sexp.h>
#include <sexp_ops.h>
#include <sexp_vis.h>

#define STRING_SIZE 512

#define car hd_sexp
#define cdr next_sexp

typedef struct tangle { 
  int size=0;
 } tangle_t;

struct qubit_map {
  int qid;
  tangle_t* tangle;
};

typedef struct qmem {
  int foo;
} qmem_t;

#define init_qmem(qmem) qmem = (qmem_t){ 42 }

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
  int qid1;
  int qid2;
 
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
  return;
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
  qmem_t qmem;
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
  
  eval( mc_program->list, &qmem );


  printf("ok\n");
  return 0;

}

