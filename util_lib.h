#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define EPS 1e-6

typedef struct {
    int *arr;
    size_t curr_capacity;
    size_t max_capacity;
} Stack;

double **dmatrix(int num_row, int num_col);
double *dvector(int nh);
int **imatrix(int num_row, int num_col);
int *ivector(int nh);

void free_dvector(double *v);
void free_dmatrix(double **m,int num_row);
void free_ivector(int *v);
void free_imatrix(int **m,int num_row);

int tr2int (char *nt);
int nt2int (char nt);
int nt2int_rc (char nt);


int trinucleotide (char a, char b, char c);
double log2(double a);
void get_protein(char *dna, char *protein, int strand, int whole_genome);
void print_usage();
  
void print_viterbi(double **matr, int len_seq, int n_state, FILE *f);
void print_path(int **matr, int len_seq, int n_state, FILE *f);
int min(int lhs, int rhs);
int max(int lhs, int rhs);

void swap(int *lhs, int *rhs);

int is_equal(double a, double b);

void fprint_imatrix(int **matr, int nrow, int ncol, char *fname);
void fprint_dmatrix(double **matr, int nrow, int ncol, char *fname);
void fprint_ivector(int *v, int n, char *fname);
void print_ivector(int *v, int n);

void get_rc_dna(char *dna, char *dna1);
void get_rc_dna_indel(char *dna, char *dna1);

Stack create_stack();
void free_stack(Stack *st);
int top(Stack* st);
int pop(Stack *st);
void push(Stack* st, int x);
