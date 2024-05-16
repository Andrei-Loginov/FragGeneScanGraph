#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "util_lib.h"

#define STRINGLEN 4096

#define A 0
#define C 1
#define G 2
#define T 3

#define NUM_STATE 29

#define NOSTATE -1
#define S_STATE 0
#define E_STATE 1
#define R_STATE 2
#define S_STATE_1 3
#define E_STATE_1 4
#define M1_STATE 5
#define M2_STATE 6
#define M3_STATE 7
#define M4_STATE 8
#define M5_STATE 9
#define M6_STATE 10
#define M1_STATE_1 11
#define M2_STATE_1 12
#define M3_STATE_1 13
#define M4_STATE_1 14
#define M5_STATE_1 15
#define M6_STATE_1 16
#define I1_STATE 17
#define I2_STATE 18
#define I3_STATE 19
#define I4_STATE 20
#define I5_STATE 21
#define I6_STATE 22
#define I1_STATE_1 23
#define I2_STATE_1 24
#define I3_STATE_1 25
#define I4_STATE_1 26
#define I5_STATE_1 27
#define I6_STATE_1 28


#define TR_MM 0
#define TR_MI 1
#define TR_MD 2
#define TR_II 3
#define TR_IM 4
#define TR_DD 5
#define TR_DM 6
#define TR_GE 7
#define TR_GG 8
#define TR_ER 9
#define TR_RS 10
#define TR_RR 11
#define TR_ES 12
#define TR_ES1 13

#define M_GROUP 0
#define I_GROUP 1
#define R_GROUP 2
#define E_GROUP 3
#define S_GROUP 4
#define M_GROUP_1 5
#define I_GROUP_1 6
#define E_GROUP_1 7
#define S_GROUP_1 8



typedef struct {

  double  pi[29];    /* pi[1..N] pi[i] is the initial state distribution. */
  int N;               /* number of state */

  double tr[14];                 /* transition probability from a (delete/insert/match) state to a state */

  double e_M_1[6][16][4];      /* transition probability from a lowest-level state  to a  lowest-level state*/
  double e_M[6][16][4];

  double tr_R_R[4][4];
  double tr_I_I[4][4];
  double tr_M_I[4][4];

  double tr_S[61][64];
  double tr_E[61][64];
  double tr_S_1[61][64];
  double tr_E_1[61][64];

  double S_dist[6];  /*sigma, mu,alpha, sigma_r, mu_r, alpha_r */
  double E_dist[6];
  double S1_dist[6];
  double E1_dist[6];
} HMM;



typedef struct {

  double trans[44][6][16][4];
  double rtrans[44][6][16][4];
  double noncoding[44][4][4];
  double start[44][61][64];
  double stop[44][61][64];
  double start1[44][61][64];
  double stop1[44][61][64];

  double S_dist[44][6];
  double E_dist[44][6];
  double S1_dist[44][6];
  double E1_dist[44][6];

} TRAIN;


typedef struct {
    double **alpha;
    int **path;
    char* O;
    int len_seq;
    int *temp_i;
    int *temp_i_1;
    int first_column_prev[NUM_STATE + 1];
    int curr_column_prev[NUM_STATE + 1];
    int tmp_curr_column_prev[NUM_STATE + 1];
    int calculated_flg;
} ViterbiResult;

typedef struct {
    size_t n_edge;
    size_t overlap;
    int** adjacency_matrix;
    int* seq_len;
    int* dead_end_flg;
    char** obs_seq;
    char** head;
    int *ind;
    int *order;
    ViterbiResult *edge_results;
} Graph;

typedef struct edge {
    int num;
    ViterbiResult res;
    int next_count;
    Graph *g;
} Edge;

typedef struct tmp_result {
    double alpha;
    double alpha2;
    int path;
    int prev_ind;
} TmpResult;

typedef struct graph_path {
    int *vpath;
    char *O;
    double *alpha;
    int seq_len;
} GraphPath;

int get_prob_from_cg(HMM *hmm, TRAIN *train, char *O); //return cg - 26 Ye April 16, 2016
int get_prob_form_cg_graph(HMM *hmm_ptr, TRAIN *train_ptr, Graph *g);
void get_train_from_file(char *filename, HMM *hmm_ptr, char *mfilename, char *mfilename1, char *nfilename, char *sfilename,char *pfilename,char *s1filename,char *p1filename, char *dfilename, TRAIN *train_ptr);
ViterbiResult viterbi(HMM *hmm_ptr, char *O, int whole_genome, ViterbiResult *prev_result, char* head);
void backtrack(HMM *hmm_ptr, TRAIN *train_ptr, FILE *fp_out, FILE *fp_aa, FILE *fp_dna,char *head, int whole_genome, int cg, int format,
               ViterbiResult *viterbi_result);

void backtrack_graph_path(HMM *hmm_ptr, TRAIN *train_ptr, FILE *fp_out, FILE *fp_aa, FILE *fp_dna,char *head, int whole_genome, int cg, int format,
               GraphPath *gp); //gene prediction using restored path

ViterbiResult viterbi_edge(HMM *hmm_ptr, Graph *g, size_t edge_index, int whole_genome);
void viterbi_graph_dag(HMM *hmm_ptr, Graph* g, int whole_genome);
GraphPath restore_path(ViterbiResult *res, Graph *g, int start, int num_state); //restores path using ViterbiMatrix


//Optimization refactoring
int state2group(int i);

double any_state_prob(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int *prev_ind, int n_prev, size_t overlap, int whole_genome);

TmpResult match_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap, int whole_genome, int to);
TmpResult match_state1_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap, int whole_genome, int to);
TmpResult insertion_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap, int to);
TmpResult insertion_state1_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap, int to);
TmpResult non_coding_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap, int to);
TmpResult end_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap);
TmpResult end_state1_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap);
TmpResult start_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap);
TmpResult start_state1_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap);
//


int count_from2(int t, char *O, int seq_len, char *prev_O, int prev_seq_len);

void free_hmm(HMM *hmm);
void free_ViterbiResult(ViterbiResult *viterbi_result);
void get_protein(char *dna, char *protein, int strand, int whole_genome);
void get_rc_dna(char *dna, char *dna1);
void get_corrected_dna(char *dna, char *dna_f);


char complementary_nucleotide(char c);
size_t get_edge_num_m(Graph *g, int ind, char strand);
size_t get_edge_num(int* ind, size_t n_edge, int num, char strand);
Graph read_gfa(FILE *f);

Graph read_graph(FILE* fp, FILE* fp_matr);
int *topological_sort(int **adj_matrix, int n_vert);
void dfs(int **adj_matrix, int n_vert, int *used, int *ans, int *curr_len, int ind);
void free_graph(Graph* g);
