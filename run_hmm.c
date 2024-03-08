#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "hmm.h"

#include <pthread.h>

#define ADD_LEN 1024
#define STRINGLEN 4096

#define graph_check
#define HMM_read_check

typedef struct thread_data
{
    FILE *out;
    FILE *aa;
    FILE *dna;
    char *obs_head;
    char *obs_seq;
    int wholegenome;
    int cg;
    int format;
    HMM *hmm;
    TRAIN *train;
    double **alpha;
    int **path;
} thread_data;


int main (int argc, char **argv)
{
    clock_t start = clock();
    int i, j, c, max;
    HMM hmm;
    char **obs_seq, **obs_head;
    TRAIN train;
    Graph g;
    int wholegenome;
    int format=0;
    FILE *fp_out, *fp_aa, *fp_dna, *fp, *fp_matr;
    char hmm_file[STRINGLEN] = "";
    char out_header[STRINGLEN] = "";
    char aa_file[STRINGLEN] = "";
    char seq_file[STRINGLEN] = "";
    char out_file[STRINGLEN] = "";
    char dna_file[STRINGLEN] = "";
    char train_file[STRINGLEN] = "";
    char mstate_file[STRINGLEN] = "";
    char rstate_file[STRINGLEN] = "";
    char nstate_file[STRINGLEN] = "";
    char sstate_file[STRINGLEN] = "";
    char pstate_file[STRINGLEN] = "";
    char s1state_file[STRINGLEN] = "";     /* stop codon of gene in - stand */
    char p1state_file[STRINGLEN] = "";
    char dstate_file[STRINGLEN] = "";
    char train_dir[STRINGLEN] = "";
    int count=0;
    int currcount = 0;
    int total = 0;
    char tmp_str[STRINGLEN] = "";
    int *obs_seq_len;
    int bp_count;  /* count the length of each line in input file */
    int head_len;
    int threadnum = 1;
    int rc;
    int edge_num;
    int cg_count;

    thread_data *threadarr;
    char **lastline, **currline;

    strncpy(train_dir, argv[0], strlen(argv[0])-12);
    strcat(train_dir, "train/");
    strcpy(mstate_file, train_dir);
    strcat(mstate_file, "gene");
    strcpy(rstate_file, train_dir);
    strcat(rstate_file, "rgene");
    strcpy(nstate_file, train_dir);
    strcat(nstate_file, "noncoding");
    strcpy(sstate_file, train_dir);
    strcat(sstate_file, "start");
    strcpy(pstate_file, train_dir);
    strcat(pstate_file, "stop");
    strcpy(s1state_file, train_dir);
    strcat(s1state_file, "stop1");
    strcpy(p1state_file, train_dir);
    strcat(p1state_file, "start1");
    strcpy(dstate_file, train_dir);
    strcat(dstate_file, "pwm");


    /* read command line argument */
    if (argc <= 8){
        fprintf(stderr, "ERROR: You missed some parameters for input\n");
        print_usage();
        exit(EXIT_FAILURE);
    }

    while ((c=getopt(argc, argv, "fs:o:w:t:p:")) != -1){
        switch (c){
        case 's':
            strcpy(seq_file, optarg);
            if (access(seq_file, F_OK)==-1){
                fprintf(stderr, "ERROR: Sequence file [%s] does not exist\n", seq_file);
                print_usage();
                exit(EXIT_FAILURE);
            }
        break;
        case 'w':
            wholegenome = atoi(optarg);
            if (wholegenome != 0 && wholegenome != 1){
                fprintf(stderr, "ERROR: An incorrect value for the option -w was entered\n");
                print_usage();
                exit(EXIT_FAILURE);
            }
        break;
        case 'p':
            threadnum = atoi(optarg);
            if (threadnum < 1){
                fprintf(stderr, "ERROR: An incorrect value [%d] for the option -p was entered\n", threadnum);
                print_usage();
                exit(EXIT_FAILURE);
            }
            printf("Using %d threads.\n", threadnum);
        break;
        case 'o':
            strcpy(out_header, optarg);
            break;
        case 't':
            strcpy(train_file, optarg);
            strcpy(hmm_file, train_dir);
            strcat(hmm_file, train_file);

            if (access(hmm_file, F_OK)==-1){
                fprintf(stderr, "ERROR: The file for model parameters [%s] does not exist\n", hmm_file);
                print_usage();
                exit(EXIT_FAILURE);
            }
        break;
        case 'f':
            format = 1;
        break;
        }
    }


    /* check whether the specified files exist */
    if (access(mstate_file, F_OK)==-1){
        fprintf(stderr, "Forward prob. file [%s] does not exist\n", mstate_file);
        exit(1);
    }
    if (access(rstate_file, F_OK)==-1){
        fprintf(stderr, "Backward prob. file [%s] does not exist\n", rstate_file);
        exit(1);
    }
    if (access(nstate_file, F_OK)==-1){
        fprintf(stderr, "noncoding prob. file [%s] does not exist\n", nstate_file);
        exit(1);
    }
    if (access(sstate_file, F_OK)==-1){
        fprintf(stderr, "start prob. file [%s] does not exist\n", sstate_file);
        exit(1);
    }
    if (access(pstate_file, F_OK)==-1){
        fprintf(stderr, "stop prob. file [%s] does not exist\n", pstate_file);
        exit(1);
    }
    if (access(s1state_file, F_OK)==-1){
        fprintf(stderr, "start1 prob. file [%s] does not exist\n", s1state_file);
        exit(1);
    }
    if (access(p1state_file, F_OK)==-1){
        fprintf(stderr, "stop1 prob. file [%s] does not exist\n", p1state_file);
        exit(1);
    }
    if (access(dstate_file, F_OK)==-1){
        fprintf(stderr, "pwm dist. file [%s] does not exist\n", dstate_file);
        exit(1);
    }
    if (access(hmm_file, F_OK)==-1){
        fprintf(stderr, "hmm file [%s] does not exist\n", hmm_file);
        exit(1);
    }
    threadnum = 1; //crutch to use only 1 thread;
    /* read all initial model */
    hmm.N=NUM_STATE;
    get_train_from_file(hmm_file, &hmm, mstate_file, rstate_file, nstate_file, sstate_file, pstate_file,s1state_file, p1state_file, dstate_file, &train);

#ifdef HMM_read_check
    FILE *fdbg = fopen("../bin/hmm_e_M_new.txt", "w");
    int k;
    for (i = 0; i < 6; ++i){
        fprintf(fdbg, "M%d matrix\n", i);
        for (j = 0; j < 16; ++j){
            for (k = 0; k < 4; ++k){
                fprintf(fdbg, "%lf ", hmm.e_M[i][j][k]);
            }
            fprintf(fdbg, "\n");
        }
    }
    fclose(fdbg);
#endif

    sprintf(tmp_str, "%s.out", out_header);
    fp_out = fopen(tmp_str, "w");

    if (!fp_out) {
        printf("Can't open out file %s\n", tmp_str);
        exit(0);
    }

    sprintf(tmp_str, "%s.faa", out_header);
    fp_aa = fopen(tmp_str, "w");

    sprintf(tmp_str, "%s.ffn", out_header);
    fp_dna = fopen(tmp_str, "w");


    fp = fopen(seq_file, "r");
    if (!fp) {
        printf("Can't open seq_file %s\n", seq_file);
        exit(0);
    }

    sprintf(tmp_str, "%s-matrix.txt", seq_file);
    fp_matr = fopen(tmp_str, "r");
    edge_num = 0;

    if (!fp_matr) {
        printf("Can't open matrix file %s\n", tmp_str);
        exit(0);
    }
    g = read_graph(fp, fp_matr);
    //ViterbiResult* results = (ViterbiResult*)malloc(g.n_edge * sizeof(ViterbiResult));
    cg_count = get_prob_form_cg_graph(&hmm, &train, &g);
    viterbi_graph(&hmm, &g, 0, wholegenome);
    fprint_imatrix(g.adjacency_matrix, g.n_edge, g.n_edge, "some_out_file.txt");
    GraphPath ans = restore_path(g.edge_results, &g, 1, NUM_STATE);
    printf("%s\n", ans.O);
    print_ivector(ans.vpath, ans.seq_len);

    free_graph(&g);

    fclose(fp);
    fclose(fp_out);
    fclose(fp_aa);
    fclose(fp_dna);
    clock_t end = clock();
    printf("Clock time used (by %d threads) = %.2f mins\n", threadnum, (end - start) / (60.0 * CLOCKS_PER_SEC));
}




int appendSeq(char *input, char **seq, int input_max)
{
    int len, inputlen, max;
    char *tmp;

    max = input_max;
    if (*seq != NULL)
    {
        len = strlen(*seq);
    }
    else
    {
        len = 0;
    }
    inputlen = strlen(input);
    if ((len + inputlen) >= max)
    {
        while ((len + inputlen) >= max)
        {
            max += ADD_LEN;
        }
        tmp = (char*)malloc(sizeof(char) * max);
        memset(tmp, '\0', sizeof(char) * max);
        if (*seq != NULL)
        {
            memcpy(tmp, *seq, len);
        }
        free(*seq);
        *seq = tmp;
    }
    strcat(*seq, input);
    return max;
}

