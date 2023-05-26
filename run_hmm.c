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
} thread_data;

void* thread_func(void *threadarr);

int main (int argc, char **argv)
{
    clock_t start = clock();
    int i, j, c, max;
    HMM hmm;
    char *obs_seq, *obs_head;
    TRAIN train;
    int wholegenome;
    int format=0;
    FILE *fp_out, *fp_aa, *fp_dna, *fp;
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

    thread_data td;
    sprintf(tmp_str, "%s.out", out_header);
    td.out = fopen(tmp_str, "w");

    sprintf(tmp_str, "%s.faa", out_header);
    td.aa = fopen(tmp_str, "w");

    sprintf(tmp_str, "%s.ffn", out_header);
    td.dna = fopen(tmp_str, "w");

    td.hmm = &hmm;
    td.train = &train;
    td.wholegenome = wholegenome;
    td.format = format;

    fp = fopen(seq_file, "r");
    bp_count = 0;
    while (fgets(tmp_str, sizeof(tmp_str), fp)){
        if (tmp_str[0] == '>'){
            head_len = strlen(tmp_str);
        } else {
            bp_count += strlen(tmp_str);
        }
    }
    ++bp_count;
    td.obs_seq = (char*)malloc(bp_count * sizeof(char));
    td.obs_head = (char*)malloc(head_len);

    rewind(fp);
    while (fgets(tmp_str, sizeof(tmp_str), fp)){
        if (tmp_str[0] == '>'){
            strcpy(td.obs_head, tmp_str);
        } else {
            for (i = 0; i < strlen(tmp_str); ++i)
                if (!isalpha(tmp_str[i]))
                    tmp_str[i] = '\0';
            strcat(td.obs_seq, tmp_str);   
        }
    }

    td.cg = get_prob_from_cg(td.hmm, td.train, td.obs_seq);
    //printf("%s\n%s\n", td.obs_head, td.obs_seq);
    if (strlen(td.obs_seq) > 70) {
        //printf("%s\n%s\nWholegenome: %d\nCg: %d\nFormat:%d\n", td.obs_head, td.obs_seq, td.wholegenome, td.cg, td.format);
        viterbi(td.hmm, td.train, td.obs_seq, td.out, td.aa, td.dna, td.obs_head, td.wholegenome, td.cg, td.format);
    }

    free(td.obs_seq);
    free(td.obs_head);

    fclose(fp);
    fclose(td.out);
    fclose(td.aa);
    fclose(td.dna);
    clock_t end = clock();
    printf("Clock time used (by %d threads) = %.2f mins\n", threadnum, (end - start) / (60.0 * CLOCKS_PER_SEC));
}


void* thread_func(void *threadarr)
{
    thread_data *d;
    d = (thread_data*)threadarr;
    d->cg = get_prob_from_cg(d->hmm, d->train, d->obs_seq); //cg - 26 Ye April 16, 2016
    if (strlen(d->obs_seq)>70){
        viterbi(d->hmm, d->train, d->obs_seq, d->out, d->aa, d->dna, d->obs_head, d->wholegenome, d->cg, d->format);
    }
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

