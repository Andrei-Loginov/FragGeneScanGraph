#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <ctype.h>
#include "hmm.h"
#define viterbi_out_flg
//#define I_state_debug
//#define I1_state_debug
//#define E_state_debug
//#define M_state_profiling
//#define profiling
//#define M_state_debug
//#define R_state_debug
//#define match_vertex
//#define I1_debug
//#define R_state_debug
//#define E1_state_debug
//#define M_rc_debug
//#define M_inf
#define crash67
//#define printf_head_indices
//#define fprint_71_prev
#define print_graph_after_reading
//#define E1_issue
//#define non_coding_issue
#define non_utf8_issue


void dump_memory(void *p, int size);


ViterbiResult viterbi_edge(HMM *hmm_ptr, Graph *g, size_t edge_index, int whole_genome) {
    ViterbiResult ans;
    ans.O = (char*)malloc((g->seq_len[edge_index] + 1) * sizeof (char));
    ans.O = strcpy(ans.O, g->obs_seq[edge_index]);
    ans.len_seq = g->seq_len[edge_index];
    ans.alpha = dmatrix(NUM_STATE, ans.len_seq);
    ans.path = imatrix(NUM_STATE, ans.len_seq);
    ans.temp_i = (int*)malloc(6 * sizeof(int));
    ans.temp_i_1 = (int*)malloc(6 * sizeof(int));
    int i, t;
    for (i = 0; i < 6; ++i){
        ans.temp_i[i] = 0;
        ans.temp_i_1[i] = 0;
    }

    for (i = -1; i < NUM_STATE; ++i) {
        ans.curr_column_prev[i + 1] = -1;
        ans.first_column_prev[i + 1] = -1;
        ans.tmp_curr_column_prev[i + 1] = -1;
    }
    size_t n_prev= 0 ;
    for (i = 0; i < g->n_edge; ++i){
        if (g->adjacency_matrix[i][edge_index] && g->edge_results[i].calculated_flg)
            ++n_prev;
    }
    int *prev_indices = ivector(n_prev);
    size_t counter = 0;
    g->dead_end_flg[edge_index] = 1;
    for (i = 0; i < g->n_edge; ++i){
        if (g->adjacency_matrix[i][edge_index] && g->edge_results[i].calculated_flg){
            prev_indices[counter] = i;
            ++counter;
        }
        if (g->adjacency_matrix[edge_index][i]){
            g->dead_end_flg[edge_index] = 0;
        }
    }
    int rbound = ans.len_seq;
    if (!g->dead_end_flg[edge_index]){
        rbound -= g->overlap;
    }

    for (t = 0; t < rbound; ++t){
        //printf("edge = %d\tt = %d; ", edge_index, t);
        for (i = 0; i < NUM_STATE; ++i){
#ifdef crash67
            char *curr_head = g->head[edge_index];
#endif
            ans.alpha[i][t] = any_state_prob(hmm_ptr, t, i, &ans, g->edge_results, prev_indices, n_prev, g->overlap, whole_genome);
        }

        for (i = -1; i < NUM_STATE; ++i){
            ans.curr_column_prev[i + 1] = ans.tmp_curr_column_prev[i + 1];
        }
    }
    free_ivector(prev_indices);
///
///Print viterbi matrix and path matrix to files
///
#ifdef viterbi_out_flg
    char fname[4096];
    sprintf(fname, "../run_result/with_graph/multiple_edge/%s-new-matrix.csv", g->head[edge_index]);
    FILE *f = fopen(fname, "w");
    if (!f) {
        printf("The file was not opened\n");
    } else
        print_viterbi(ans.alpha, ans.len_seq, NUM_STATE, f);
    fclose(f);

    sprintf(fname, "../run_result/with_graph/multiple_edge/%s-new-path.csv", g->head[edge_index]);
    FILE *f_path = fopen(fname, "w");
    if (!f_path){
        printf("The file was not opened\n");
    } else {
        print_path(ans.path, ans.len_seq, NUM_STATE, f_path);
    }
    fclose(f_path);
#endif
#ifdef crash67
            char *curr_head = g->head[edge_index];
#endif
    return ans;
}

void viterbi_graph_dag(HMM *hmm_ptr, Graph* g, int whole_genome) {
    //assuming that adjacency matrix is upper triangular
    g->order = topological_sort(g->adjacency_matrix, g->n_edge);
    g->edge_results = (ViterbiResult*)malloc(g->n_edge * sizeof (ViterbiResult));
    size_t i;
    for (i = 0; i < g->n_edge; ++i) {
        g->edge_results[i].calculated_flg = 0;
    }
    for (i = 0; i < g->n_edge; ++i){
        g->edge_results[g->order[i]] = viterbi_edge(hmm_ptr, g, g->order[i], whole_genome);
        g->edge_results[g->order[i]].calculated_flg = 1;
#ifdef crash67
        ViterbiResult current_res = g->edge_results[g->order[i]];
#endif
    }
}

GraphPath restore_path(ViterbiResult *res, Graph *g, int start, int num_state){
    /*Find end of the most likely path*/
    int i, j, t, curr_edge, vertex_path;
    GraphPath ans;
    ans.seq_len = 0;// = res[start].len_seq;
    int **vpath, *edge_path, edge_counter = 0, curr_len_seq;
    double **alpha;
    curr_edge = start;
    vpath = (int**)malloc(g->n_edge * sizeof(int*));
    edge_path = (int*)malloc(g->n_edge * sizeof(int));
    alpha = (double**)malloc(g->n_edge * sizeof(double*));
#ifdef crash67
    // /home/andrei/Documents/Projects/Masters diploma/run_result/with_graph/multiple_edge
    char fpath_name[10000];
    sprintf(fpath_name, "/home/andrei/Documents/Projects/Masters diploma/run_result/with_graph/multiple_edge/trajectory/%s_ep.txt", g->head[start]);
    FILE *f_ep = fopen(fpath_name, "w");
    fclose(f_ep);
#endif

    while (curr_edge != -1){
        g->used_backtrack[curr_edge] = 1;

#ifdef crash67
        FILE *f_app = fopen(fpath_name, "a");
        fprintf(f_app, "%s\n", g->head[curr_edge]);
        fclose(f_app);
#endif
        int fully_calc_flg = is_fully_computed(res->alpha, NUM_STATE, res->len_seq, g->overlap);
        ans.seq_len += (g->dead_end_flg[curr_edge]) ? res[curr_edge].len_seq : res[curr_edge].len_seq - g->overlap;
        vpath[edge_counter] = (int*)malloc(res[curr_edge].len_seq * sizeof(int));
        alpha[edge_counter] = (double*)malloc(res[curr_edge].len_seq * sizeof(double));
        edge_path[edge_counter] = curr_edge;
        curr_len_seq = (g->dead_end_flg[curr_edge]) ? res[curr_edge].len_seq : res[curr_edge].len_seq - g->overlap;

        if (edge_counter == 0){ //the last edge
            vpath[edge_counter][curr_len_seq - 1] = 0;
            for (i = 0; i < num_state; ++i){
                if (res[curr_edge].alpha[i][curr_len_seq - 1] < res[curr_edge].alpha[vpath[edge_counter][curr_len_seq - 1]][curr_len_seq - 1]) {
                    vpath[edge_counter][curr_len_seq - 1] = i;
                }
            }
        } else {
            vpath[edge_counter][curr_len_seq - 1] = vertex_path;
        }
        int prev_vpath, curr_vpath;
        //char *head = g->head[curr_edge];
        //double prev_alpha;
        alpha[edge_counter][curr_len_seq - 1] = res[curr_edge].alpha[vpath[edge_counter][curr_len_seq - 1]][curr_len_seq - 1];
        for (t = curr_len_seq - 2; t >= 0; --t){
            //prev_vpath = curr_vpath;
            vpath[edge_counter][t] = res[curr_edge].path[vpath[edge_counter][t + 1]][t + 1];
            curr_vpath = vpath[edge_counter][t];
            alpha[edge_counter][t] = res[curr_edge].alpha[vpath[edge_counter][t]][t];
        }
        vertex_path = res[curr_edge].path[vpath[edge_counter][0]][0];
        curr_edge = res[curr_edge].curr_column_prev[vpath[edge_counter][curr_len_seq - 1] + 1];
        //printf("New curr_edge = %d\n", curr_edge);
        ++edge_counter;
    }
    ans.vpath = (int*)malloc(ans.seq_len * sizeof(int));
    ans.alpha = (double*)malloc(ans.seq_len * sizeof(double));
    ans.O = (char*)malloc(ans.seq_len + 1);
    ans.O[ans.seq_len] = '\0';
    int counter = 0;
    for (i = edge_counter - 1; i >= 0; --i){
        for (j = 0; j < (g->dead_end_flg[edge_path[i]] ? res[edge_path[i]].len_seq : res[edge_path[i]].len_seq - g->overlap); ++j, ++counter){
            ans.vpath[counter] = vpath[i][j];
            ans.alpha[counter] = alpha[i][j];
            ans.O[counter] = res[edge_path[i]].O[j];
        }
    }
    for (i = 0; i < edge_counter; ++i){
        free(vpath[i]);
        free(alpha[i]);
    }
    free(vpath);
    free(alpha);
    free(edge_path);
    printf("%s\tSeq_len = %zu\n", g->head[start], strlen(ans.O));
    return ans;
}

void backtrack_graph_path(HMM *hmm_ptr, TRAIN *train_ptr, FILE *fp_out, FILE *fp_aa, FILE *fp_dna,char *head, int whole_genome, int cg, int format,
                          GraphPath *gp){
    char *head_short=NULL;
    char delimi[] = " ";
    double max_dbl = 10000000000.0;
    double prob = max_dbl, final_score;
    double *alpha = gp->alpha;
    //int **path = viterbi_result->path;
    int i, j, t, kk, len_seq = gp->seq_len;
    char *O = gp->O;
    int *vpath = gp->vpath;
    int print_save, codon_start;
    int start_t, dna_start_t;

    char dna_tmp[300000];
    char dna[300010];
    char dna1[300010];
    char dna_f[300010];
    char dna_f1[300010];
    char protein[100000];
    int dna_id=0, dna_f_id=0;

    int insert[100];
    int delete[100];
    int insert_id, delete_id;

    int start_orf, prev_match, frame;
    int end_t, temp_t;
    int out_nt;

    int gene_len, refine = 0;
    if (whole_genome==1){
        gene_len = 120;
        refine = 1;
    } else {
        gene_len = 60;
    }

    head_short = strtok(head, delimi);
    fprintf(fp_out, "%s\n", head_short); //use head_short, Ye, April 22, 2016



    print_save = 0;
    codon_start=0;
    start_t=-1;

    int glen = strlen(O);
    char codon[4], utr[65];
    for (t=0; t<len_seq; t++){

        if (codon_start==0 && start_t < 0 &&
           ((vpath[t]>=M1_STATE && vpath[t]<=M6_STATE) ||
           (vpath[t]>=M1_STATE_1 && vpath[t]<=M6_STATE_1) ||
           vpath[t] == S_STATE || vpath[t] == S_STATE_1 ))
        {
            start_t=t+1;
            //printf("Note assign start_t %d (t=%d vpath %d)\n", start_t, t, vpath[t]);
        }

        if (codon_start==0 &&
           (vpath[t]==M1_STATE || vpath[t]==M4_STATE ||
           vpath[t]==M1_STATE_1 || vpath[t]==M4_STATE_1)) {

            memset(dna,0,300000);
            memset(dna1,0,300000);
            memset(dna_f,0,300000);
            memset(dna_f1,0,300000);
            memset(protein,0, 100000);
            memset(insert,0,100);
            memset(delete,0,100);

            insert_id = 0;
            delete_id = 0;
            dna_id = 0;
            dna_f_id = 0;
            dna[dna_id] = O[t];
            dna_start_t = t + 1; //Ye April 21, 2016
            dna_f[dna_f_id] = O[t];
            //printf("Note start dna: t = %d, dna_id %d, dna_f_id %d, add %c\n", t, dna_id, dna_f_id, O[t]);
            start_orf= t + 1;
            prev_match = vpath[t];

            if (vpath[t] < M6_STATE){
                codon_start=1;
            }else{
                codon_start=-1;
            }

        } else if (codon_start!=0 && (vpath[t]==E_STATE || vpath[t]==E_STATE_1 || t==len_seq-1)){

            if (vpath[t]==E_STATE || vpath[t]==E_STATE_1){
                end_t = t + 3;
            } else {
                end_t = t + 1;

                /* FGS1.12 start: remove incomplete codon */
                temp_t = t;
                while(vpath[temp_t] != M1_STATE && vpath[temp_t] != M4_STATE  &&
                        vpath[temp_t] != M1_STATE_1  && vpath[temp_t] != M4_STATE_1){

                    dna_f[dna_f_id] = '\0';
                    dna_f_id--;

                    dna[dna_id] = '\0';
                    dna_id--;

                    temp_t--;
                }
            /* FGS1.12 end: remove incomplete codon */
            }
            final_score = (alpha[end_t-4]- alpha[start_t+2] ) / (end_t-start_t-5);
            frame = start_orf % 3;
            if (frame == 0){
                frame = 3;
            }

            if (dna_id > gene_len  ){
                if (codon_start == 1){
                    if(start_t == dna_start_t - 3) { //add complete start codon to dna, Ye April 21, 2016
                        strcpy(dna_tmp, dna);
                        sprintf(dna, "%c%c%c%s", O[start_t - 1], O[start_t], O[start_t + 1], dna_tmp);
                        //printf("add start codon to dna: %c%c%c\n", O[start_t-1], O[start_t], O[start_t+1]);
                        //printf("old dna %d %s\n", strlen(dna_tmp), dna_tmp);
                        //printf("new dna %d %s\n", strlen(dna), dna);
                    }
                    if( refine ) { //add refinement of the start codons here, Ye, April 16, 2016
                        int start_old = start_t;
                        codon[0] = 0;
                        strncpy(codon, O + start_old - 1, 3);
                        codon[3] = 0;
                        int s = 0;
                        //find the optimal start codon within 30bp up- and downstream of start codon
                        double e_save;
                        int s_save;
                        while( ( ! (!strcmp(codon, "TAA") || !strcmp(codon, "TAG") || !strcmp(codon, "TGA")) ) && (start_old-1-s-35>=0)) {
                            if(!strcmp(codon, "ATG") || !strcmp(codon, "GTG") || !strcmp(codon, "TTG")) {
                                utr[0] = 0;
                                strncpy(utr, O+start_old-1-s-30,63);
                                utr[63] = 0;
                                //printf("check s=%d, codon %s\n", s, codon);
                                double freq_sum = 0;
                                for(j = 0; j < strlen(utr) - 2; j ++) {
                                    int idx = trinucleotide(utr[j], utr[j+1], utr[j+2]);
                                    freq_sum -= train_ptr->start[cg][j][idx];
                                    //printf("j=%d, key=%c%c%c %d, start %lf\n", j, utr[j], utr[j+1], utr[j+2], idx, train_ptr->start[cg][j][idx]);
                                }
                                if(s == 0) {
                                    e_save = freq_sum;
                                    s_save = s;
                                }
                                else if(freq_sum < e_save) {
                                    e_save = freq_sum;
                                    s_save = -1 * s;
                                }
                                //printf("s=%d freq_sum %lf\n", s, freq_sum);
                                //getchar();
                            }
                            s += 3;
                            codon[0] = 0;
                            strncpy(codon, O+start_old-1-s, 3);
                            codon[3] = 0;
                        }
                        start_t = start_old+s_save;
                        //update dna
                        if(s_save != 0) {
                            //printf("start refined + %d -> %d\n", start_old, start_t);
                            dna[0] = 0;
                            strncpy(dna, O + start_t - 1, end_t - start_t + 1);
                            dna[end_t - start_t + 1] = 0;
                        }
                    }
                    if (check_dna_symbols(dna)) {
                        fprintf(fp_out, "%d\t%d\t+\t%d\t%lf\t", start_t, end_t, frame, final_score);
                        fprintf(fp_out, "I:");
                        for (i=0; i<insert_id; i++){
                            fprintf(fp_out, "%d,", insert[i]);
                        }
                        fprintf(fp_out, "\tD:");
                        for (i=0; i<delete_id; i++){
                            fprintf(fp_out, "%d,", delete[i]);
                        }
                        fprintf(fp_out, "\n");

                        fprintf(fp_aa, "%s_%d_%d_+\n", head_short, start_t, end_t);
                        fprintf(fp_dna, "%s_%d_%d_+\n", head_short, start_t, end_t);

                        get_protein(dna,protein,1, whole_genome);
                        fprintf(fp_aa, "%s\n", protein);
                        if (format==0){
                            fprintf(fp_dna, "%s\n", dna);
                        }else if (format==1){
                            fprintf(fp_dna, "%s\n", dna_f);
                        }
#ifdef non_utf8_issue
                        if(!check_dna_symbols(dna) || !check_dna_symbols(dna_f)){
                            printf("Bad gene\n");
                        }
#endif
                    }

                } else if (codon_start==-1){
                    //printf("reverse strand dna-len %d, start_t %d, dna_start %d, add-up %d, end_t %d\n", strlen(dna), start_t, dna_start_t, dna_start_t + strlen(dna), end_t);
                    //getchar();
                    if(dna_start_t + strlen(dna) == end_t - 2) { //add complete start codon (on reverse strand) to dna, Ye April 21, 2016
                        strcpy(dna_tmp, dna);
                        sprintf(dna, "%s%c%c%c", dna_tmp, O[end_t-3], O[end_t-2], O[end_t-1]);
                        //printf("add start codon on the reverse strand to dna: %c%c%c\n", O[end_t-3], O[end_t-2], O[end_t-1]);
                    }
                    if(refine) { //add refinement of the start codons here, Ye, April 16, 2016
                        int end_old = end_t; //reverse
                        codon[0] = 0;
                        strncpy(codon, O + end_t-1-2, 3);
                        codon[3] = 0;
                        int s = 0;
                        //find the optimal start codon within 30bp up- and downstream of start codon
                        double e_save;
                        int s_save;
                        while((!(!strcmp(codon, "TTA") || !strcmp(codon, "CTA") || !strcmp(codon, "TCA"))) && (end_old-2+s+35 < glen)) {
                            if(!strcmp(codon, "CAT") || !strcmp(codon, "CAC") || !strcmp(codon, "CAA")) {
                                utr[0] = 0;
                                strncpy(utr, O+end_old-1-2+s-30,63);
                                utr[63] = 0;
                                //printf("check s=%d, codon %s\n", s, codon);
                                double freq_sum = 0;
                                for(j = 0; j < strlen(utr) - 2; j ++) {
                                    int idx = trinucleotide(utr[j], utr[j+1], utr[j+2]);
                                    freq_sum -= train_ptr->stop1[cg][j][idx]; //stop1?? Ye, April 18, 2016
                                    //printf("j=%d, key=%c%c%c %d, stop1 %lf\n", j, utr[j], utr[j+1], utr[j+2], idx, train_ptr->stop1[cg][j][idx]);
                                }
                                if(s == 0) {
                                    e_save = freq_sum;
                                    s_save = s;
                                } else if( freq_sum < e_save ) {
                                    e_save = freq_sum;
                                    s_save = s;
                                }
                                //printf("s=%d freq_sum %lf\n", s, freq_sum);
                                //getchar();
                            }
                            s += 3;
                            codon[0] = 0;
                            strncpy(codon, O+end_old-1-2+s, 3);
                            codon[3] = 0;
                        }
                        end_t = end_old+s_save;
                        //update dna
                        if(s_save != 0) {
                            //printf("start refined - end %d -> %d\n", end_old, end_t);
                            dna[0] = 0;
                            strncpy(dna, O + start_t - 1, end_t - start_t + 1);
                            dna[end_t - start_t + 1] = 0;
                        }
                    }

                    fprintf(fp_out, "%d\t%d\t-\t%d\t%lf\t", start_t, end_t, frame, final_score);
                    fprintf(fp_out, "I:");
                    for (i=0; i<insert_id; i++){
                        fprintf(fp_out, "%d,", insert[i]);
                    }
                    fprintf(fp_out, "\tD:");
                    for (i=0; i<delete_id; i++){
                        fprintf(fp_out, "%d,", delete[i]);
                    }
                    fprintf(fp_out, "\n");

                    fprintf(fp_aa, "%s_%d_%d_-\n", head_short, start_t, end_t);
                    fprintf(fp_dna, "%s_%d_%d_-\n", head_short, start_t, end_t);

                    get_protein(dna,protein,-1, whole_genome);
                    get_rc_dna(dna, dna1);
                    get_rc_dna_indel(dna_f, dna_f1);
                    fprintf(fp_aa, "%s\n", protein);
                    if (format==0){
                        fprintf(fp_dna, "%s\n", dna1);
                    }else if (format==1){
                        fprintf(fp_dna, "%s\n", dna_f1);
                    }
                }
            }
            codon_start=0;
            start_t = -1;
            end_t = -1;
            dna_id = 0;
            dna_f_id = 0;

        } else if (codon_start!=0 &&
          ((vpath[t]>=M1_STATE && vpath[t]<=M6_STATE) ||
           (vpath[t]>=M1_STATE_1 && vpath[t]<=M6_STATE_1)) &&
           vpath[t]-prev_match<6) {

            if (vpath[t] < prev_match){
                out_nt = vpath[t]+6-prev_match;
            } else {
                out_nt = vpath[t]-prev_match;
            }
            for (kk=0; kk<out_nt; kk++){   /* for deleted nt in reads */
                dna_id ++;
                dna[dna_id] = 'N';
                //printf("dna_id %d, dna-len %d\n", dna_id, strlen(dna));
                dna_f_id ++;
                dna_f[dna_f_id] = 'x';
                if (kk>0){
                    delete[delete_id]=t+1;
                    delete_id++;
                }
            }
            dna[dna_id]=O[t];
            //printf("dna_id %d, add %d %c dna-len %d\n", dna_id, t, O[t], strlen(dna));
            dna_f[dna_f_id]=O[t];
            prev_match = vpath[t];

        } else if (codon_start!=0 &&
          ((vpath[t]>=I1_STATE && vpath[t]<=I6_STATE) ||
           (vpath[t]>=I1_STATE_1 && vpath[t]<=I6_STATE_1))) {

            dna_f_id ++;
            dna_f[dna_f_id] = tolower(O[t]);
            insert[insert_id]=t+1;
            insert_id++;

        } else if (codon_start!=0 && vpath[t]==R_STATE) {
            /* for long NNNNNNNNN, pretend R state */
            codon_start=0;
            start_t=-1;
            end_t = -1;
            dna_id=0;
            dna_f_id=0;

        }
    }
    free_ivector(vpath);
    return;
}

void backtrack(HMM *hmm_ptr, TRAIN *train_ptr, FILE *fp_out, FILE *fp_aa, FILE *fp_dna,char *head, int whole_genome, int cg, int format,
               ViterbiResult *viterbi_result){
    char *head_short=NULL;
    char delimi[] = " ";
    double max_dbl = 10000000000.0;
    double prob = max_dbl, final_score;
    double **alpha = viterbi_result->alpha;
    int **path = viterbi_result->path;
    int i, j, t, kk, len_seq = strlen(viterbi_result->O);
    char* O = viterbi_result->O;
    int *  vpath = (int *)ivector(len_seq);
    int print_save, codon_start;
    int start_t, dna_start_t;

    char dna_tmp[300000];
    char dna[300010];
    char dna1[300010];
    char dna_f[300010];
    char dna_f1[300010];
    char protein[100000];
    int dna_id=0, dna_f_id=0;

    int insert[100];
    int delete[100];
    int insert_id, delete_id;

    int start_orf, prev_match, frame;
    int end_t, temp_t;
    int out_nt;

    int gene_len, refine = 0;
    if (whole_genome==1){
        gene_len = 120;
        refine = 1;
    } else {
        gene_len = 60;
    }

    head_short = strtok(head, delimi);
    fprintf(fp_out, "%s\n", head_short); //use head_short, Ye, April 22, 2016

    /* find the state for O[N] with the highest probability */
    prob = max_dbl;
    //������� ����� ����, ������� ���� ������������ �������������.
    for (i = 0; i < hmm_ptr->N; i++){
        if (alpha[i][len_seq-1] < prob){
            prob = alpha[i][len_seq-1];
            vpath[len_seq-1] = i;
        }
    }

    /* backtrack the optimal path */
    //��������������� ���� ����� ����, ������� ������������ ������������ �������������
    for(t=len_seq-2; t>=0; t--){
        vpath[t] = path[vpath[t+1]][t+1];
    }

    print_save = 0;
    codon_start=0;
    start_t=-1;

    int glen = strlen(O);
    char codon[4], utr[65];
    for (t=0; t<len_seq; t++){

        if (codon_start==0 && start_t < 0 &&
           ((vpath[t]>=M1_STATE && vpath[t]<=M6_STATE) ||
           (vpath[t]>=M1_STATE_1 && vpath[t]<=M6_STATE_1) ||
           vpath[t] == S_STATE || vpath[t] == S_STATE_1 ))
        {
            start_t=t+1;
            //printf("Note assign start_t %d (t=%d vpath %d)\n", start_t, t, vpath[t]);
        }

        if (codon_start==0 &&
           (vpath[t]==M1_STATE || vpath[t]==M4_STATE ||
           vpath[t]==M1_STATE_1 || vpath[t]==M4_STATE_1)) {

            memset(dna,0,300000);
            memset(dna1,0,300000);
            memset(dna_f,0,300000);
            memset(dna_f1,0,300000);
            memset(protein,0, 100000);
            memset(insert,0,100);
            memset(delete,0,100);

            insert_id = 0;
            delete_id = 0;
            dna_id = 0;
            dna_f_id = 0;
            dna[dna_id] = O[t];
            dna_start_t = t + 1; //Ye April 21, 2016
            dna_f[dna_f_id] = O[t];
            //printf("Note start dna: t = %d, dna_id %d, dna_f_id %d, add %c\n", t, dna_id, dna_f_id, O[t]);
            start_orf= t + 1;
            prev_match = vpath[t];

            if (vpath[t] < M6_STATE){
                codon_start=1;
            }else{
                codon_start=-1;
            }

        } else if (codon_start!=0 && (vpath[t]==E_STATE || vpath[t]==E_STATE_1 || t==len_seq-1)){

            if (vpath[t]==E_STATE || vpath[t]==E_STATE_1){
                end_t = t + 3;
            } else {
                end_t = t + 1;

                /* FGS1.12 start: remove incomplete codon */
                temp_t = t;
                while(vpath[temp_t] != M1_STATE && vpath[temp_t] != M4_STATE  &&
                        vpath[temp_t] != M1_STATE_1  && vpath[temp_t] != M4_STATE_1){

                    dna_f[dna_f_id] = '\0';
                    dna_f_id--;

                    dna[dna_id] = '\0';
                    dna_id--;

                    temp_t--;
                }
            /* FGS1.12 end: remove incomplete codon */
            }
            final_score = (alpha[vpath[end_t-4]][end_t-4]- alpha[vpath[start_t+2]][start_t+2] ) / (end_t-start_t-5);
            frame = start_orf % 3;
            if (frame == 0){
                frame = 3;
            }

            if (dna_id > gene_len  ){
                if (codon_start == 1){
                    if(start_t == dna_start_t - 3) { //add complete start codon to dna, Ye April 21, 2016
                        strcpy(dna_tmp, dna);
                        sprintf(dna, "%c%c%c%s", O[start_t - 1], O[start_t], O[start_t + 1], dna_tmp);
                        //printf("add start codon to dna: %c%c%c\n", O[start_t-1], O[start_t], O[start_t+1]);
                        //printf("old dna %d %s\n", strlen(dna_tmp), dna_tmp);
                        //printf("new dna %d %s\n", strlen(dna), dna);
                    }
                    if( refine ) { //add refinement of the start codons here, Ye, April 16, 2016
                        int start_old = start_t;
                        codon[0] = 0;
                        strncpy(codon, O + start_old - 1, 3);
                        codon[3] = 0;
                        int s = 0;
                        //find the optimal start codon within 30bp up- and downstream of start codon
                        double e_save;
                        int s_save;
                        while( ( ! (!strcmp(codon, "TAA") || !strcmp(codon, "TAG") || !strcmp(codon, "TGA")) ) && (start_old-1-s-35>=0)) {
                            if(!strcmp(codon, "ATG") || !strcmp(codon, "GTG") || !strcmp(codon, "TTG")) {
                                utr[0] = 0;
                                strncpy(utr, O+start_old-1-s-30,63);
                                utr[63] = 0;
                                //printf("check s=%d, codon %s\n", s, codon);
                                double freq_sum = 0;
                                for(j = 0; j < strlen(utr) - 2; j ++) {
                                    int idx = trinucleotide(utr[j], utr[j+1], utr[j+2]);
                                    freq_sum -= train_ptr->start[cg][j][idx];
                                    //printf("j=%d, key=%c%c%c %d, start %lf\n", j, utr[j], utr[j+1], utr[j+2], idx, train_ptr->start[cg][j][idx]);
                                }
                                if(s == 0) {
                                    e_save = freq_sum;
                                    s_save = s;
                                }
                                else if(freq_sum < e_save) {
                                    e_save = freq_sum;
                                    s_save = -1 * s;
                                }
                                //printf("s=%d freq_sum %lf\n", s, freq_sum);
                                //getchar();
                            }
                            s += 3;
                            codon[0] = 0;
                            strncpy(codon, O+start_old-1-s, 3);
                            codon[3] = 0;
                        }
                        start_t = start_old+s_save;
                        //update dna
                        if(s_save != 0) {
                            //printf("start refined + %d -> %d\n", start_old, start_t);
                            dna[0] = 0;
                            strncpy(dna, O + start_t - 1, end_t - start_t + 1);
                            dna[end_t - start_t + 1] = 0;
                        }
                    }
                    fprintf(fp_out, "%d\t%d\t+\t%d\t%lf\t", start_t, end_t, frame, final_score);
                    fprintf(fp_out, "I:");
                    for (i=0; i<insert_id; i++){
                        fprintf(fp_out, "%d,", insert[i]);
                    }
                    fprintf(fp_out, "\tD:");
                    for (i=0; i<delete_id; i++){
                        fprintf(fp_out, "%d,", delete[i]);
                    }
                    fprintf(fp_out, "\n");

                    fprintf(fp_aa, "%s_%d_%d_+\n", head_short, start_t, end_t);
                    fprintf(fp_dna, "%s_%d_%d_+\n", head_short, start_t, end_t);

                    //printf("dna-start %c%c%c exp %c%c%c (%c%c%c)\n", dna[0], dna[1], dna[2], O[start_t-1], O[start_t], O[start_t+1], O[start_t+2], O[start_t+3], O[start_t+4]);
                    //printf("dna-len %d start_t %d end_t %d exp-len %d diff %d\n", strlen(dna), start_t, end_t, end_t - start_t + 1, end_t - start_t + 1 - strlen(dna));
                    get_protein(dna,protein,1, whole_genome);
                    fprintf(fp_aa, "%s\n", protein);
                    if (format==0){
                        fprintf(fp_dna, "%s\n", dna);
                    }else if (format==1){
                        fprintf(fp_dna, "%s\n", dna_f);
                    }
                } else if (codon_start==-1){
                    //printf("reverse strand dna-len %d, start_t %d, dna_start %d, add-up %d, end_t %d\n", strlen(dna), start_t, dna_start_t, dna_start_t + strlen(dna), end_t);
                    //getchar();
                    if(dna_start_t + strlen(dna) == end_t - 2) { //add complete start codon (on reverse strand) to dna, Ye April 21, 2016
                        strcpy(dna_tmp, dna);
                        sprintf(dna, "%s%c%c%c", dna_tmp, O[end_t-3], O[end_t-2], O[end_t-1]);
                        //printf("add start codon on the reverse strand to dna: %c%c%c\n", O[end_t-3], O[end_t-2], O[end_t-1]);
                    }
                    if(refine) { //add refinement of the start codons here, Ye, April 16, 2016
                        int end_old = end_t; //reverse
                        codon[0] = 0;
                        strncpy(codon, O + end_t-1-2, 3);
                        codon[3] = 0;
                        int s = 0;
                        //find the optimal start codon within 30bp up- and downstream of start codon
                        double e_save;
                        int s_save;
                        while((!(!strcmp(codon, "TTA") || !strcmp(codon, "CTA") || !strcmp(codon, "TCA"))) && (end_old-2+s+35 < glen)) {
                            if(!strcmp(codon, "CAT") || !strcmp(codon, "CAC") || !strcmp(codon, "CAA")) {
                                utr[0] = 0;
                                strncpy(utr, O+end_old-1-2+s-30,63);
                                utr[63] = 0;
                                //printf("check s=%d, codon %s\n", s, codon);
                                double freq_sum = 0;
                                for(j = 0; j < strlen(utr) - 2; j ++) {
                                    int idx = trinucleotide(utr[j], utr[j+1], utr[j+2]);
                                    freq_sum -= train_ptr->stop1[cg][j][idx]; //stop1?? Ye, April 18, 2016
                                    //printf("j=%d, key=%c%c%c %d, stop1 %lf\n", j, utr[j], utr[j+1], utr[j+2], idx, train_ptr->stop1[cg][j][idx]);
                                }
                                if(s == 0) {
                                    e_save = freq_sum;
                                    s_save = s;
                                } else if( freq_sum < e_save ) {
                                    e_save = freq_sum;
                                    s_save = s;
                                }
                                //printf("s=%d freq_sum %lf\n", s, freq_sum);
                                //getchar();
                            }
                            s += 3;
                            codon[0] = 0;
                            strncpy(codon, O+end_old-1-2+s, 3);
                            codon[3] = 0;
                        }
                        end_t = end_old+s_save;
                        //update dna
                        if(s_save != 0) {
                            //printf("start refined - end %d -> %d\n", end_old, end_t);
                            dna[0] = 0;
                            strncpy(dna, O + start_t - 1, end_t - start_t + 1);
                            dna[end_t - start_t + 1] = 0;
                        }
                    }

                    fprintf(fp_out, "%d\t%d\t-\t%d\t%lf\t", start_t, end_t, frame, final_score);
                    fprintf(fp_out, "I:");
                    for (i=0; i<insert_id; i++){
                        fprintf(fp_out, "%d,", insert[i]);
                    }
                    fprintf(fp_out, "\tD:");
                    for (i=0; i<delete_id; i++){
                        fprintf(fp_out, "%d,", delete[i]);
                    }
                    fprintf(fp_out, "\n");

                    fprintf(fp_aa, "%s_%d_%d_-\n", head_short, start_t, end_t);
                    fprintf(fp_dna, "%s_%d_%d_-\n", head_short, start_t, end_t);

                    get_protein(dna,protein,-1, whole_genome);
                    get_rc_dna(dna, dna1);
                    get_rc_dna_indel(dna_f, dna_f1);
                    fprintf(fp_aa, "%s\n", protein);
                    if (format==0){
                        fprintf(fp_dna, "%s\n", dna1);
                    }else if (format==1){
                        fprintf(fp_dna, "%s\n", dna_f1);
                    }
                }
            }
            codon_start=0;
            start_t = -1;
            end_t = -1;
            dna_id = 0;
            dna_f_id = 0;

        } else if (codon_start!=0 &&
          ((vpath[t]>=M1_STATE && vpath[t]<=M6_STATE) ||
           (vpath[t]>=M1_STATE_1 && vpath[t]<=M6_STATE_1)) &&
           vpath[t]-prev_match<6) {

            if (vpath[t] < prev_match){
                out_nt = vpath[t]+6-prev_match;
            } else {
                out_nt = vpath[t]-prev_match;
            }
            for (kk=0; kk<out_nt; kk++){   /* for deleted nt in reads */
                dna_id ++;
                dna[dna_id] = 'N';
                //printf("dna_id %d, dna-len %d\n", dna_id, strlen(dna));
                dna_f_id ++;
                dna_f[dna_f_id] = 'x';
                if (kk>0){
                    delete[delete_id]=t+1;
                    delete_id++;
                }
            }
            dna[dna_id]=O[t];
            //printf("dna_id %d, add %d %c dna-len %d\n", dna_id, t, O[t], strlen(dna));
            dna_f[dna_f_id]=O[t];
            prev_match = vpath[t];

        } else if (codon_start!=0 &&
          ((vpath[t]>=I1_STATE && vpath[t]<=I6_STATE) ||
           (vpath[t]>=I1_STATE_1 && vpath[t]<=I6_STATE_1))) {

            dna_f_id ++;
            dna_f[dna_f_id] = tolower(O[t]);
            insert[insert_id]=t+1;
            insert_id++;

        } else if (codon_start!=0 && vpath[t]==R_STATE) {
            /* for long NNNNNNNNN, pretend R state */
            codon_start=0;
            start_t=-1;
            end_t = -1;
            dna_id=0;
            dna_f_id=0;

        }
    }
#ifdef viterbi_out_flg
    printf("End of viterbi\n");
#endif
    free_ivector(vpath);
    return;
}

int get_prob_from_cg(HMM *hmm_ptr, TRAIN *train_ptr, char *O){ //change from void to int, Ye, April 18, 2016
    int cg_id = -1;
    int cg_count=0;
    int len_seq;
    int i,j,k;

    len_seq = strlen(O);
    for (i=0; i<len_seq; i++){
        if ((O[i] == 'C'||O[i] =='c') || (O[i] == 'G'||O[i] == 'g') ){
            ++cg_count;
        }
    }
    cg_count = floor((cg_count*1.0/len_seq)*100)-26;
    if (cg_count < 0){
        cg_count = 0;
    }else if (cg_count > 43){
        cg_count = 43;
    }

    memcpy(hmm_ptr->e_M, train_ptr->trans[cg_count], sizeof(hmm_ptr->e_M));
    memcpy(hmm_ptr->e_M_1, train_ptr->rtrans[cg_count], sizeof(hmm_ptr->e_M_1));
    memcpy(hmm_ptr->tr_R_R, train_ptr->noncoding[cg_count], sizeof(hmm_ptr->tr_R_R));
    memcpy(hmm_ptr->tr_S, train_ptr->start[cg_count], sizeof(hmm_ptr->tr_S));
    memcpy(hmm_ptr->tr_E, train_ptr->stop[cg_count], sizeof(hmm_ptr->tr_E));
    memcpy(hmm_ptr->tr_S_1, train_ptr->start1[cg_count], sizeof(hmm_ptr->tr_S_1));
    memcpy(hmm_ptr->tr_E_1, train_ptr->stop1[cg_count], sizeof(hmm_ptr->tr_E_1));
    memcpy(hmm_ptr->S_dist, train_ptr->S_dist[cg_count], sizeof(hmm_ptr->S_dist));
    memcpy(hmm_ptr->E_dist, train_ptr->E_dist[cg_count], sizeof(hmm_ptr->E_dist));
    memcpy(hmm_ptr->S1_dist, train_ptr->S1_dist[cg_count], sizeof(hmm_ptr->S1_dist));
    memcpy(hmm_ptr->E1_dist, train_ptr->E1_dist[cg_count], sizeof(hmm_ptr->E1_dist));

    return cg_count;
}


int get_prob_form_cg_graph(HMM *hmm_ptr, TRAIN *train_ptr, Graph *g){
#ifdef viterbi_out_flg
    printf("get_prob_from_cg(HMM*,TRAIN*, Grpah*))\n");
#endif

    unsigned long long cg_count = 0, total_len = 0;
    int i, j;
    for (i = 0; i < g->n_edge; ++i){
        total_len += g->seq_len[i];
        //printf("i = %d, seq_len = %d\n%s\n", i, g->seq_len[i], g->obs_seq[i]);
        for (j = 0; j < g->seq_len[i]; ++j){
            if (g->obs_seq[i][j] == 'C' || g->obs_seq[i][j] =='c' || g->obs_seq[i][j] == 'G' || g->obs_seq[i][j] == 'g'){
                ++cg_count;
            }
        }
    }
    cg_count = floor((double)cg_count / total_len * 100.0) - 26;
    if (cg_count < 0) {
        cg_count = 0;
    }
    else if (cg_count > 43) {
        cg_count = 43;
    }

    memcpy(hmm_ptr->e_M, train_ptr->trans[cg_count], sizeof(hmm_ptr->e_M));
    memcpy(hmm_ptr->e_M_1, train_ptr->rtrans[cg_count], sizeof(hmm_ptr->e_M_1));
    memcpy(hmm_ptr->tr_R_R, train_ptr->noncoding[cg_count], sizeof(hmm_ptr->tr_R_R));
    memcpy(hmm_ptr->tr_S, train_ptr->start[cg_count], sizeof(hmm_ptr->tr_S));
    memcpy(hmm_ptr->tr_E, train_ptr->stop[cg_count], sizeof(hmm_ptr->tr_E));
    memcpy(hmm_ptr->tr_S_1, train_ptr->start1[cg_count], sizeof(hmm_ptr->tr_S_1));
    memcpy(hmm_ptr->tr_E_1, train_ptr->stop1[cg_count], sizeof(hmm_ptr->tr_E_1));
    memcpy(hmm_ptr->S_dist, train_ptr->S_dist[cg_count], sizeof(hmm_ptr->S_dist));
    memcpy(hmm_ptr->E_dist, train_ptr->E_dist[cg_count], sizeof(hmm_ptr->E_dist));
    memcpy(hmm_ptr->S1_dist, train_ptr->S1_dist[cg_count], sizeof(hmm_ptr->S1_dist));
    memcpy(hmm_ptr->E1_dist, train_ptr->E1_dist[cg_count], sizeof(hmm_ptr->E1_dist));


    return cg_count;
}


void get_train_from_file(char *filename, HMM *hmm_ptr, char *mfilename, char *mfilename1, char *nfilename,
			 char *sfilename,char *pfilename,char *s1filename,char *p1filename,char *dfilename, TRAIN *train_ptr){

  int i, j, k, p;
  double prob;
  FILE *fp, *fpm, *fpm1, *fpn, *fps, *fpp, *fps1, *fpp1, *fpd;

  char name[10];
  char head[20];
  char start[10];
  char end[10];

  /* probabilities saved in log Ye April 18, 2016 */
  /* start <- ./train/start (start in forward)
     stop <- ./train/stop (stop in forward)
     start1 <- ./train/stop1 (start in reverse)
     stop1 <- ./train/start1 (stop in reverse)
  */

  /****************************************************/
  /* transition                                       */
  /****************************************************/
  fp = fopen (filename , "r");

  /* Transition */
  fscanf(fp, "%s", head);
  for (i=0; i<14; i++){
    fscanf(fp, "%s %lf", name, &prob);
    hmm_ptr->tr[tr2int(name)] = log(prob);
  }

  /* TransitionMI */
  fscanf(fp, "%s", head);
  for (i=0; i<16; i++){
    fscanf(fp, "%s %s %lf\n", start, end, &prob);
    hmm_ptr->tr_M_I[nt2int(start[0])][nt2int(end[0])] = log(prob);
  }

  /* TransitionII */
  fscanf(fp, "%s", head);
  for (i=0; i<16; i++){
    fscanf(fp, "%s %s %lf", start, end, &prob);
    hmm_ptr->tr_I_I[nt2int(start[0])][nt2int(end[0])] = log(prob);
 }

  /* PI */
  fscanf(fp, "%s", head);
  for (i=0; i<NUM_STATE; i++){
    fscanf(fp, "%s %lf", name, &prob);
    hmm_ptr->pi[i] = log(prob);
  }
  fclose(fp);


  /****************************************************/
  /* M state transition                               */
  /****************************************************/
  fpm = fopen (mfilename , "r");
  for (p=0; p<44; p++){                        /* cg */
    fscanf(fpm, "%s", head);
    for (i=0; i<6; i++){                       /* period */
      for (j=0; j<16; j++){                    /* condition */
	for (k=0; k<4; k++){                   /* emission */
	  fscanf(fpm, "%lf", &prob);
	  train_ptr->trans[p][i][j][k] = log(prob);
	}
      }
    }
  }
  fclose(fpm);


  /****************************************************/
  /* M state_1 transition                             */
  /****************************************************/
  fpm1 = fopen (mfilename1 , "r");
  for (p=0; p<44; p++){
    fscanf(fpm1, "%s", head);
    for (i=0; i<6; i++){
      for (j=0; j<16; j++){
	for (k=0; k<4; k++){
	  fscanf(fpm1, "%lf", &prob);
	  train_ptr->rtrans[p][i][j][k] = log(prob);
	}
      }
    }
  }
  fclose(fpm1);


  /****************************************************/
  /* noncoding state  transition                      */
  /****************************************************/
  fpn = fopen (nfilename, "r");
  for (p=0; p<44; p++){
    fscanf(fpn, "%s", head);
    for (j=0; j<4; j++){
      for (k=0; k<4; k++){
	fscanf(fpn, "%lf", &prob);
	train_ptr->noncoding[p][j][k] = log(prob);
      }
    }
  }
  fclose(fpn);


  /****************************************************/
  /* start                                            */
  /****************************************************/
  fps = fopen (sfilename, "r");
  for (p=0; p<44; p++){
    fscanf(fps, "%s", head);
    for (j=0; j<61; j++){
      for (k=0; k<64; k++){
	fscanf(fps, "%lf", &prob);
	train_ptr->start[p][j][k] = log(prob);
      }
    }
  }
  fclose(fps);


  /****************************************************/
  /* stop                                             */
  /****************************************************/
  fpp = fopen (pfilename, "r"); //sfilename->pfilename, Ye, April 18, 2016
  for (p=0; p<44; p++){
    fscanf(fpp, "%s", head);
    //for (j=0; j<58; j++){ //58->61, Ye, April 18, 2016
    for (j=0; j<61; j++){
      for (k=0; k<64; k++){
	fscanf(fpp, "%lf", &prob);
	train_ptr->stop[p][j][k] = log(prob);
      }
    }
  }
  fclose(fpp);


  /****************************************************/
  /* start1                                           */
  /****************************************************/
  fps1 = fopen (s1filename, "r");
  for (p=0; p<44; p++){
    fscanf(fps1, "%s", head);
    for (j=0; j<61; j++){ //58->61 Ye, April 18, 2016
      for (k=0; k<64; k++){
	fscanf(fps1, "%lf", &prob);
	train_ptr->start1[p][j][k] = log(prob);
      }
    }
  }
  fclose(fps1);


  /****************************************************/
  /* stop1                                            */
  /****************************************************/
  fpp1 = fopen (p1filename, "r");
  for (p=0; p<44; p++){
    fscanf(fpp1, "%s", head);
    for (j=0; j<61; j++){
      for (k=0; k<64; k++){
	fscanf(fpp1, "%lf", &prob);
	train_ptr->stop1[p][j][k] = log(prob);
      }
    }
  }
  fclose(fpp1);


  /****************************************************/
  /* pwm distribution                                 */
  /* S_dist, E_dist, S1_dist, E1_dist NOT in log      */
  /****************************************************/
  fpd = fopen (dfilename, "r");
  for (p=0; p<44; p++){
    fscanf(fpd, "%s", head);
    for (k=0; k<6; k++){
      fscanf(fpd, "%lf", &prob);
      train_ptr->S_dist[p][k] = prob;
    }
    for (k=0; k<6; k++){
      fscanf(fpd, "%lf", &prob);
      train_ptr->E_dist[p][k] = prob;
    }
    for (k=0; k<6; k++){
      fscanf(fpd, "%lf", &prob);
      train_ptr->S1_dist[p][k] = prob;
    }
    for (k=0; k<6; k++){
      fscanf(fpd, "%lf", &prob);
      train_ptr->E1_dist[p][k] = prob;
    }
  }
  fclose(fpd);

}

void free_hmm(HMM *hmm_ptr){

  free_dvector(hmm_ptr->pi);
}

void dump_memory(void *p, int size)
{
	int i, s;
	double *c;
	c = (double*)p;
	s = size / sizeof(double);

	printf("Dump size %d\n", size);
	for (i = 0; i < s; i++)
	{
		if (i > 0 && i % 10 == 0)
		{
			printf("\n");
		}
		printf("%f ", *c);
		c++;
	}
	printf("\n");
}

void free_ViterbiResult(ViterbiResult* res){
    if (res->alpha) {
        free_dmatrix(res->alpha, NUM_STATE);
    }
    if (res->path) {
        free_imatrix(res->path, NUM_STATE);
    }
    free(res->O);
}


Graph read_graph(FILE *fp, FILE *fp_matr){
    Graph res;
    size_t i, j;

    fscanf(fp_matr, "%zd %zd", &res.n_edge, &res.overlap);
    //res.adjacency_matrix = (int**)malloc(res.n_edge * sizeof(int*));
    res.adjacency_matrix = imatrix(res.n_edge, res.n_edge);
    for (i = 0; i < res.n_edge; ++i){
        //res.adjacency_matrix[i] = (int*)malloc(res.n_edge * sizeof(int));
        for (j = 0; j < res.n_edge; ++j)
            fscanf(fp_matr, "%d", &res.adjacency_matrix[i][j]);
    }

    res.seq_len = (int*)malloc(res.n_edge * sizeof (int));
    res.dead_end_flg = (int*)malloc(res.n_edge * sizeof(int));
    res.head = (char**)malloc(res.n_edge * sizeof(char*));
    res.obs_seq = (char**)malloc(res.n_edge * sizeof(char*));

    char tmp_str[STRINGLEN] = "";

    fgets(tmp_str, sizeof(tmp_str), fp);
    for (i = 0; i < res.n_edge; ++i){
        if (tmp_str[strlen(tmp_str) - 1] == 10) tmp_str[strlen(tmp_str) - 1] = '\0';
        res.head[i] = (char*)malloc((strlen(tmp_str) + 1) * sizeof(char));
        strcpy(res.head[i], tmp_str);

        res.seq_len[i] = 0;
        while (fgets(tmp_str, sizeof(tmp_str), fp)){
            if (tmp_str[0] == '>')
                break;
            res.seq_len[i] += strlen(tmp_str) - 1;
        }
        res.obs_seq[i] = (char*)malloc((res.seq_len[i] + 1) * sizeof(char));
    }
    rewind(fp);

    fgets(tmp_str, sizeof(tmp_str), fp);
    for (i = 0; i < res.n_edge; ++i){
        while(fgets(tmp_str, sizeof(tmp_str), fp)){
            if (tmp_str[0] == '>')
                break;
            for (j = 0; j < strlen(tmp_str); ++j)
                if (!isalpha(tmp_str[j]))
                    tmp_str[j] = '\0';
            strcat(res.obs_seq[i], tmp_str);
        }
    }

    for (i = 0; i < res.n_edge; ++i){
        for (j = 0; j < strlen(res.obs_seq[i]); ++j)
            res.obs_seq[i][j] = toupper(res.obs_seq[i][j]);
    }

    return res;
}

void free_graph(Graph *g){
    int i;
    if (g->seq_len)
        free(g->seq_len);
    if (g->obs_seq){
        for (i = 0; i < g->n_edge; ++i)
            if (g->obs_seq[i])
                free(g->obs_seq[i]);
        free(g->obs_seq);
    }
    if (g ->adjacency_matrix){
        for (i = 0; i < g->n_edge; ++i)
            if (g->adjacency_matrix[i])
                free(g->adjacency_matrix[i]);
        free(g->adjacency_matrix);
    }

    if (g->head){
        for (i = 0; i < g->n_edge; ++i){
            if (g->head[i])
                free(g->head[i]);
        }
        free(g->head);
    }

    if (g->edge_results){
        for (i = 0; i < g->n_edge; ++i){
            free_ViterbiResult(&g->edge_results[i]);
        }
    }

    return;
}





Edge *create_raw_edge(Graph* g){
    Edge *ans = (Edge*)malloc(sizeof(Edge));
    ans->g = g;
    return ans;
}



double any_state_prob(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int *prev_ind, int n_prev, size_t overlap, int whole_genome){
    TmpResult ans_res, temp;
    double max_dbl = 10000000000.0;
    int to = nt2int(curr_res->O[t]), group = state2group(i);

    if (n_prev == 0) {
        switch (group){
            case M_GROUP:
                ans_res = match_state_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, overlap, whole_genome, to);
            break;
            case I_GROUP:
                ans_res = insertion_state_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, overlap, to);
            break;
            case R_GROUP:
                ans_res = non_coding_state_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, overlap, to);
            break;
            case E_GROUP:
                ans_res = end_state_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, overlap);
            break;
            case S_GROUP:
                ans_res = start_state_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, overlap);
                break;
            case M_GROUP_1:
                ans_res = match_state1_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, overlap, whole_genome, to);
            break;
            case I_GROUP_1:
                ans_res = insertion_state1_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, overlap, to);
            break;
            case E_GROUP_1:
                ans_res = end_state1_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, overlap);
            break;
            case S_GROUP_1:
                ans_res = start_state1_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, overlap);
            break;
            default:
                return max_dbl + 1;
            break;
        }
    } else {
        if (t == 0) {
            int j, curr_prev_ind;
            for (j = 0; j < n_prev; ++j){
                curr_prev_ind = prev_ind[j];
                switch (group){
                    case M_GROUP:
                        temp = match_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap, whole_genome, to);
                    break;
                    case I_GROUP:
                        temp = insertion_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap, to);
                    break;
                    case R_GROUP:
                        temp = non_coding_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap, to);
                    break;
                    case E_GROUP:
                        temp = end_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap);

                    break;
                    case S_GROUP:
                        temp = start_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap);
                    break;
                    case M_GROUP_1:
                        temp = match_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap, whole_genome, to);
                    break;
                    case I_GROUP_1:
                        temp = insertion_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap, to);
                    break;
                    case E_GROUP_1:
                        temp = end_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap);
                    break;
                    case S_GROUP_1:
                        temp = start_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap);
                    break;
                    default:
                        return max_dbl + 1;
                }

                if (group == E_GROUP || group == S_GROUP || group == E_GROUP_1 || group == S_GROUP_1){
                    if (j == 0 || temp.alpha2 < ans_res.alpha2) {
                        ans_res = temp;
                    }
                } else {
                    if (j == 0 || temp.alpha < ans_res.alpha) {
                        ans_res = temp;
                    }
                }
            }

            if (group == I_GROUP){
                if (ans_res.path >= M1_STATE && ans_res.path <= M6_STATE){
                    curr_res->temp_i[i - I1_STATE] = -(prev_res[ans_res.prev_ind].len_seq - 1);
                }
            } else if (group == I_GROUP_1) {
                if (ans_res.path >= M1_STATE_1 && ans_res.path <= M6_STATE_1) {
                    curr_res->temp_i_1[i - I1_STATE_1] = -(prev_res[ans_res.prev_ind].len_seq - 1);
                }
            }
        } else {
            int curr_prev_ind = curr_res->curr_column_prev[i];
            switch (group){
                case M_GROUP:
                    ans_res = match_state_prob_eval(hmm_ptr, t, i, curr_res,prev_res, curr_prev_ind, n_prev, overlap, whole_genome, to);
                break;
                case M_GROUP_1:
                    ans_res = match_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap, whole_genome, to);
                break;
                case I_GROUP:
                    ans_res = insertion_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap, to);
                break;
                case I_GROUP_1:
                    ans_res = insertion_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap, to);
                break;
                case R_GROUP:
                    ans_res = non_coding_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap, to);
                break;
                case E_GROUP:
                    ans_res = end_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap);
                break;
                case S_GROUP:
                    ans_res = start_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap);
                break;
                case E_GROUP_1:
                    ans_res = end_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap);
                break;
                case S_GROUP_1:
                    ans_res = start_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, overlap);
                break;
                default:
                    return max_dbl + 1;

            }
        }
    }
    if ((i == E_STATE || i == E_STATE_1 || i == S_STATE_1 || i == S_STATE) && (t < curr_res->len_seq - 2) && curr_res->alpha[i][t] == 0 && ans_res.alpha2 != 0){
        curr_res->alpha[i][t + 2] = ans_res.alpha2;

        curr_res->alpha[i][t + 1] = max_dbl;
        curr_res->path[i][t + 1] = i;
        curr_res->path[i][t + 2] = i;

        if (i == E_STATE) {
            curr_res->alpha[M6_STATE][t + 2] = max_dbl;
            curr_res->alpha[M5_STATE][t + 1] = max_dbl;
            curr_res->alpha[M4_STATE][t] = max_dbl;
            curr_res->alpha[M3_STATE][t + 2] = max_dbl;
            curr_res->alpha[M2_STATE][t + 1] = max_dbl;
            curr_res->alpha[M1_STATE][t] = max_dbl;
        }

    }

    curr_res->tmp_curr_column_prev[i + 1] = ans_res.prev_ind;
    //curr_res->curr_column_prev[i + 1] = ans_res.prev_ind;
    curr_res->path[i][t] = ans_res.path;
    if (t == 0) {
        curr_res->first_column_prev[i + 1] = ans_res.prev_ind;
    }
    //curr_res->alpha[i][t] = ans_res.alpha;
    return ans_res.alpha;
}

int  count_from2(int t, char *O, int seq_len, char *prev_O, int prev_seq_len){
    char prev1, prev2;
    prev1 = (t >= 1) ? O[t - 1] : (prev_O ? prev_O[prev_seq_len - 1] : -1);
    prev2 = (t >= 2) ? O[t - 2] : (prev_O ? prev_O[prev_seq_len + (t - 2)] : -1);
    return (prev2 == -1 ? 2 : nt2int(prev2)) * 4 + nt2int(prev1);
}

TmpResult match_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_ind, int n_prev, size_t overlap, int whole_genome, int to) {

    double max_dbl = 10000000000.0;
    TmpResult ans_res, temp;
    int j, num_d, prev_seq_len;// = prev_res ? prev_res[prev_ind].len_seq : -1;
    double **alpha = curr_res->alpha, **prev_alpha, alpha1;// = prev_res ? prev_res[prev_ind].alpha : NULL, alpha1;
    char* O = curr_res->O, *prev_O;// = prev_res ? prev_res[prev_ind].O : NULL;
    int *temp_i = curr_res->temp_i, *prev_temp_i; //= prev_res ? prev_res[prev_ind].temp_i : NULL;


    if (t == 0 && n_prev == 0){
        ans_res.alpha = -hmm_ptr->pi[i];
        ans_res.alpha2 = 0;
        ans_res.path = 0;
        ans_res.prev_ind = -1;
    } else if (curr_res->alpha[i][t] < max_dbl ){ //t > 0 || prev_res
        int from2;
        char prev1, prev2;

        if (i == M1_STATE) {
            //from M
            j = M6_STATE;
            if (t == 0) {
                prev_alpha = prev_res[prev_ind].alpha;
                prev_seq_len = prev_res[prev_ind].len_seq - overlap;
                prev_O = prev_res[prev_ind].O;
            } else {
                prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq - overlap : -1;
                prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
            }
            from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);
            alpha1 = (t == 0) ? prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];

            ans_res.alpha = alpha1 - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[0][from2][to];
            ans_res.path = j;
            ans_res.prev_ind = (t == 0) ? prev_ind : curr_res->curr_column_prev[j + 1];

            //from D
            if (whole_genome == 0) {
                for (j = M5_STATE; j> M1_STATE; --j){
                    num_d = i - j + 6;

                    if (t == 0) {
                        prev_alpha = prev_res[prev_ind].alpha;
                        prev_seq_len = prev_res[prev_ind].len_seq - overlap;
                        prev_O = prev_res[prev_ind].O;
                    } else {
                        prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                        prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq - overlap : -1;
                        prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
                    }
                    from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);

                    alpha1 = (t == 0) ? prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
                    temp.alpha = alpha1 - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[0][from2][to] -
                            log(0.25) * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                    if (temp.alpha < ans_res.alpha) {
                        ans_res.alpha = temp.alpha;
                        ans_res.path = j;
                        ans_res.prev_ind = (t == 0) ? prev_ind : curr_res->curr_column_prev[j + 1];
                    }
                }
            }

            //from S
            j = S_STATE;
            int x = curr_res->curr_column_prev[j + 1];
            if (t == 0) {
                prev_alpha = prev_res[prev_ind].alpha;
                prev_seq_len = prev_res[prev_ind].len_seq - overlap;
                prev_O = prev_res[prev_ind].O;
            } else {
                prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq - overlap : -1;
                prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
            }
            from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);

            alpha1 = (t == 0) ? prev_alpha[S_STATE][prev_seq_len - 1] : alpha[S_STATE][t - 1];
            temp.alpha = alpha1 - hmm_ptr->e_M[0][from2][to];

            if (temp.alpha < ans_res.alpha) {
                ans_res.alpha = temp.alpha;
                ans_res.path = S_STATE;
                ans_res.prev_ind = (t == 0) ? prev_ind : curr_res->curr_column_prev[j + 1];
            }
        } else {//i = M2, ..., M6
            //from M
            j = i - 1;

            if (t == 0) {
                prev_alpha = prev_res[prev_ind].alpha;
                prev_seq_len = prev_res[prev_ind].len_seq - overlap;
                prev_O = prev_res[prev_ind].O;
            } else {
                prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq - overlap : -1;
                prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
            }
            from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);

            alpha1 = (t == 0) ? prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
            ans_res.alpha = alpha1 - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[i - M1_STATE][from2][to];
            ans_res.path = j;
            ans_res.prev_ind = (t == 0) ? prev_ind : curr_res->curr_column_prev[j + 1];
#ifdef crash67
            if (t == 1 && (ans_res.prev_ind == 17 || ans_res.prev_ind == 19)){
                 printf("sequence 71, t = 1, M_STATE\n");
            }
#endif

            //from D
            if (whole_genome == 0) {
                for (j = M6_STATE; j >= M1_STATE; --j) {
                    if (j >= i) {
                        num_d = i - j + 6;
                    } else if (j + 1 < i) {
                        num_d = i - j;
                    } else {
                        num_d = -10;
                    }
                    if (num_d > 0) {
                        if (t == 0) {
                            prev_alpha = prev_res[prev_ind].alpha;
                            prev_seq_len = prev_res[prev_ind].len_seq - overlap;
                            prev_O = prev_res[prev_ind].O;
                        } else {
                            prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                            prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq - overlap : -1;
                            prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
                        }
                        from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);


                        alpha1 = (t == 0) ? prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
                        temp.alpha = alpha1 - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[i - M1_STATE][from2][to] -
                                log(0.25) * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                        if (temp.alpha < ans_res.alpha) {
                            ans_res.alpha = temp.alpha;
                            ans_res.path = j;
                            ans_res.prev_ind = (t == 0) ? prev_ind : curr_res->curr_column_prev[j + 1];
                        }
                    }
                }
            }
        }
        //from I
        j = (i == M1_STATE) ? I6_STATE : I1_STATE + (i - M1_STATE - 1);

        if (t == 0) {
            prev_alpha = prev_res[prev_ind].alpha;
            prev_seq_len = prev_res[prev_ind].len_seq - overlap;
            prev_O = prev_res[prev_ind].O;
        } else {
            prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
            prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq - overlap : -1;
            prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
        }
        from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);
        prev_temp_i = (t == 0 && curr_res->curr_column_prev[j + 1] >= 0) ? prev_res[curr_res->curr_column_prev[j + 1]].temp_i : NULL;

        if (t < 2 && n_prev == 0){

        } else { //t >= 2 || prev_res
            prev1 = (curr_res->curr_column_prev[j + 1] < 0) ? O[temp_i[j - I1_STATE]] : // если нет предыдущего
                    ( (t == 0) ? prev_O[prev_temp_i[j - I1_STATE]] : // если есть предыдущее ребро и мы на 1 позиции текущего ребра
                        (temp_i[j - I1_STATE] >= 0 ? O[temp_i[j - I1_STATE]] : prev_O[-temp_i[j - I1_STATE]]) ); // есть предыдущее ребро; отлельно смотрим, на текущем или предыдущем ребре был последний match

            if ( (i == M2_STATE || i == M5_STATE) && prev1 == 'T' && (t + 1 < curr_res->len_seq) &&
                 ( (O[t] == 'A' && O[t + 1] == 'A') ||
                   (O[t] == 'A' && O[t + 1] == 'G') ||
                   (O[t] == 'G' && O[t + 1] == 'A'))) {
            } else {
                prev2 = (curr_res->curr_column_prev[j + 1] < 0) ? O[temp_i[j - I1_STATE] - 1] : // если нет предыдущего
                        ( (t == 0) ? prev_O[prev_temp_i[j - I1_STATE] - 1] : // если есть предыдущее ребро и мы в начале текущего
                            ( temp_i[j - I1_STATE] > 0 ? O[temp_i[j - I1_STATE] - 1] :
                                (temp_i[j - I1_STATE] < 0 ? prev_O[-temp_i[j - I1_STATE] - 1] : //проверяем, где был последний match
                                    prev_O[prev_seq_len - 1])));
                if ((i == M3_STATE || i == M6_STATE) && prev2 == 'T' &&
                        ((prev1 == 'A' && O[t] == 'A') ||
                         (prev1 == 'A' && O[t] == 'G') ||
                         (prev1 == 'G' && O[t] == 'A'))){

                } else {
                    alpha1 = (t == 0) ?  prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
                    temp.alpha = alpha1 - hmm_ptr->tr[TR_IM] - log(0.25);
                    if (temp.alpha < ans_res.alpha) {
                        ans_res.alpha = temp.alpha;
                        ans_res.path = j;
                        ans_res.prev_ind = (t == 0) ? prev_ind : curr_res->curr_column_prev[j + 1];
                    }
                }

            }
        }
    } else {
        ans_res.alpha = curr_res->alpha[i][t];
        ans_res.path = i;
        ans_res.prev_ind = curr_res->curr_column_prev[i + 1];//???;
        ans_res.alpha2 = (t + 2 < curr_res->len_seq) ? curr_res->alpha[i][t + 2] : 0;
    }
    return ans_res;
}


TmpResult match_state1_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap, int whole_genome, int to) {
    double max_dbl = 10000000000.0, temp_alpha;
    TmpResult ans_res;
    int j, num_d;
    double **alpha = curr_res->alpha, alpha1, **prev_alpha;
    char* O = curr_res->O, *prev_O;
    int *temp_i_1 = curr_res->temp_i_1, *prev_temp_i_1;
    int prev_seq_len;

    if (t == 0 && n_prev == 0) {
        ans_res.alpha = -hmm_ptr->pi[i];
        ans_res.alpha2 = 0;
        ans_res.path = 0;
        ans_res.prev_ind = -1;
    } else {
        int from2;
        char prev1, prev2, prev3;

        //from S
        j = S_STATE_1;
        if (t == 0) {
            prev_alpha = prev_res[prev_index].alpha;
            prev_seq_len = prev_res[prev_index].len_seq - overlap;
            prev_O = prev_res[prev_index].O;
        } else {
            prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
            prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq - overlap : -1;
            prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
        }
        from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);
        prev1 = (t >= 1) ? O[t - 1] : (curr_res->curr_column_prev[j + 1] >= 0 ? prev_O[prev_seq_len - 1] : -1);
        prev2 = (t >= 2) ? O[t - 2] : (curr_res->curr_column_prev[j + 1] >= 0 ? prev_O[prev_seq_len + (t - 2)] : -1);
        prev3 = (t >= 3) ? O[t - 3] : (curr_res->curr_column_prev[j + 1] >= 0 ? prev_O[prev_seq_len + (t - 3)] : -1);

        if ((i == M1_STATE_1 || i == M4_STATE_1) && (t >= 3 || prev_res) &&
                    ( (prev3 == 'T' && prev2 == 'T' && prev1 == 'A') ||
                      (prev3 == 'C' && prev2 == 'T' && prev1 == 'A') ||
                      (prev3 == 'T' && prev2 == 'C' && prev1 == 'A') )) {
            //from S
            alpha1 = (t > 0) ? alpha[S_STATE_1][t - 1] : prev_alpha[S_STATE_1][prev_seq_len - 1];
            ans_res.alpha = alpha1 - hmm_ptr->e_M_1[i - M1_STATE_1][from2][to];
            ans_res.path = S_STATE_1;
            ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[S_STATE_1 + 1];
        } else {
            if (i == M1_STATE_1) {
                //from M
                j = M6_STATE_1;
                if (t == 0) {
                    prev_alpha = prev_res[prev_index].alpha;
                    prev_seq_len = prev_res[prev_index].len_seq - overlap;
                    prev_O = prev_res[prev_index].O;
                } else {
                    prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                    prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq - overlap : -1;
                    prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
                }

                alpha1 = (t == 0) ? prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
                ans_res.alpha = alpha1 - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M_1[0][from2][to];
                ans_res.path = j;
                ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[j + 1];
                //from D
                if (!whole_genome){
                    for (j = M5_STATE_1; j >= M1_STATE_1; --j){
                        num_d = i - j + 6;
                        if (t == 0) {
                            prev_alpha = prev_res[prev_index].alpha;
                            prev_seq_len = prev_res[prev_index].len_seq - overlap;
                            prev_O = prev_res[prev_index].O;
                        } else {
                            prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                            prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq - overlap : -1;
                            prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
                        }
                        from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);

                        alpha1 = (t == 0) ? prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
                        temp_alpha = alpha1 - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M_1[0][from2][to] -
                                        log(0.25) * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                        if (temp_alpha < ans_res.alpha) {
                            ans_res.alpha = temp_alpha;
                            ans_res.path = j;
                            ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[j + 1];
                        }
                    }
                }
            } else {
                //from M
                j = i - 1;
                if (t == 0) {
                    prev_alpha = prev_res[prev_index].alpha;
                    prev_seq_len = prev_res[prev_index].len_seq - overlap;
                    prev_O = prev_res[prev_index].O;
                } else {
                    prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                    prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq - overlap: -1;
                    prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
                }
                from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);

                alpha1 = (t == 0) ?  prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
                ans_res.alpha = alpha1 - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M_1[i - M1_STATE_1][from2][to];
                ans_res.path = j;
                ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[j + 1];

                //from D
                if (!whole_genome) {
                    for (j = M6_STATE_1; j >= M1_STATE_1; --j) {
                        if (t == 0) {
                            prev_alpha = prev_res[prev_index].alpha;
                            prev_seq_len = prev_res[prev_index].len_seq - overlap;
                            prev_O = prev_res[prev_index].O;
                        } else {
                            prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                            prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq - overlap : -1;
                            prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
                        }
                        from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);

                        if (j >= i) {
                            num_d = i - j + 6;
                        } else if (j + 1 < i) {
                            num_d = i - j;
                        } else {
                            num_d = -10;
                        }

                        if (num_d > 0) {
                            alpha1 = (t == 0) ?  prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
                            temp_alpha = alpha1 - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M_1[i - M1_STATE_1][from2][to] -
                                            log(0.25) * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];

                            if (temp_alpha < ans_res.alpha) {
                                ans_res.alpha = temp_alpha;
                                ans_res.path = j;
                                ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[j + 1];
                            }
                        }
                    }
                }
            }
            //from I
            j = (i == M1_STATE_1) ? I6_STATE_1 : (I1_STATE_1 + (i - M1_STATE_1 - 1));
            if (t == 0) {
                prev_alpha = prev_res[prev_index].alpha;
                prev_seq_len = prev_res[prev_index].len_seq - overlap;
                prev_O = prev_res[prev_index].O;
            } else {
                int x = curr_res->curr_column_prev[j + 1];
                prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq - overlap : -1;
                prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
            }
            prev_temp_i_1 = (/*t == 0 &&*/ curr_res->curr_column_prev[j + 1] >= 0) ? prev_res[curr_res->curr_column_prev[j + 1]].temp_i_1 : NULL;
            from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);
            if (t < 2 && n_prev == 0) {

            } else {
                prev1 = (curr_res->curr_column_prev[j + 1] < 0) ? O[temp_i_1[j - I1_STATE_1]] :
                        ( (t == 0) ? prev_O[prev_temp_i_1[j - I1_STATE_1]] :
                            (temp_i_1[j - I1_STATE_1] >= 0 ? O[temp_i_1[j - I1_STATE_1]] : prev_O[-temp_i_1[j - I1_STATE_1]]));
                if ((i == M2_STATE_1 || i == M5_STATE_1) && (t + 1 < curr_res->len_seq) && O[t + 1] == 'A' &&
                       ( ( prev1 == 'T' && O[t] == 'T') ||
                         ( prev1 == 'C' && O[t] == 'T') ||
                         ( prev1 == 'T' && O[t] == 'C'))) {

                } else {

                    //now prev_temp_i_1 is null;
                    prev2 =(curr_res->curr_column_prev[j + 1] < 0) ? O[temp_i_1[j - I1_STATE_1] - 1] :
                            ((t == 0) ? prev_O[prev_temp_i_1[j - I1_STATE_1] - 1] :
                                (temp_i_1[j - I1_STATE_1] > 0 ? O[temp_i_1[j - I1_STATE_1] - 1] :
                                    (temp_i_1[j - I1_STATE_1] < 0 ? prev_O[prev_temp_i_1[j - I1_STATE_1] - 1] :
                                        prev_O[prev_seq_len - 1]) ));
                    if ((i == M3_STATE_1 || i == M6_STATE_1) && (temp_i_1[j - I1_STATE_1] > 0 || prev_res) && O[t] == 'A' &&
                        ( (prev2 == 'T' && prev1 == 'T') ||
                          (prev2 == 'C' && prev1 == 'T') ||
                          (prev2 == 'T' && prev1 == 'C') )){
                        //nothing
                    } else {
                        alpha1 = (t == 0) ?  prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
                        temp_alpha = alpha1 - hmm_ptr->tr[TR_IM] - log(0.25);
                        if (temp_alpha < ans_res.alpha) {
                            ans_res.alpha = temp_alpha;
                            ans_res.path = j;
                            ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[j + 1];
                        }
                    }
                }
            }
        }
    }
    return ans_res;
}

TmpResult insertion_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap, int to){
    double max_dbl = 10000000000.0;
    TmpResult ans_res, temp;
    int j, from;
    if (n_prev == 0 && t == 0){
        ans_res.alpha = -hmm_ptr->pi[i];
        ans_res.path = 0;
        ans_res.prev_ind = -1;
    } else if (t == 0 && n_prev > 0) {
        ans_res.prev_ind = prev_index;
        int prev_seq_len = prev_res[prev_index].len_seq - overlap;
        //from I
        j = i;
        from = nt2int(prev_res[prev_index].O[prev_seq_len - 1]);
        ans_res.alpha = prev_res[prev_index].alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
        ans_res.path = j;
        ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[j + 1];

        //from M
        j = i - I1_STATE + M1_STATE;
        if (i == I6_STATE) {
            temp.alpha = prev_res[prev_index].alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
        } else {
            temp.alpha = prev_res[prev_index].alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
        }
        if (temp.alpha < ans_res.alpha) {
            ans_res.alpha = temp.alpha;
            ans_res.path = j;
            ans_res.prev_ind = prev_index;
            //temp_i reassignment is in insertion_state_prob_function, because we take argmax over all previous edges
        }
    } else if (t > 0) {
        //from I
        j = i;
        from = nt2int(curr_res->O[t - 1]);
#ifdef I_inf_error
        double prev_alpha = curr_res->alpha[j][t - 1], tr_ii = hmm_ptr->tr[TR_II], tr_i_i = hmm_ptr->tr_I_I[from][to];
        int tyr = 0;
        ++tyr;
#endif
        ans_res.alpha = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
        ans_res.path = j;
        ans_res.prev_ind = curr_res->curr_column_prev[j];
        //from M
        j = i - I1_STATE + M1_STATE;
        if (i == I6_STATE) {
            temp.alpha = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
        } else {
            temp.alpha = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
        }
        if (temp.alpha < ans_res.alpha) {
            ans_res.alpha = temp.alpha;
            ans_res.path = j;
            ans_res.prev_ind = curr_res->curr_column_prev[j];
            curr_res->temp_i[i - I1_STATE] = t - 1;
        }
    } else {
        ans_res.alpha = max_dbl + 1;
    }
    return ans_res;
}

TmpResult insertion_state1_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap, int to){
    double max_dbl = 10000000000.0, alpha1, temp_alpha;
    double **alpha = curr_res->alpha, **prev_alpha ;
    TmpResult ans_res;
    char *prev_O;
    int j, from, prev_seq_len;
    int **path = curr_res->path, **prev_path;
    ans_res.alpha2 = 0;
    if (n_prev == 0 && t == 0){
        ans_res.alpha = -hmm_ptr->pi[i];
        ans_res.path = 0;
        ans_res.alpha2 = 0;
        ans_res.prev_ind = -1;
    } else {
        //from I
        j = i;
        if (t == 0) {
            prev_O = prev_res[prev_index].O;
            prev_seq_len = prev_res[prev_index].len_seq - overlap;
            from = nt2int(prev_O[prev_seq_len - 1]);
            alpha1 = prev_res[prev_index].alpha[j][prev_seq_len - 1];
        } else {
            from = nt2int(curr_res->O[t - 1]);
            alpha1 = curr_res->alpha[j][t - 1];
        }
        ans_res.alpha = alpha1 - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
        ans_res.path = j;
        ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[j + 1];

        //from M
        j = i - I1_STATE_1 + M1_STATE_1;
        int prev3_flg, prev4_flg, prev5_flg;
        if (t == 0) {
            alpha1 = prev_res[prev_index].alpha[j][prev_seq_len - 1];
            prev_path = prev_res[prev_index].path;
            prev3_flg = (prev_path[S_STATE_1][prev_seq_len - 3] != R_STATE);
            prev4_flg = (prev_path[S_STATE_1][prev_seq_len - 4] != R_STATE);
            prev5_flg = (prev_path[S_STATE_1][prev_seq_len - 5] != R_STATE);
        } else {
            alpha1 = curr_res->alpha[j][t - 1];
            prev3_flg = (t < 3) ? 1 : (curr_res->path[S_STATE_1][t - 3] != R_STATE);
            prev4_flg = (t < 4) ? 1 : (curr_res->path[S_STATE_1][t - 4] != R_STATE);
            prev5_flg = (t < 5) ? 1 : (curr_res->path[S_STATE_1][t - 5] != R_STATE);
        }
        if (prev3_flg && prev4_flg && prev5_flg){
            temp_alpha = alpha1 - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
            if (i == I6_STATE_1){
                temp_alpha -= hmm_ptr->tr[TR_GG];
            }

            if (temp_alpha < ans_res.alpha){
                ans_res.alpha = temp_alpha;
                ans_res.path = j;
                ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[j + 1];
            }
            if (t > 0) {
                curr_res->temp_i[i - I1_STATE_1] = t - 1;
            }
        }
    }
    return ans_res;
}

TmpResult non_coding_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap, int to) {
    double max_dbl = 10000000000.0;
    TmpResult ans_res;
    double temp_alpha, **prev_alpha;
    char *prev_O;
    int j, from, prev_seq_len;
    if (n_prev == 0 && t == 0){
        ans_res.alpha = -hmm_ptr->pi[i];
        ans_res.path = 0;
        ans_res.prev_ind = -1;
    } else if (t == 0 && n_prev > 0) {
        j = R_STATE;
        prev_alpha = prev_res[prev_index].alpha;
        prev_seq_len = prev_res[prev_index].len_seq - overlap;
        prev_O = prev_res[prev_index].O;

        from = nt2int(prev_O[prev_seq_len - 1]);

        ans_res.alpha = prev_alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_RR] - hmm_ptr->tr_R_R[from][to];
        ans_res.path = j;
        ans_res.prev_ind = prev_index;
#ifdef crash67
        double tr_rr = hmm_ptr->tr[TR_RR], tr_r_r = hmm_ptr->tr_R_R[from][to], last_alpha  = prev_alpha[j][prev_seq_len- 1];
        tr_rr *= 1;
        int musor;
        ++musor;
#endif
        //from E
        j = E_STATE;

        temp_alpha = prev_res[prev_index].alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_ER];
        if (temp_alpha < ans_res.alpha) {
            ans_res.alpha = temp_alpha;
            ans_res.path = j;
        }

        //from E'
        j = E_STATE_1;

        temp_alpha = prev_res[prev_index].alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_ER];
        if (temp_alpha < ans_res.alpha) {
            ans_res.alpha = temp_alpha;
            ans_res.path = j;
        }
        ans_res.alpha -= log(0.95);
    } else if (t > 0){
        from = nt2int(curr_res->O[t - 1]);
        //from R
        j = R_STATE;
        ans_res.alpha = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_RR] - hmm_ptr->tr_R_R[from][to];
        ans_res.path = j;
        ans_res.prev_ind = curr_res->curr_column_prev[j + 1];
#ifdef non_coding_issue
        if (curr_res->len_seq > 16000 && t == 1) {
            double tr_rr = hmm_ptr->tr[TR_RR], tr_r_r = hmm_ptr->tr_R_R[from][to], last_alpha  = curr_res->alpha[j][t - 1];
            tr_rr *= 1;
            int musor;
            ++musor;
        }
#endif

        //from E
        j = E_STATE;
        temp_alpha = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_ER];

        if (temp_alpha < ans_res.alpha) {
            ans_res.alpha = temp_alpha;
            ans_res.path = j;
            ans_res.prev_ind = curr_res->curr_column_prev[j + 1];
        }

        //from E1
        j = E_STATE_1;
        temp_alpha = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_ER];
        if (temp_alpha < ans_res.alpha) {
            ans_res.alpha = temp_alpha;
            ans_res.path = j;
            ans_res.prev_ind = curr_res->curr_column_prev[j + 1];
        }

        ans_res.alpha -= log(0.95);
    } else {
        ans_res.alpha = max_dbl + 1;
        ans_res.path = -2;
        ans_res.prev_ind = -2;
    }

    return ans_res;
}

TmpResult end_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap) {
    double max_dbl = 10000000000.0;
    TmpResult ans_res;
    double temp_alpha, **alpha = curr_res->alpha;
    char *O = curr_res->O;

    if(n_prev == 0 && t == 0){
        ans_res.alpha = -hmm_ptr->pi[E_STATE];
        ans_res.path = 0;
        ans_res.prev_ind = -1;
        ans_res.alpha2 = 0;
    } else {
        double **prev_alpha;
        int prev_seq_len;
        char *prev_O;
        if (alpha[E_STATE][t] == 0){
            ans_res.alpha = max_dbl;
            //ans_res.alpha2 = max_dbl + 2;
            ans_res.path = NOSTATE;
            ans_res.alpha2 = 0;
            ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[i + 1];

            if (t + 2 < curr_res->len_seq && O[t] == 'T' &&
                ((O[t + 1] == 'A' && O[t + 2] == 'A') ||
                 (O[t + 1] == 'A' && O[t + 2] == 'G') ||
                 (O[t + 1] == 'G' && O[t + 2] == 'A'))){

                ans_res.alpha2 = 2 * max_dbl;
                //transition from M6
                if(t == 0) {
                    //since t == 0 we now the edge which we consider to be previos
                    prev_alpha = prev_res[prev_index].alpha;
                    prev_seq_len = prev_res[prev_index].len_seq - overlap;
                    temp_alpha = prev_alpha[M6_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_GE];
                } else {
                    temp_alpha = alpha[M6_STATE][t - 1] - hmm_ptr->tr[TR_GE];
                }
                if(t == 2){
                    int x = 5;
                    double palpha = alpha[M6_STATE][t - 1], tr_ge = hmm_ptr->tr[TR_GE];
                    ++x;
                }
                if (temp_alpha < ans_res.alpha2) {
                    ans_res.alpha2 = temp_alpha;
                    ans_res.path = M6_STATE;
                    ans_res.prev_ind = curr_res->curr_column_prev[M6_STATE + 1];
                }
                //transition from M3
                if (t == 0){
                    prev_alpha = prev_res[prev_index].alpha;
                    prev_seq_len = prev_res[prev_index].len_seq - overlap;
                    temp_alpha = prev_alpha[M3_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_GE];
                } else {
                    temp_alpha = alpha[M3_STATE][t - 1] - hmm_ptr->tr[TR_GE];
                }
                if(t == 2){
                    int y = 5;
                    ++y;
                }
                if (temp_alpha < ans_res.alpha2) {
                    ans_res.alpha2 = temp_alpha;
                    ans_res.path = M3_STATE;
                    ans_res.prev_ind = curr_res->curr_column_prev[M3_STATE + 1];
                }

                //assignment of M1, ... M6 must be in main function

                if (O[t + 1] == 'A' && O[t + 2] == 'A'){
                    ans_res.alpha2 -= log(0.54);
                } else if (O[t + 1] == 'A' && O[t + 2] == 'G'){
                    ans_res.alpha2 -= log(0.16);
                } else if (O[t + 1] == 'G' && O[t + 2] == 'A'){
                    ans_res.alpha2 -= log(0.3);
                }

                //adjustment based on probability distribution
                double start_freq = 0, h_kd, p_kd, r_kd;
                int lbound, i;
                if (t >= 60 || !prev_res){
                    lbound = min(60, t);
                    for (i = -lbound; i <= -3; ++i){
                        start_freq -= hmm_ptr->tr_E[i + 60][trinucleotide(O[t + i], O[t + i + 1], O[t + i + 2])];
                    }
                } else { // t < 60 && prev_res
                    prev_O = (t == 0) ? prev_res[prev_index].O : prev_res[curr_res->curr_column_prev[ans_res.path + 1]].O;
                    prev_seq_len = (t == 0) ? prev_res[prev_index].len_seq - overlap : prev_res[curr_res->curr_column_prev[ans_res.path + 1]].len_seq - overlap;
                    lbound = min(60, t + prev_seq_len);
                    char nt1, nt2, nt3;
                    for (i = -lbound; i <= -3; ++i){
                        nt1 = (t + i >= 0) ? O[t + i] : prev_O[prev_seq_len + t + i];
                        nt2 = (t + i + 1 >= 0) ? O[t + i + 1] : prev_O[prev_seq_len + t + i + 1];
                        nt3 = (t + i + 2 >= 0) ? O[t + i + 2] : prev_O[prev_seq_len + t + i + 2];
                        start_freq -= hmm_ptr->tr_E[i + 60][trinucleotide(nt1, nt2, nt3)];
                    }
                }
                start_freq *= 58.0 / (lbound - 2);
                h_kd = hmm_ptr->E_dist[2] * exp(-1*pow(start_freq-hmm_ptr->E_dist[1],2)/(2*pow(hmm_ptr->E_dist[0],2)));
                r_kd = hmm_ptr->E_dist[5] * exp(-1*pow(start_freq-hmm_ptr->E_dist[4],2)/(2*pow(hmm_ptr->E_dist[3],2)));
                p_kd = h_kd / (h_kd + r_kd);
                if (p_kd < 0.01) {
                    p_kd = 0.01;
                } else if (p_kd > 0.99){
                    p_kd = 0.99;
                }
                ans_res.alpha2 -= log(p_kd);
            }
        } else {
            ans_res.alpha = curr_res->alpha[E_STATE][t];
            ans_res.path = curr_res->path[E_STATE][t];
            ans_res.prev_ind = curr_res->curr_column_prev[E_STATE + 1];
            ans_res.alpha2 = (t < curr_res->len_seq - 2) ? curr_res->alpha[E_STATE][t + 2] : 0;
        }
    }
    return ans_res;
}

TmpResult start_state1_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap){
    double max_dbl = 10000000000.0, alpha1;
    TmpResult ans_res;
    double temp_alpha, **alpha = curr_res->alpha, **prev_alpha;// = prev_res ? prev_res[prev_index].alpha : NULL;
    char *O = curr_res->O;
    int len_seq = curr_res->len_seq, prev_seq_len;// = prev_res ? prev_res[prev_index].len_seq : -1;
    if (n_prev == 0 && t == 0){
        ans_res.alpha = -hmm_ptr->pi[i];
        ans_res.alpha2 = 0;
        ans_res.path = 0;
        ans_res.prev_ind = -1;
    } else {
        if (alpha[S_STATE_1][t] == 0) {
            ans_res.alpha = max_dbl;
            ans_res.path = NOSTATE;
            ans_res.alpha2 = 0;
            ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[i + 1];

            if (t < len_seq - 2 && O[t + 2] == 'A' &&
                    ( (O[t] == 'T' && O[t + 1] == 'T') ||
                      (O[t] == 'C' && O[t + 1] == 'T') ||
                      (O[t] == 'T' && O[t + 1] == 'C'))) {

                //from R

                alpha1 = (t > 0) ? alpha[R_STATE][t - 1] : prev_res[prev_index].alpha[R_STATE][prev_res[prev_index].len_seq - overlap - 1];
                ans_res.alpha2 = alpha1 - hmm_ptr->tr[TR_RS];
                ans_res.path = R_STATE;
                ans_res.prev_ind = curr_res->curr_column_prev[R_STATE + 1];
                //не забыть в main переназначить alpha и path для кодона

                //from E'
                alpha1 = (t > 0) ? alpha[E_STATE_1][t - 1] : prev_res[prev_index].alpha[E_STATE_1][prev_res[prev_index].len_seq - overlap - 1];
                temp_alpha = alpha1 - hmm_ptr->tr[TR_ES];
                if (temp_alpha < ans_res.alpha2) {
                    ans_res.alpha2 = temp_alpha;
                    ans_res.path = E_STATE_1;
                    ans_res.prev_ind = curr_res->curr_column_prev[E_STATE_1 + 1];
                }

                //from E
                alpha1 = (t > 0) ? alpha[E_STATE][t - 1] : prev_res[prev_index].alpha[E_STATE][prev_res[prev_index].len_seq - overlap - 1];
                temp_alpha = alpha1 - hmm_ptr->tr[TR_ES1];
                if (temp_alpha < ans_res.alpha2){
                    ans_res.alpha2 = temp_alpha;
                    ans_res.path = E_STATE;
                    ans_res.prev_ind = curr_res->curr_column_prev[E_STATE + 1];
                }

                if (O[t] == 'T' && O[t+1] == 'T'){
                    ans_res.alpha2 -= log(0.54);
                } else if (O[t] == 'C' && O[t + 1] == 'T') {
                    ans_res.alpha2 -= log(0.16);
                } else if (O[t] == 'T' && O[t + 1] == 'C'){
                    ans_res.alpha2 -= log(0.30);
                }

                //adjustment based on probability distribution
                double start_freq = 0;
                for (i = 3; i <= 60; ++i){
                    if (t + i + 2 < len_seq){
                        start_freq -= hmm_ptr->tr_S_1[i - 3][trinucleotide(O[t + i], O[t + i + 1], O[t + i + 2])];
                    }
                }
                double h_kd = hmm_ptr->S1_dist[2] * exp(-1*pow(start_freq-hmm_ptr->S1_dist[1],2)/(2*pow(hmm_ptr->S1_dist[0],2)));
                double r_kd = hmm_ptr->S1_dist[5] * exp(-1*pow(start_freq-hmm_ptr->S1_dist[4],2)/(2*pow(hmm_ptr->S1_dist[3],2)));
                double p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01){
                    p_kd=0.01;
                }else if (p_kd>0.99){
                    p_kd=0.99;
                }
                ans_res.alpha2 -= log(p_kd);
            }

        } else {
            ans_res.alpha = alpha[S_STATE_1][t];
            ans_res.alpha2 = 0;
            ans_res.path = curr_res->path[S_STATE_1][t];
            ans_res.prev_ind = curr_res->curr_column_prev[S_STATE_1];
        }
    }

    return ans_res;
}


TmpResult end_state1_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap){
    double max_dbl = 10000000000.0, alpha1;
    TmpResult ans_res;
    double temp_alpha, **alpha = curr_res->alpha, **prev_alpha;// = prev_res ? prev_res[prev_index].alpha : NULL;
    char *O = curr_res->O, *prev_O;// = prev_res ? prev_res[prev_index].O : NULL;
    int len_seq = curr_res->len_seq, prev_seq_len;// = prev_res ? prev_res[prev_index].len_seq : -1;
    if (t == 0 && n_prev == 0){
        ans_res.alpha = -hmm_ptr->pi[i];
        ans_res.alpha2 = 0;
        ans_res.path = 0;
        ans_res.prev_ind = -1;
    } else {
#ifdef E1_issue
      if (t > 16955 && t < 16965){
          double curr_alpha = alpha[E_STATE_1][t];
          char c_0 = curr_res->O[t], c_1 = curr_res->O[t + 1], c_2 = curr_res->O[t + 2], c_3 = curr_res->O[t + 3];
          int empty_var = -1;
          printf("Empty_var %d\n", empty_var);
      }
#endif
        if (alpha[E_STATE_1][t] == 0) {
            ans_res.alpha = max_dbl;
            ans_res.path = NOSTATE;
            ans_res.alpha2 = 0;
            ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[i + 1];

            if (t < len_seq - 2 && O[t] == 'C' && O[t + 1] == 'A' &&
                    ( O[t + 2] == 'T' || O[t + 2] == 'A' || O[t + 2] == 'C')) {
                //transition from frame 6
                alpha1 = (t > 0) ? alpha[M6_STATE_1][t - 1] : prev_res[prev_index].alpha[M6_STATE_1][prev_res[prev_index].len_seq - overlap - 1]; //here only if t >0 || t == 0 && prev_res => no null
                ans_res.alpha2 = alpha1 - hmm_ptr->tr[TR_GE];
                ans_res.path = M6_STATE_1;
                ans_res.prev_ind = curr_res->curr_column_prev[M6_STATE_1 + 1];

                //не забыть в main сделать переназначение для первых 2 нуклеотидов стоп-кодона

                if (O[t + 2] == 'T') {
                    ans_res.alpha2 -= log(0.83);
                } else if (O[t + 2] == 'C') {
                    ans_res.alpha2 -= log(0.10);
                } else if (O[t + 2] == 'A'){
                    ans_res.alpha2 -= log(0.07);
                }

               //adjustment based on probability distribution
                double start_freq = 0;
                int lbound, i;
                if (n_prev == 0 || t >= 30){
                    lbound = min(t, 30);
                    for (i = -lbound; i<= 30; ++i) {
                        if (t + i + 2 < len_seq) {
                            start_freq -= hmm_ptr->tr_E_1[i + 30][trinucleotide(O[t + i], O[t + i + 1], O[t + i + 2])];

                        }
                    }
                } else {
                    prev_O = (t == 0) ? prev_res[prev_index].O : prev_res[curr_res->curr_column_prev[ans_res.path + 1]].O;
                    prev_seq_len = (t == 0) ? prev_res[prev_index].len_seq - overlap : prev_res[ans_res.prev_ind].len_seq - overlap;
                    lbound = min(t + prev_seq_len, 30);
                    char nt1, nt2, nt3;
                    for (i = -lbound; i <= 30; ++i){
                        nt1 = (t + i >= 0) ? O[t + i] : prev_O[prev_seq_len + t + i];
                        nt2 = (t + i + 1 >= 0) ? O[t + i + 1] : prev_O[prev_seq_len + t + i + 1];
                        nt3 = (t + i + 2 >= 0) ? O[t + i + 2] : prev_O[prev_seq_len + t + i + 2];
                        start_freq -= hmm_ptr->tr_E_1[i + 30][trinucleotide(nt1, nt2, nt3)];
                    }
                }
                start_freq *= 61.0 / (31 + lbound);

                double h_kd = hmm_ptr->E1_dist[2] * exp(-1*pow(start_freq-hmm_ptr->E1_dist[1],2)/(2*pow(hmm_ptr->E1_dist[0],2)));
                double r_kd = hmm_ptr->E1_dist[5] * exp(-1*pow(start_freq-hmm_ptr->E1_dist[4],2)/(2*pow(hmm_ptr->E1_dist[3],2)));
                double p_kd = h_kd / (h_kd + r_kd);

                if (p_kd < 0.01) {
                    p_kd = 0.01;
                } else if (p_kd > 0.99) {
                    p_kd = 0.99;
                }
                ans_res.alpha2 -= log(p_kd);
#ifdef  E1_state_debug
                if (t == 9 || t == 107){
                    char codone[4];
                    codone[3] = '\0';
                    codone[2] = curr_res->O[t + 2];
                    codone[1] = curr_res->O[t + 1];
                    codone[0] = curr_res->O[t + 0];
                    int xxl = 12;
                    ++xxl;
                    --xxl;
              }
#endif
            }
        } else {
            ans_res.alpha = alpha[E_STATE_1][t];
            ans_res.path = curr_res->path[E_STATE_1][t];
            ans_res.alpha2 = 0;
            ans_res.prev_ind = curr_res->curr_column_prev[S_STATE_1 + 1];
        }
    }

    return ans_res;
}

TmpResult start_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, size_t overlap) {
    double max_dbl = 10000000000.0;
    TmpResult ans_res;
    double temp_alpha, **alpha = curr_res->alpha;
    char *O = curr_res->O;
    if (t == 0 && n_prev == 0){
        ans_res.alpha = -hmm_ptr->pi[i];
        ans_res.path = 0;
        ans_res.prev_ind = -1;
        ans_res.alpha2 = 0;
    } else {
        double **prev_alpha;// = (prev_res) ? prev_res[prev_index].alpha : NULL;
        int prev_seq_len;// = (prev_res) ? prev_res[prev_index].len_seq : 0;
        char *prev_O;// = (prev_res) ? prev_res[prev_index].O : NULL;
        if (alpha[S_STATE][t] == 0){
            ans_res.alpha = max_dbl;
            ans_res.path = NOSTATE;
            ans_res.alpha2 = 0;
            ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[i + 1];

            if (t < curr_res->len_seq - 2 &&
                O[t + 1] == 'T' && O[t + 2] == 'G' &&
                (O[t] == 'A' || O[t] == 'A' || O[t] == 'G')) {

                //не забыть в 'main' переназначить alpha[S_STATE][t + 1] и path[S_STATE][t+1]
                //from R
                if (t == 0) {
                    ans_res.alpha2 = prev_res[prev_index].alpha[R_STATE][prev_res[prev_index].len_seq - overlap - 1] - hmm_ptr->tr[TR_RS];
                } else {
                    ans_res.alpha2 = alpha[R_STATE][t - 1] - hmm_ptr->tr[TR_RS];
                }
                ans_res.path = R_STATE;
                ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[R_STATE + 1];

                //from E state
                if (t == 0) {
                    temp_alpha = prev_res[prev_index].alpha[E_STATE][prev_res[prev_index].len_seq - overlap - 1] - hmm_ptr->tr[TR_ES];
                } else {
                    temp_alpha = alpha[E_STATE][t - 1] - hmm_ptr->tr[TR_ES];
                }
                if (temp_alpha < ans_res.alpha2){
                    ans_res.alpha2 = temp_alpha;
                    ans_res.path = E_STATE;
                    ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[E_STATE + 1];
                }

                //from E1
                if (t == 0) {
                    temp_alpha = prev_res[prev_index].alpha[E_STATE_1][prev_res[prev_index].len_seq -overlap  - 1] - hmm_ptr->tr[TR_ES1];
                } else {
                    temp_alpha = alpha[E_STATE_1][t - 1] - hmm_ptr->tr[TR_ES1];
                }
                if (temp_alpha < ans_res.alpha2) {
                    ans_res.alpha2 = temp_alpha;
                    ans_res.path = E_STATE_1;
                    ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[E_STATE_1 + 1];
                }

                if (O[t] == 'A'){
                    ans_res.alpha2 -= log(0.83);
                } else if (O[t] == 'G') {
                    ans_res.alpha2 -= log(0.10);
                } else if (O[t] == 'T'){
                    ans_res.alpha2 -= log(0.07);
                }

                //adjustment based on probability distribution
                double start_freq = 0;
                int lbound;
                if (!prev_res){
                    lbound = min(t, 30);

                    for (i = -lbound; i <= 30; ++i){
                        if (t + i + 2 < curr_res->len_seq) {
                            start_freq -= hmm_ptr->tr_S[i + 30][trinucleotide(O[t + i], O[t + i + 1], O[t + i + 2])];
                        }
                    }
                    start_freq *= 61.0 / (30 + lbound + 1);
                } else if (prev_res){
                    prev_O = (t == 0) ? prev_res[prev_index].O : prev_res[ans_res.prev_ind].O;
                    prev_seq_len = (t == 0) ? prev_res[prev_index].len_seq - overlap : prev_res[ans_res.prev_ind].len_seq - overlap;
                    lbound = min(30, t + prev_seq_len);
                    for (i = -lbound; i <= 30 ; ++i){
                        if (t + i + 2 < curr_res->len_seq) {
                            int cd1 = (t + i < 0) ? prev_O[prev_seq_len + t + i] : O[t + i];
                            int cd2 = (t + i + 1 < 0) ? prev_O[prev_seq_len + t + i + 1] : O[t + i + 1];
                            int cd3 = (t + i + 2 < 0) ? prev_O[prev_seq_len + t + i + 2] : O[t + i + 2];
                            start_freq -= hmm_ptr->tr_S[i + 30][trinucleotide(cd1, cd2, cd3)];
                        }
                    }
                    start_freq *= 61.0 / (30 + lbound + 1);
                }
                double h_kd = hmm_ptr->S_dist[2] * exp(-1*pow(start_freq-hmm_ptr->S_dist[1],2)/(2*pow(hmm_ptr->S_dist[0],2)));
                double r_kd = hmm_ptr->S_dist[5] * exp(-1*pow(start_freq-hmm_ptr->S_dist[4],2)/(2*pow(hmm_ptr->S_dist[3],2)));
                double p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01){
                    p_kd=0.01;
                }else if (p_kd>0.99){
                    p_kd=0.99;
                }
                ans_res.alpha2 -= log(p_kd);
            }
        } else {
            ans_res.alpha = alpha[S_STATE][t];
            ans_res.alpha2 = (t < curr_res->len_seq - 2) ? alpha[S_STATE][t + 2] : 0;
            ans_res.path = curr_res->path[S_STATE][t];
            ans_res.prev_ind = curr_res->curr_column_prev[S_STATE + 1];
        }
    }
    return ans_res;
}

int state2group(int i) {
    if (i >= M1_STATE && i <= M6_STATE)
        return M_GROUP;
    if (i >= I1_STATE && i <= I6_STATE)
        return I_GROUP;
    if (i == R_STATE)
        return R_GROUP;
    if (i == E_STATE)
        return E_GROUP;
    if (i == S_STATE)
        return S_GROUP;
    if (i >= M1_STATE_1 && i <= M6_STATE_1)
        return M_GROUP_1;
    if (i >= I1_STATE_1 && i <= I6_STATE_1)
        return I_GROUP_1;
    if (i == E_STATE_1)
        return E_GROUP_1;
    if (i == S_STATE_1)
        return S_GROUP_1;
    return -1;
}

void dfs(int **adj_matrix, int n_vert, int *used, int *ans, int *curr_len, int ind){
    used[ind] = 1;
    int i;
    for (i = 0; i < n_vert; ++i){
        //printf("ind = %d, i = %d, adj = %d\n", ind, i, adj_matrix[ind][i]);
        if (adj_matrix[ind][i] && !used[i]){
            //printf("ind = %d, i = %d, adj = %d\n", ind, i, adj_matrix[ind][i]);
            dfs(adj_matrix, n_vert, used, ans, curr_len, i);
        }
    }
    ans[*curr_len] = ind;
    (*curr_len)++;
    //printf("Push %d\n", ind + 1);
}

int *topological_sort(int **adj_matrix, int n_vert){
    int *used = (int*)malloc(n_vert * sizeof(int));
    int *ans = (int*)malloc(n_vert * sizeof(int));
    int i;
    for (i = 0; i < n_vert; ++i){
        used[i] = 0;
        ans[i] = 0;
    }
    int curr_len = 0;
    for (i = 0; i < n_vert; ++i){
        if (!used[i])
            dfs(adj_matrix, n_vert, used, ans, &curr_len, i);
    }
    for (i = 0; i < n_vert / 2; ++i){
        int tmp = ans[i];
        ans[i] = ans[n_vert - i - 1];
        ans[n_vert - i - 1] = tmp;
    }
    free(used);
    return ans;
}

char complementary_nucleotide(char c){
    if (c == 'A') return 'T';
    if (c == 'T') return 'A';
    if (c == 'G') return 'C';
    if (c == 'C') return 'G';
    return -1;
}

size_t get_edge_num_m(Graph *g, int ind, char strand){
   ind = (strand == '+') ? ind : -ind;
   size_t i;
   size_t ans = 0;
   for (i = 0; i < g->n_edge && !(ans); ++i){
       if (g->ind[i] == ind)
           return i;
   }
   return -1;
}

size_t get_edge_num(int* ind, size_t n_edge, int num, char strand){
    size_t i;
    for (i = 0; i < n_edge; ++i){
        if (ind[i] == num){
            if (strand == '-')
                return 2 * i + 1;
            return 2 * i;
        }
    }
    return -1;
}


Graph read_gfa(FILE *f){
    Graph ans;
    ans.n_edge = 0;
    //FILE *f = fopen(fname, "r");
    char line[500000], c;
    size_t i, j, seq_len, n_segm = 0;
    int num;
    /*
     * First pass. Count number of segments;
     */
    fgets(line, 500000, f);
    while (line[0] == 'S'){
        //sscanf(line, "%c\t%d\t%s", &c, &num, )
        ++n_segm;
        fgets(line, 500000, f);
    }
    ans.n_edge = n_segm * 2;

    ans.ind = ivector(ans.n_edge);
    ans.adjacency_matrix = imatrix(ans.n_edge, ans.n_edge);
    ans.seq_len = ivector(ans.n_edge);
    ans.dead_end_flg = ivector(ans.n_edge);
    ans.obs_seq = (char**)malloc(ans.n_edge * sizeof(char*));
    ans.head = (char**)malloc(ans.n_edge * sizeof(char*));
    rewind(f);

    //Read and save sequences
    char seq[500000], strand_from, strand_to;
    int num_from ,num_to;
    for (i = 0; i < n_segm; ++i){
        fgets(line, 500000, f);
        sscanf(line, "%c\t%d\t%s", &c, &num, seq);
        //+strand
        ans.ind[2 * i] = num;
        ans.seq_len[2 * i] = strlen(seq);
        ans.obs_seq[2 * i] = (char*)malloc(ans.seq_len[2 * i] + 1);
        strcpy(ans.obs_seq[2 * i], seq);
        ans.head[2 * i] = (char*)malloc(20);
        sprintf(ans.head[2 * i], "%d+", ans.ind[2 * i]);
        //-strand
        ans.ind[2 * i + 1] = -num;
        ans.seq_len[2 * i + 1] = ans.seq_len[2 * i];
        ans.obs_seq[2 * i  + 1] = (char*)malloc(ans.seq_len[2 * i + 1] + 1);
        ans.obs_seq[2 * i + 1][ans.seq_len[2 * i  +1]] = '\0';
        for (j = 0; j < ans.seq_len[2 * i + 1]; ++j){
            ans.obs_seq[2 * i + 1][j] = complementary_nucleotide(seq[ans.seq_len[2 * i + 1] - j - 1]);
        }
        ans.head[2 * i + 1] = (char*)malloc(20);
        sprintf(ans.head[2 * i + 1], "%d-", ans.ind[2 * i + 1]);
    }
    //read links and construct adjacency matrix
    while (fgets(line, 500000, f)){
        sscanf(line, "%c\t%d\t%c\t%d\t%c\t%zuM", &c, &num_from, &strand_from, &num_to, &strand_to, &ans.overlap);
        //straight
        size_t ind_from = get_edge_num_m(&ans, num_from, strand_from);
        size_t ind_to = get_edge_num_m(&ans, num_to, strand_to);
        //printf("%d%c; %d%c; %zu; %zu\n", num_from, strand_from, num_to, strand_to, ind_from, ind_to);
        ans.adjacency_matrix[ind_from][ind_to] = 1;
        //reverse-complementary
        //if the straight strand was '+' then
        ind_from += (strand_from == '+') ? 1 : -1;
        ind_to += (strand_to == '+') ? 1 : -1;
        //printf("%d%c; %d%c; %zu; %zu\n", num_to, (strand_to == '+') ? '-' : '+', num_from, (strand_from == '+') ? '-' : '+', ind_to, ind_from);
        ans.adjacency_matrix[ind_to][ind_from] = 1;
    }
    ans.used_backtrack = ivector(ans.n_edge);
#ifdef print_graph_after_reading
    FILE *outf = fopen("/home/andrei/Documents/Projects/Masters diploma/run_result/debug/read_summary.txt", "w");
    for (i = 0; i < ans.n_edge; ++i){
        fprintf(outf, "i = %d\t%s\t%zu\n", i, ans.head[i], strlen(ans.obs_seq[i]));
    }
    fclose(outf);
#endif
    return ans;
}


int check_dna_symbols(char* str){
    size_t i;
    for (i = 0; i < strlen(str); ++i)
        if (str[i] != 'A' && str[i] != 'T' && str[i] != 'G' && str[i] != 'C' && str[i] != '\0' && str[i] != '\n')
            return 0;
    return 1;
}

int is_fully_computed(double **viterbi_matrix, int n_state, int seq_len, int overlap){
    int i;
    for (i = 0; i < n_state; ++i){
        if (viterbi_matrix[i][seq_len - overlap] != 0)
            return 1;
    }
    return 0;
}

void tarjan_scc(int **adj_matrix, int n_vert, int v, Stack *st, int *index, int *curr_index, int *lowlink, int *onStack, Stack *rev_order){
    index[v] = *curr_index;
    lowlink[v] = *curr_index;
    ++(*curr_index);

    push(st, v);
    onStack[v] = 1;

    int i, w;
    for (i = 0; i < n_vert; ++i){
        if (adj_matrix[v][i] == 1){
            w = i;
            if (index[w] == -1) {
                tarjan_scc(adj_matrix, n_vert, w, st, index, curr_index, lowlink, onStack, rev_order);
                lowlink[v] = min(lowlink[v], lowlink[w]);
            } else if (onStack[w]){
                lowlink[v] = min(lowlink[v], index[w]);
            }
        }
    }

    if (index[v] == lowlink[v]){
        //printf("scc: ");
        while (st->curr_capacity > 0){
            w = pop(st);
            onStack[w] = 0;
            push(rev_order, w);
            //printf("%d ", w);
            if (w == v){
            //    printf("\n");
                break;
            }
        }
    }
}

int *ordering(int **adj_matrix, int n_vert){
    int *index = (int*)malloc(n_vert * sizeof(int));
    int *lowlink = (int*)malloc(n_vert * sizeof(int));
    int *onStack = (int*)malloc(n_vert * sizeof(int));
    int *ans = (int*)malloc(n_vert * sizeof(int));
    Stack reverse_order = create_stack();
    int curr_index = 0;
    int i;
    for (i = 0; i < n_vert; ++i){
        index[i] = -1;
        lowlink[i] = -1;
        onStack[i] = 0;
    }

    Stack st = create_stack();
    for (i = 0; i < n_vert; ++i){
        if (index[i] == -1){
            tarjan_scc(adj_matrix, n_vert, i, &st, index, &curr_index, lowlink, onStack, &reverse_order);
        }
    }
    for (i = 0; i < n_vert; ++i){
        ans[i] = pop(&reverse_order);
    }
    free(index);
    free(lowlink);
    free(onStack);
    return  ans;
}
