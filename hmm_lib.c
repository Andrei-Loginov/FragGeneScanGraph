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
#define I1_debug

void dump_memory(void *p, int size);

ViterbiResult viterbi(HMM *hmm_ptr, char *O, int whole_genome, ViterbiResult *prev_result, char *head) {
#ifdef viterbi_out_flg
    printf("Viterbi start\n");
#endif
    double max_dbl = 10000000000.0;

    double **alpha;                      /* viterbi prob array */
    int **path;                          /* viterbi path array */
    int i, j, t;

    double start_freq;

    int from, from0, to;   /*from0: i-2 position, from: i-1 position */
    int from2;             /* from2: i-2, i-1 for condition in emission probability */

    double temp_alpha;
    int len_seq;
    int prev_seq_len;
    int gene_len;

    int num_d;          /* the number of delete */
    int freq_id;
    double h_kd, r_kd, p_kd;

    int *temp_i = (int*)malloc(6 * sizeof(int));
    int *temp_i_1 = (int*)malloc(6 * sizeof(int));

    for (i = 0; i < 6; ++i){
        temp_i[i] = 0;
        temp_i_1[i] = 0;
    }

    int num_N=0;
    ViterbiResult res;

    char *prev_O = NULL;


    /***************************************************************/
    /* initialize                                                  */
    /***************************************************************/

    double log53 = log(0.53);
    double log16 = log(0.16);
    double log30 = log(0.30);
    double log25 = log(0.25);
    double log95 = log(0.95);
    double log54 = log(0.54);
    double log83 = log(0.83);
    double log07 = log(0.07);

    len_seq = strlen(O);
    res.alpha = (double **)dmatrix(hmm_ptr->N, len_seq);
    res.path = (int **)imatrix(hmm_ptr->N, len_seq);
    alpha = res.alpha;//29 x length(sequece)
    path = res.path;

    res.O = (char*)malloc(len_seq + 1 * sizeof(char));
    res.O = strcpy(res.O, O);
    res.len_seq = len_seq;
    res.temp_i = temp_i;;
    res.temp_i_1 = temp_i_1;

    temp_i = res.temp_i;
    temp_i_1 = res.temp_i;

    int prev_null_flg = (!prev_result) ? 1 : 0;

    if (!prev_result) {
        for (i=0; i<hmm_ptr->N; ++i){
            alpha[i][0] = -1 * hmm_ptr->pi[i];
        }

        /* stop state */
        if ((O[0] == 'T'|| O[0] == 't')  && //It would be much simpler to use toupper() function in the beginning instead of checking complex conditions.
            (((O[1] == 'A'|| O[1] == 'a') && (O[2] == 'A'|| O[2] == 'a')) ||
            ((O[1] == 'A'|| O[1] == 'a') && (O[2] == 'G'|| O[2] == 'g')) ||
            ((O[1] == 'G'|| O[1] == 'g') && (O[2] == 'A'|| O[2] == 'a')))) {

            alpha[E_STATE][0] = max_dbl;
            alpha[E_STATE][1] = max_dbl;
            path[E_STATE][1] = E_STATE;
            path[E_STATE][2] = E_STATE;

            alpha[M6_STATE][2] = max_dbl;
            alpha[M5_STATE][1] = max_dbl;
            alpha[M4_STATE][0] = max_dbl;
            alpha[M3_STATE][2] = max_dbl;
            alpha[M2_STATE][1] = max_dbl;
            alpha[M1_STATE][0] = max_dbl;


            if ((O[1] == 'A'|| O[1] == 'a') && (O[2] == 'A'|| O[2] == 'a')){
                alpha[E_STATE][2] = alpha[E_STATE][2] - log53;
            }else if ((O[1] == 'A'|| O[1] == 'a') && (O[2] == 'G'|| O[2] == 'g')){
                alpha[E_STATE][2] = alpha[E_STATE][2] - log16;
            }else if((O[1] == 'G'|| O[1] == 'g') && (O[2] == 'A'|| O[2] == 'a')){
                alpha[E_STATE][2] = alpha[E_STATE][2] - log30;
            }
        }
        //reverse-complimentary stop-codons
        if ((O[2] == 'A'|| O[0] == 'a')  &&
            (((O[0] == 'T'|| O[0] == 't') && (O[1] == 'T'|| O[1] == 't')) ||
            ((O[0] == 'C'|| O[0] == 'c') && (O[1] == 'T'|| O[1] == 't')) ||
            ((O[0] == 'T'|| O[0] == 't') && (O[1] == 'C'|| O[1] == 'c')))) {

            alpha[S_STATE_1][0] = max_dbl;
            alpha[S_STATE_1][1] = max_dbl;
            alpha[S_STATE_1][2] = alpha[S_STATE][0];
            path[S_STATE_1][1] = S_STATE_1;
            path[S_STATE_1][2] = S_STATE_1;

            alpha[M3_STATE_1][2] = max_dbl;
            alpha[M6_STATE_1][2] = max_dbl;

            if ((O[0] == 'T'|| O[0] == 't') && (O[1] == 'T'|| O[1] == 't')){
                alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - log53;
            }else if ((O[0] == 'C'|| O[0] == 'c') && (O[1] == 'T'|| O[1] == 't')){
                alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - log16;
            }else if((O[0] == 'T'|| O[0] == 't') && (O[1] == 'C'|| O[1] == 'c')){
                alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - log30;
            }
        }
    } else {
        prev_seq_len = (prev_result == NULL) ? 0 : strlen(prev_result->O);
        prev_O = prev_result->O;

        from = nt2int(prev_O[prev_seq_len - 1]);
        from0 = nt2int(prev_O[prev_seq_len - 2]);
        to = nt2int(O[0]);

        if (from==4){ from=2; } //Why do we change other nt to G? G is 2 according to defines and nt2int.
        if (from0==4){ from0=2; }
        if (to==4){
            to = 2;
            num_N += 1;
        }else{
            num_N = 0;
        }
        from2 = from0 * 4 + from;

        for (i = 0; i < 6; ++i){
            temp_i[i] = -prev_result->temp_i[i];
            temp_i_1[i] = -prev_result->temp_i_1[i];
        }

        ///
        ///M state
        ///
        for (i = M1_STATE; i <= M6_STATE; ++i){
            alpha[i][0] = match_state_prob_evaluation(0, i, hmm_ptr, &res, prev_result, whole_genome, from, from2, to);
        }

        ///
        /// I state
        ///
        for (i = I1_STATE; i <= I6_STATE; ++i){
            alpha[i][0] = insertion_state_prob_evaluation(0, i, hmm_ptr, &res, prev_result, from, from2, to);
        }

        ///
        /// M' state
        ///
        for (i = M1_STATE_1; i <= M6_STATE_1; ++i){
            alpha[i][0] = match1_state_prob_evaluation(0, i, hmm_ptr, &res, prev_result, whole_genome, from, from2, to);
        }
        ///
        /// I' state
        ///
        for (i = I1_STATE_1; i <= I6_STATE_1; ++i){
            alpha[i][0] = insertion1_state_prob_evaluation(0, i, hmm_ptr, &res, prev_result, from, from2, to);
        }

        ///
        ///Non-coding state
        ///

        alpha[R_STATE][0] = non_coding_state_prob_evaluation(0, hmm_ptr, &res, prev_result, from, to);

        ///
        ///END state
        ///
        alpha[E_STATE][0] = end_state_prob_evaluation(0, hmm_ptr,&res, prev_result);

        ///
        /// START' state
        ///
        if (alpha[S_STATE_1][0] == 0) {
            alpha[S_STATE_1][0] = max_dbl;
            path[S_STATE_1][0] = NOSTATE;

            if (len_seq > 2 && O[2] == 'A' &&
                ( (O[0] == 'T' && O[1] == 'T') ||
                  (O[0] == 'C' && O[1] == 'T') ||
                  (O[0] == 'T' && O[1] == 'C') )) {

                alpha[S_STATE_1][0] = max_dbl;
                path[S_STATE_1][0] = R_STATE;
                alpha[S_STATE_1][1] = max_dbl;
                path[S_STATE_1][1] = S_STATE_1;
                alpha[S_STATE_1][2] = prev_result->alpha[R_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_RS];
                path[S_STATE_1][2] = S_STATE_1;

                temp_alpha = prev_result->alpha[E_STATE_1][prev_seq_len - 1] - hmm_ptr->tr[TR_ES];
                if (temp_alpha < alpha[S_STATE_1][2]){
                    alpha[S_STATE_1][2] = temp_alpha;
                    path[S_STATE_1][2] = E_STATE_1;
                }

                temp_alpha = prev_result->alpha[E_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_ES1];
                if (temp_alpha < alpha[S_STATE_1][2]) {
                    alpha[S_STATE_1][2] = temp_alpha;
                    path[S_STATE_1][2] = E_STATE;
                }

                alpha[M3_STATE_1][2] = max_dbl;
                alpha[M6_STATE_1][2] = max_dbl;

                if (O[0] == 'T' && O[1] == 'T') {
                    alpha[S_STATE_1][2] -= log54;
                } else if (O[0] == 'C' && O[1] == 'T') {
                    alpha[S_STATE_1][2] -= log16;
                } else if (O[0] == 'T' && O[1] == 'C') {
                    alpha[S_STATE_1][2] -= log30;
                }


                //adjustment based on probability distribution
                start_freq = 0;
                for (i = 3; i <= 60; ++i){
                    if (i + 2 < len_seq) {
                        start_freq -= hmm_ptr->tr_S_1[i - 3][trinucleotide(O[i], O[i + 1], O[i + 2])];
                    }
                }

                h_kd = hmm_ptr->S1_dist[2] * exp(-1 * pow(start_freq - hmm_ptr->S1_dist[1], 2) / (2 * pow(hmm_ptr->S1_dist[0], 2)));
                r_kd = hmm_ptr->S1_dist[5] * exp(-1 * pow(start_freq - hmm_ptr->S1_dist[4], 2) / (2 * pow(hmm_ptr->S1_dist[3], 2)));
                p_kd = h_kd / (h_kd + r_kd);

                if (p_kd < 0.01){
                    p_kd = 0.01;
                } else if (p_kd > 0.99) {
                    p_kd = 0.99;
                }
                alpha[S_STATE_1][2] = alpha[S_STATE_1][2] - log(p_kd);
            }
        }

        ///
        ///START state
        ///
        if (alpha[S_STATE][0] == 0){
            alpha[S_STATE][0] = max_dbl;
            path[S_STATE][0] = NOSTATE;

            if (len_seq > 2 && O[1] == 'T' && O[2] == 'G' &&
                    (O[0] == 'A' || O[0] == 'G' || O[0] == 'T')){
                alpha[S_STATE][0] = max_dbl;
                alpha[S_STATE][1] = max_dbl;
                alpha[S_STATE][2] = prev_result->alpha[R_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_RS];
                path[S_STATE][0] = R_STATE;
                path[S_STATE][1] = S_STATE;
                path[S_STATE][2] = S_STATE;

                temp_alpha = prev_result->alpha[E_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_ES];
                if (temp_alpha < alpha[S_STATE][2]){
                    alpha[S_STATE][2] = temp_alpha;
                    path[S_STATE][0] = E_STATE;
                }

                temp_alpha = prev_result->alpha[E_STATE_1][prev_seq_len - 1] - hmm_ptr->tr[TR_ES1];
                if (temp_alpha < alpha[S_STATE][2]){
                    alpha[S_STATE][2] = temp_alpha;
                    path[S_STATE][0] = E_STATE_1;
                }

                if (O[0] == 'A') {
                    alpha[S_STATE][2] -= log83;
                } else if (O[0] == 'G') {
                    alpha[S_STATE][2] -= log(0.10);
                } else if (O[0] == 'T') {
                    alpha[S_STATE][2] -= log07;
                }

                //adjustment based on probability distribution

                start_freq = 0;
                int lower = min(30, prev_seq_len);
                for (i = -lower; i < 0; ++i) {
                    start_freq -= hmm_ptr->tr_S[i + 30][trinucleotide(prev_O[prev_seq_len + i], prev_O[prev_seq_len + i + 1], prev_O[prev_seq_len + i + 2])];
                }
                for (i = 0; i <= 30; ++i) {
                    if (i + 2 < len_seq) {
                        start_freq -= hmm_ptr->tr_S[i + 30][trinucleotide(O[i], O[i + 1], O[i + 2])];
                    }
                }
                start_freq *= 61.0 / (30 + lower + 1);

                h_kd = hmm_ptr->S_dist[2] * exp(-1*pow(start_freq-hmm_ptr->S_dist[1],2)/(2*pow(hmm_ptr->S_dist[0],2)));
                r_kd = hmm_ptr->S_dist[5] * exp(-1*pow(start_freq-hmm_ptr->S_dist[4],2)/(2*pow(hmm_ptr->S_dist[3],2)));
                p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01){
                    p_kd=0.01;
                }else if (p_kd>0.99){
                    p_kd=0.99;
                }
                alpha[S_STATE][2] -= log(p_kd);
            }
        }
        ///
        ///END' state
        ///

        if (alpha[E_STATE_1][0] == 0){
            alpha[E_STATE_1][0] = max_dbl;
            path[E_STATE_1][0] = NOSTATE;

            if (len_seq > 2 && O[0] == 'C' && O[1] == 'A' &&
                    (O[2] == 'A' || O[0] == 'C' || O[0] == 'T')){

                //transition from frame6
                alpha[E_STATE_1][2] = prev_result->alpha[M6_STATE_1][prev_seq_len - 1] - hmm_ptr->tr[TR_GE];
                path[E_STATE_1][0] = M6_STATE_1;
                alpha[E_STATE_1][0] = max_dbl;
                alpha[E_STATE_1][1] = max_dbl;
                path[E_STATE_1][1] = E_STATE_1;
                path[E_STATE_1][2] = E_STATE_1;

                if (O[2] == 'T') {
                    alpha[E_STATE_1][2] -= log83;
                } else if (O[2] == 'C') {
                    alpha[E_STATE_1][2] -= log(0.10);
                } else if (O[2] == 'A') {
                    alpha[E_STATE_1][2] -= log07;
                }

                //adjustment based on probability distribution

                start_freq = 0;
                int lower = min(30, prev_seq_len);
                for (i = -lower; i < 0; ++i) {
                    start_freq -= hmm_ptr->tr_E_1[i + 30][trinucleotide(prev_O[prev_seq_len + i], prev_O[prev_seq_len + i + 1], prev_O[prev_seq_len + i + 2])];
                }
                for (i = 0; i <= 30; ++i) {
                    if (i + 2 < len_seq) {
                        start_freq -= hmm_ptr->tr_E_1[i + 30][trinucleotide(O[i], O[i + 1], O[i + 2])];
                    }
                }
                start_freq *= 61.0 / (30 + lower + 1);

                h_kd = hmm_ptr->E1_dist[2] * exp(-1*pow(start_freq-hmm_ptr->E1_dist[1],2)/(2*pow(hmm_ptr->E1_dist[0],2)));
                r_kd = hmm_ptr->E1_dist[5] * exp(-1*pow(start_freq-hmm_ptr->E1_dist[4],2)/(2*pow(hmm_ptr->E1_dist[3],2)));
                p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01){
                    p_kd=0.01;
                }else if (p_kd>0.99){
                    p_kd=0.99;
                }
                alpha[E_STATE_1][2] -= log(p_kd);
            }
        }
        if (num_N > 9){
            for (i = 0; i < NUM_STATE; ++i){
                if (i != R_STATE){
                    alpha[i][0] = max_dbl;
                    path[i][0] = R_STATE;
                }
            }

        }

    }
    /******************************************************************/
    /*  fill out the rest of the columns                              */
    /******************************************************************/
    for (t = 1; t < len_seq; t++) {
        from = nt2int(O[t-1]); //from and to are emitted symbols
        if (t>1){
            from0 = nt2int(O[t-2]);
        }else{
            if (prev_result) {
                from0 = nt2int(prev_O[prev_seq_len - 1]);
            } else {
                from0 = 2;
            }
        }
        to = nt2int(O[t]);

        /* if DNA is other than ACGT, do it later */
        if (from==4){ from=2; } //Why do we change other nt to G? G is 2 according to defines and nt2int.
        if (from0==4){ from0=2; }
        if (to==4){
            to = 2;
            num_N += 1;
        }else{
            num_N = 0;
        }
        from2 = from0 * 4 + from;
#ifdef profiling
    FILE *time_f = fopen("group_times.txt", "a");
    time_t start_t, end_t;
    time(&start_t);
#endif
        /******************/
        /* M state        */
        /******************/
        for (i = M1_STATE; i <= M6_STATE; i++) {
            if (alpha[i][t] < max_dbl){
                alpha[i][t] = match_state_prob_evaluation(t, i, hmm_ptr, &res, prev_result, whole_genome, from, from2, to);
            }
        }
#ifdef profiling
    time(&end_t);
    fprintf(time_f, "t = %d, match duration = %f\n", t, difftime(end_t, start_t));
    time(&start_t);
#endif

        /******************/
        /* I state        */
        /******************/
        for (i=I1_STATE; i<=I6_STATE; i++) {
            alpha[i][t] = insertion_state_prob_evaluation(t, i, hmm_ptr, &res, prev_result, from, from2, to);
        }
#ifdef profiling
    time(&end_t);
    fprintf(time_f, "t = %d, insertion duration = %f\n", t, difftime(end_t, start_t));
    time(&start_t);
#endif
        /******************/
        /* M' state        */
        /******************/
        for (i = M1_STATE_1; i <= M6_STATE_1; ++i) {
            alpha[i][t] = match1_state_prob_evaluation(t, i, hmm_ptr, &res, prev_result, whole_genome, from, from2, to);
        }
#ifdef profiling
    time(&end_t);
    fprintf(time_f, "t = %d, match' duration = %f\n", t, difftime(end_t, start_t));
    time(&start_t);
#endif
    /******************/
    /* I' state        */
    /******************/

    for (i=I1_STATE_1; i<=I6_STATE_1; i++) {
        alpha[i][t] = insertion1_state_prob_evaluation(t, i, hmm_ptr, &res, prev_result, from, from2, to);
    }

#ifdef profiling
    time(&end_t);
    fprintf(time_f, "t = %d, insertion' duration = %f\n", t, difftime(end_t, start_t));
    time(&start_t);
#endif
    /***********************/
    /* Non_coding state    */
    /***********************/

    alpha[R_STATE][t] = non_coding_state_prob_evaluation(t, hmm_ptr, &res, prev_result, from, to);
#ifdef profiling
    time(&end_t);
    fprintf(time_f, "t = %d, non coding duration = %f\n", t, difftime(end_t, start_t));
    time(&start_t);
#endif
    /******************/
    /* END state      */
    /******************/
    alpha[E_STATE][t] = end_state_prob_evaluation(t, hmm_ptr, &res, prev_result);

#ifdef profiling
    time(&end_t);
    fprintf(time_f, "t = %d, end duration = %f\n", t, difftime(end_t, start_t));
    time(&start_t);
#endif
    /*************************************************/
    /* START' state                                  */
    /* origianlly stop codon of genes in - strand    */
    /*************************************************/
    if (alpha[S_STATE_1][t] == 0){

        alpha[S_STATE_1][t] = max_dbl;
        path[S_STATE_1][t] = NOSTATE;


      if (t<len_seq-2 && (O[t+2] == 'A'||O[t+2] == 'a') &&
				 (((O[t] == 'T'||O[t] =='t') && (O[t+1] == 'T'|| O[t+1] == 't')) ||
					((O[t] == 'C'||O[t] =='c') && (O[t+1] == 'T'|| O[t+1] == 't')) ||
					((O[t] == 'T'||O[t] =='t') && (O[t+1] == 'C'|| O[t+1] == 'c')))) {

				alpha[S_STATE_1][t] = max_dbl;
				path[S_STATE_1][t] = R_STATE;
				alpha[S_STATE_1][t+1] = max_dbl;
				alpha[S_STATE_1][t+2] = alpha[R_STATE][t-1] - hmm_ptr->tr[TR_RS];
				path[S_STATE_1][t+1] = S_STATE_1;
				path[S_STATE_1][t+2] = S_STATE_1;

				temp_alpha = alpha[E_STATE_1][t-1] - hmm_ptr->tr[TR_ES];
				if (temp_alpha < alpha[S_STATE_1][t+2]){
					alpha[S_STATE_1][t+2] = temp_alpha;
					path[S_STATE_1][t] = E_STATE_1;
				}

				temp_alpha = alpha[E_STATE][t-1] - hmm_ptr->tr[TR_ES1];
				if (temp_alpha < alpha[S_STATE_1][t+2]){
					alpha[S_STATE_1][t+2] = temp_alpha;
					path[S_STATE_1][t] = E_STATE;
				}

                alpha[M3_STATE_1][t+2] = max_dbl;
                alpha[M6_STATE_1][t+2] = max_dbl;

				if ((O[t] == 'T'||O[t] == 't') && (O[t+1] == 'T'||O[t+1] == 't')){
					alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log54;
				}else if ((O[t] == 'C'||O[t] =='c') && (O[t+1] == 'T'||O[t+1]=='t')){
					alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log16;
				}else if((O[t] == 'T'||O[t] =='t') && (O[t+1] == 'C'||O[t+1] =='c')){
					alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log30;
				}

				/* adjustment based on probability distribution */
				start_freq=0;
				freq_id = 0;
                for(i=3; i<=60; i++){ //problems with looking to the next edge
					if (t+i+2 < len_seq){
						start_freq -= hmm_ptr->tr_S_1[i-3][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
					}
				}
				h_kd = hmm_ptr->S1_dist[2] * exp(-1*pow(start_freq-hmm_ptr->S1_dist[1],2)/(2*pow(hmm_ptr->S1_dist[0],2)));
				r_kd = hmm_ptr->S1_dist[5] * exp(-1*pow(start_freq-hmm_ptr->S1_dist[4],2)/(2*pow(hmm_ptr->S1_dist[3],2)));
				p_kd = h_kd / (h_kd + r_kd);
				if (p_kd<0.01){
					p_kd=0.01;
				}else if (p_kd>0.99){
					p_kd=0.99;
				}
				alpha[S_STATE_1][t+2] = alpha[S_STATE_1][t+2] - log(p_kd);
      }
    }
#ifdef profiling
    time(&end_t);
    fprintf(time_f, "t = %d, start' duration = %f\n", t, difftime(end_t, start_t));
    time(&start_t);
#endif

    /************************/
    /* START state          */
    /************************/
    if (alpha[S_STATE][t] == 0){

        alpha[S_STATE][t] = max_dbl;
        path[S_STATE][t] = NOSTATE;

        if (t<len_seq-2 &&  (O[t+1] == 'T'||O[t+1] =='t') && (O[t+2] == 'G'||O[t+2] =='g')&&
					((O[t] == 'A'||O[t] =='a') || (O[t] == 'G'||O[t] =='g') ||  (O[t] == 'T'||O[t] =='t'))) {

            alpha[S_STATE][t] = max_dbl;
            alpha[S_STATE][t+1] = max_dbl;
            alpha[S_STATE][t+2] = alpha[R_STATE][t-1] - hmm_ptr->tr[TR_RS];
            path[S_STATE][t] = R_STATE;
            path[S_STATE][t+1] = S_STATE;
            path[S_STATE][t+2] = S_STATE;

            temp_alpha = alpha[E_STATE][t-1] - hmm_ptr->tr[TR_ES];
            if (temp_alpha < alpha[S_STATE][t+2]){
                alpha[S_STATE][t+2] = temp_alpha;
                path[S_STATE][t] = E_STATE;
            }

            temp_alpha = alpha[E_STATE_1][t-1] - hmm_ptr->tr[TR_ES1];
            if (temp_alpha < alpha[S_STATE][t+2]){
                alpha[S_STATE][t+2] = temp_alpha;
                path[S_STATE][t] = E_STATE_1;
            }


            if ((O[t] == 'A'||O[t] =='a') ){
                alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log83;
            }else if ((O[t] == 'G'||O[t] =='g')){
                alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log(0.10);
            }else if((O[t] == 'T'||O[t] == 't')) {
                alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log07;
            }

            /* adjustment based on probability distribution */
            start_freq=0;
            freq_id = 0;
            double sub_sum = 0, lbound;
            int sub_count = 0;
            if (!prev_result /*|| t >= 30*/){
                lbound = min(30, t);

                for (i = -lbound; i <= 30; ++i){
                    if (t + i + 2 < len_seq)
                        sub_sum += hmm_ptr->tr_S[i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
                }
                sub_sum *= 61.0 / (30 + lbound + 1);
                start_freq -= sub_sum;
            } else if (prev_result) {
                lbound = min(30, t + prev_seq_len);
                for (i = -lbound; i <= 30; ++i){
                    if (t + i + 2 < len_seq) {
                        int cd1 = (t + i < 0) ? prev_O[prev_seq_len + t + i] : O[t + i];
                        int cd2 = (t + i + 1 < 0) ? prev_O[prev_seq_len + t + i + 1] : O[t + i + 1];
                        int cd3 = (t + i + 2 < 0) ? prev_O[prev_seq_len + t + i + 2] : O[t + i + 2];
                        start_freq -= hmm_ptr->tr_S[i + 30][trinucleotide(cd1, cd2, cd3)];
                    }
                }
                start_freq *= 61.0 / (30 + lbound + 1);
            }


            h_kd = hmm_ptr->S_dist[2] * exp(-1*pow(start_freq-hmm_ptr->S_dist[1],2)/(2*pow(hmm_ptr->S_dist[0],2)));
            r_kd = hmm_ptr->S_dist[5] * exp(-1*pow(start_freq-hmm_ptr->S_dist[4],2)/(2*pow(hmm_ptr->S_dist[3],2)));
            p_kd = h_kd / (h_kd + r_kd);
            if (p_kd<0.01){
                p_kd=0.01;
            }else if (p_kd>0.99){
                p_kd=0.99;
            }
            alpha[S_STATE][t+2] = alpha[S_STATE][t+2] - log(p_kd);

        }
    }
#ifdef profiling
    time(&end_t);
    fprintf(time_f, "t = %d, start duration = %f\n", t, difftime(end_t, start_t));
    time(&start_t);
#endif

    /**********************************************/
    /* END' state                                 */
    /* originally start codon of genes in - strand */
    /**********************************************/
    if (alpha[E_STATE_1][t] == 0){

        alpha[E_STATE_1][t] = max_dbl;
        path[E_STATE_1][t] = NOSTATE;

        if (t < len_seq - 2 && (O[t] == 'C'||O[t] =='c') && (O[t+1] == 'A'||O[t+1] == 'a') &&
				 ((O[t+2] == 'T'||O[t+2] =='t') || (O[t+2] == 'C'||O[t+2] =='c') || (O[t+2] == 'A'||O[t+2] =='a'))) {

            /* transition from frame6 */
            alpha[E_STATE_1][t+2] = alpha[M6_STATE_1][t-1] - hmm_ptr->tr[TR_GE];
            path[E_STATE_1][t] = M6_STATE_1;
            alpha[E_STATE_1][t] = max_dbl;
            alpha[E_STATE_1][t+1] = max_dbl;
            path[E_STATE_1][t+1] = E_STATE_1;
            path[E_STATE_1][t+2] = E_STATE_1;

            if ((O[t+2] == 'T'||O[t+2] == 't') ){
                alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log83;
            }else if ((O[t+2] == 'C'||O[t+2] =='c') ){
                alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log(0.10);
            }else if((O[t+2] == 'A'||O[t+2] =='a') ){
                alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log07;
            }

            /* adjustment based on probability distribution */
            start_freq=0;
            freq_id = 0;

            double sub_sum = 0;
            int sub_count = 0, lbound;
            if (!prev_result || t >= 30) {
                lbound = min(t, 30);
                for (i = -lbound; i <= 30; ++i){
                    if (t+i+2 < len_seq){
                        sub_sum += hmm_ptr->tr_E_1[i+30][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
                    }
                }
                sub_sum *= 61.0 / (31 + lbound);
                start_freq -= sub_sum;
            } else if (prev_result) {
                lbound = min(30, t + prev_seq_len);
                for (i = -lbound; i <= 30; ++i) {
                    if (t + i + 2 < len_seq) {
                        int cd1 = (t + i < 0) ? prev_O[prev_seq_len + t + i] : O[t + i];
                        int cd2 = (t + i + 1 < 0) ? prev_O[prev_seq_len + t + i + 1] : O[t + i + 1];
                        int cd3 = (t + i + 2 < 0) ? prev_O[prev_seq_len + t + i + 2] : O[t + i + 2];
                        start_freq -= hmm_ptr->tr_S[i + 30][trinucleotide(cd1, cd2, cd3)];
                    }
                }
                start_freq *= 61.0 / (31 + lbound);
            }
            h_kd = hmm_ptr->E1_dist[2] * exp(-1*pow(start_freq-hmm_ptr->E1_dist[1],2)/(2*pow(hmm_ptr->E1_dist[0],2)));
            r_kd = hmm_ptr->E1_dist[5] * exp(-1*pow(start_freq-hmm_ptr->E1_dist[4],2)/(2*pow(hmm_ptr->E1_dist[3],2)));
            p_kd = h_kd / (h_kd + r_kd);

            if (p_kd<0.01){
                p_kd=0.01;
            }else if (p_kd>0.99){
                p_kd=0.99;
            }
            alpha[E_STATE_1][t+2] = alpha[E_STATE_1][t+2] - log(p_kd);
        }
    }

#ifdef profiling
    time(&end_t);
    fprintf(time_f, "t = %d, end' duration = %f\n", t, difftime(end_t, start_t));
    fclose(time_f);
#endif

    if (num_N>9){
			for (i=0; i<NUM_STATE; i++){
				if (i!=R_STATE){
					alpha[i][t] = max_dbl;
					path[i][t] = R_STATE;
				}
      }
    }
  }
#ifdef viterbi_out_flg
    printf("The Viterbi matrix was evaluated\n");
    //char *fname = "../run_result/viterbi_matrix.csv";
    char fname[4096];
    sprintf(fname, "../run_result/with_graph/multiple_edge/%s-matrix.csv", head);
    FILE *f = fopen(fname, "w");
    if (!f) {
        printf("The file was not opened\n");
    }
    print_viterbi(alpha, len_seq, NUM_STATE, f);
    fclose(f);

    sprintf(fname, "../run_result/with_graph/multiple_edge/%s-path.csv", head);
    FILE *f_path = fopen(fname, "w");
    if(!f_path) {
        printf("Path file was not opened\n");
    }
    print_path(path, len_seq, NUM_STATE, f_path);
    fclose(f_path);
#endif
    res.alpha = alpha;
    res.path = path;

    res.len_seq = len_seq;
    for (i = 0; i < 6; ++i){
        res.temp_i[i] = temp_i[i];
        res.temp_i_1[i] = temp_i_1[i];
    }
    return res;
}

ViterbiResult viterbi_edge(HMM *hmm_ptr, Graph *g, size_t edge_index, int whole_genome) {
    ViterbiResult ans;
    ans.O = (char*)malloc((g->seq_len[edge_index] + 1) * sizeof (char));
    ans.O = strcpy(ans.O, g->obs_seq[edge_index]);
    ans.len_seq = g->seq_len[edge_index];
    ans.alpha = dmatrix(NUM_STATE, ans.len_seq);
    ans.path = imatrix(NUM_STATE, ans.len_seq);
    ans.temp_i = ivector(6);
    ans.temp_i_1 = ivector(6);
    int i, t;

    for (i = -1; i < NUM_STATE; ++i) {
        ans.curr_column_prev[i + 1] = -1;
        ans.first_column_prev[i + 1] = -1;
    }

    size_t n_prev= 0 ;
    for (i = 0; i < g->n_edge; ++i){
        if (g->adjacency_matrix[i][edge_index])
            ++n_prev;
    }
    int *prev_indices = ivector(n_prev);
    size_t counter = 0;
    for (i = 0; i < g->n_edge; ++i){
        if (g->adjacency_matrix[i][edge_index]){
            prev_indices[counter] = i;
            ++counter;
        }
    }
/*
#ifdef prev_index_debug
    FILE *f = fopen("../bin/prev_indices.csv", "w");
    fprintf(f, "edge_index,t");
    for (i = -1; i < NUM_STATE; ++i)
        fprintf(f, ",i=%d", (int)i);
    fprintf(f, "\n");
    fclose(f);
#endif
*/
    for (t = 0; t < ans.len_seq; ++t){

        for (i = 0; i < NUM_STATE; ++i){
            ans.alpha[i][t] = any_state_prob(hmm_ptr, t, i, &ans, g->edge_results, prev_indices, n_prev, whole_genome);
        }
    }
    free_ivector(prev_indices);
#ifdef viterbi_out_flg
    char fname[4096];
    sprintf(fname, "../run_result/with_graph/multiple_edge/%s-new-matrix.csv", g->head[edge_index]);
    FILE *f = fopen(fname, "w");
    if (!f) {
        printf("The file was not opened\n");
    }
    print_viterbi(ans.alpha, ans.len_seq, NUM_STATE, f);
    fclose(f);
#endif
    return ans;
}

void viterbi_graph(HMM *hmm_ptr, Graph* g, size_t start_index, int whole_genome) {
    //assuming that adjacency matrix is upper triangular
    g->edge_results = (ViterbiResult*)malloc(g->n_edge * sizeof (ViterbiResult));
    size_t i;
    for (i = 0; i < g->n_edge; ++i){
        g->edge_results[i] = viterbi_edge(hmm_ptr, g, i, whole_genome);
    }
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
    }else{
      gene_len = 60;
    }

    head_short = strtok(head, delimi);
    fprintf(fp_out, "%s\n", head_short); //use head_short, Ye, April 22, 2016

    /* find the state for O[N] with the highest probability */
    prob = max_dbl;
    for (i = 0; i < hmm_ptr->N; i++){

        if (alpha[i][len_seq-1] < prob){
            prob = alpha[i][len_seq-1];
            vpath[len_seq-1] = i;
        }
    }

  /* backtrack the optimal path */
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
           vpath[t]==M1_STATE_1 || vpath[t]==M4_STATE_1)){

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
            dna[dna_id]=O[t];
            dna_start_t = t + 1; //Ye April 21, 2016
            dna_f[dna_f_id]=O[t];
            //printf("Note start dna: t = %d, dna_id %d, dna_f_id %d, add %c\n", t, dna_id, dna_f_id, O[t]);
            start_orf=t+1;
            prev_match = vpath[t];

            if (vpath[t] < M6_STATE){
                codon_start=1;
            }else{
                codon_start=-1;
            }

        }else if (codon_start!=0 && (vpath[t]==E_STATE || vpath[t]==E_STATE_1 || t==len_seq-1)){

            if (vpath[t]==E_STATE || vpath[t]==E_STATE_1){
                end_t=t+3;
            }else{
                end_t=t+1;

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
            final_score = (alpha[vpath[end_t-4]][end_t-4]- alpha[vpath[start_t+2]][start_t+2] )/(end_t-start_t-5);
            frame = start_orf%3;
            if (frame==0){
                frame=3;
            }

            if (dna_id > gene_len  ){
                if (codon_start==1){
                    if(start_t == dna_start_t - 3) { //add complete start codon to dna, Ye April 21, 2016
                        strcpy(dna_tmp, dna);
                        sprintf(dna, "%c%c%c%s", O[start_t-1], O[start_t], O[start_t+1], dna_tmp);
                        //printf("add start codon to dna: %c%c%c\n", O[start_t-1], O[start_t], O[start_t+1]);
                        //printf("old dna %d %s\n", strlen(dna_tmp), dna_tmp);
                        //printf("new dna %d %s\n", strlen(dna), dna);
                    }
                    if(refine) { //add refinement of the start codons here, Ye, April 16, 2016
                        int start_old = start_t;
                        codon[0] = 0;
                        strncpy(codon, O + start_old-1, 3);
                        codon[3] = 0;
                        int s = 0;
                        //find the optimal start codon within 30bp up- and downstream of start codon
                        double e_save;
                        int s_save;
                        while((!(!strcmp(codon, "TAA") || !strcmp(codon, "TAG") || !strcmp(codon, "TGA"))) && (start_old-1-s-35>=0)) {
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
        if(s == 0) { e_save = freq_sum; s_save = s; }
        else if(freq_sum < e_save) { e_save = freq_sum; s_save = s; }
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

    }else if (codon_start!=0 &&
          ((vpath[t]>=M1_STATE && vpath[t]<=M6_STATE) ||
           (vpath[t]>=M1_STATE_1 && vpath[t]<=M6_STATE_1)) &&
          vpath[t]-prev_match<6){

      if (vpath[t] < prev_match){
    out_nt = vpath[t]+6-prev_match;
      }else{
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

    }else if (codon_start!=0 &&
          ((vpath[t]>=I1_STATE && vpath[t]<=I6_STATE) ||
           (vpath[t]>=I1_STATE_1 && vpath[t]<=I6_STATE_1))){
      dna_f_id ++;
      dna_f[dna_f_id] = tolower(O[t]);
      insert[insert_id]=t+1;
      insert_id++;

    }
    else if (codon_start!=0 && vpath[t]==R_STATE){
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
    //free_dmatrix(alpha, hmm_ptr->N);
    //free_imatrix(path, hmm_ptr->N);
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
        printf("i = %d, seq_len = %d\n%s\n", i, g->seq_len[i], g->obs_seq[i]);
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

    fscanf(fp_matr, "%zd", &res.n_edge);
    //res.adjacency_matrix = (int**)malloc(res.n_edge * sizeof(int*));
    res.adjacency_matrix = imatrix(res.n_edge, res.n_edge);
    for (i = 0; i < res.n_edge; ++i){
        //res.adjacency_matrix[i] = (int*)malloc(res.n_edge * sizeof(int));
        for (j = 0; j < res.n_edge; ++j)
            fscanf(fp_matr, "%d", &res.adjacency_matrix[i][j]);
    }

    res.seq_len = (int*)malloc(res.n_edge * sizeof (int));
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

double match_state_prob_evaluation(int t, int i, HMM *hmm_ptr, ViterbiResult *curr_res, ViterbiResult *prev_result, int whole_genome,
                        int from, int from2, int to) {
    double max_dbl = 10000000000.0;
    double ans, temp_alpha;
    int j, tmp_path, num_d;
    double **alpha = curr_res->alpha, **prev_alpha = (!prev_result) ? NULL : prev_result->alpha;
    char* O = curr_res->O, *prev_O = (!prev_result) ? NULL : prev_result->O;
    int **path = curr_res->path;
    int *temp_i = curr_res->temp_i, *prev_temp_i = (!prev_result) ? NULL : prev_result->temp_i;

    if (!prev_result && t == 0){
        alpha[i][0] = -hmm_ptr->pi[i];
        return alpha[i][0];
    }
    if (!prev_result && t > 0){
        if (i == M1_STATE) {
            //from M
            j = M6_STATE;
            ans = alpha[j][t - 1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[0][from2][to];
            tmp_path = j;

            //from D
            if (whole_genome == 0) {
                for (j = M5_STATE; j>= M1_STATE; --j){
                    num_d = i - j + 6;

                    temp_alpha = alpha[j][t - 1] -hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[0][from2][to] -
                            log(0.25) * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                    if (temp_alpha < ans) {
                        ans = temp_alpha;
                        tmp_path = j;
                    }
                }
            }

            //from S
            temp_alpha = alpha[S_STATE][t - 1] - hmm_ptr->e_M[0][from2][to];
            if (temp_alpha < ans) {
                ans = temp_alpha;
                tmp_path = S_STATE;
            }
        } else { //i = M2, ..., M6
            //from M
            j = i - 1;
            ans = alpha[j][t - 1] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[i - M1_STATE][from2][to];
            tmp_path = j;
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
                        temp_alpha = alpha[j][t - 1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[i - M1_STATE][from2][to] -
                                log(0.25) * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                        if (temp_alpha < ans) {
                            ans = temp_alpha;
                            tmp_path = j;
                        }
                    }
                }
            }

        }
        //from I
        j = (i == M1_STATE) ? I6_STATE : I1_STATE + (i - M1_STATE - 1);
        if (t < 2) {

        } else if ((i == M2_STATE || i == M5_STATE) && (O[temp_i[j - I1_STATE]] == 'T') &&
                    ((O[t] == 'A' && O[t + 1] == 'A') ||
                     (O[t] == 'A' && O[t + 1] == 'G') ||
                     (O[t] == 'G' && O[t + 1] == 'A'))) {

            } else if ((i == M3_STATE || i == M6_STATE) && temp_i[j - I1_STATE] > 0 &&
                       O[temp_i[j - I1_STATE] - 1] == 'T' &&
                       ( (O[temp_i[j - I1_STATE]] == 'A' && O[t] == 'A') ||
                         (O[temp_i[j - I1_STATE]] == 'A' && O[t] == 'G') ||
                         (O[temp_i[j - I1_STATE]] == 'G' && O[t] == 'A') )) {

            } else  {
                temp_alpha = alpha[j][t - 1] - hmm_ptr->tr[TR_IM] - log(0.25);
                if (temp_alpha < ans) {
                    ans = temp_alpha;
                    tmp_path = j;
                }
            }

        path[i][t] = tmp_path;
        alpha[i][t] = ans;
        return alpha[i][t];
    }

    if (prev_result && t == 0){
        int prev_seq_len = prev_result->len_seq;
        if (i == M1_STATE) {
            //from M
            j = M6_STATE;
            ans = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[0][from2][to];
            tmp_path = j;
            printf("From M, i = %d, j = %d\n", i, j);

            //from D
            if (!whole_genome){
                for (j = M5_STATE; j >= M1_STATE; --j){
                    num_d = i - j + 6;
                    temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[0][from2][to] -
                                log(0.25) * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                    if (temp_alpha < ans) {
                        ans = temp_alpha;
                        tmp_path = j;
                        printf("From D, i = %d, j = %d\n", i, j);

                    }
                }
            }
            //from S
            j = S_STATE;
            temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->e_M[0][from2][to];
            if (temp_alpha < ans) {
                ans = temp_alpha;
                tmp_path = j;
                printf("From S, i = %d, j = %d\n", i, j);

            }
        } else{
            /// i = M2, ..., M6
            // from M
            j = i - 1;
            ans = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[i-M1_STATE][from2][to];
            tmp_path = j;
            printf("From M, i = %d, j = %d\n", i, j);

            // from D
            if (!whole_genome) {
                for (j = M6_STATE; j >= M1_STATE; --j) {
                    num_d = (j >= i) ? (i - j + 6) : ( (j + 1 < i) ? i - j : -10 );
                    if (num_d > 0) {
                        temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[i- M1_STATE][from2][to] -
                                log(0.25) * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                        if (temp_alpha < ans) {
                            ans = temp_alpha;
                            tmp_path = j;
                            printf("From D, i = %d, j = %d\n", i, j);
                        }
                    }
                }
            }
        }
        //from I // still t = 0
        j = (i == M1_STATE) ? I6_STATE : (I1_STATE + i - M1_STATE - 1);//I1_STATE + (i - M1_STATE - 1);
        curr_res->temp_i[j - I1_STATE] = -prev_temp_i[j - I1_STATE];
        if ((i == M2_STATE || i == M5_STATE) && prev_O[prev_temp_i[j - I1_STATE]] == 'T' &&
                ((O[0] == 'A' && O[1] == 'A') ||
                 (O[0] == 'A' && O[1] == 'G') ||
                 (O[0] == 'G' && O[1] == 'A' ))){

        } else if ((i == M3_STATE || i == M6_STATE) && prev_O[prev_temp_i[j - I1_STATE] - 1] &&
                   ( (prev_O[prev_temp_i[j - I1_STATE]] == 'A' && O[0] == 'A') ||
                     (prev_O[prev_temp_i[j - I1_STATE]] == 'A' && O[0] == 'G') ||
                     (prev_O[prev_temp_i[j - I1_STATE]] == 'G' && O[0] == 'A') )) {

        } else {
            temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_IM] - log(0.25);
            if (temp_alpha < ans) {
                ans = temp_alpha;
                tmp_path = j;
                printf("From I, i = %d, j = %d\n", i, j);
            }
        }
        curr_res->path[i][0] = tmp_path;
        curr_res->alpha[i][0] = ans;
        return ans;
    }

    if (prev_result && t > 0){
        if (i == M1_STATE) {
            //from M
            j = M6_STATE;
            ans = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[0][from2][to];
            tmp_path = j;

            //from D
            if (!whole_genome) {
                for (j = M5_STATE; j >= M1_STATE; --j){
                    num_d = i - j + 6;
                    temp_alpha = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[0][from2][to] -
                            (num_d - 1) * log(0.25) - (num_d - 2) * hmm_ptr->tr[TR_DD] - hmm_ptr->tr[TR_DM];
                    if (temp_alpha < ans) {
                        ans = temp_alpha;
                        tmp_path = j;
                    }
                }
            }
            //from S
            j = S_STATE;
            temp_alpha = curr_res->alpha[j][t - 1] - hmm_ptr->e_M[0][from2][to];
            if (temp_alpha < ans) {
                ans = temp_alpha;
                tmp_path = j;
            }
        } else { //i = M2, ..., M6
            j = i - 1;
            ans = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[i - M1_STATE][from2][to];
            tmp_path = j;
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
                        temp_alpha = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[i - M1_STATE][from2][to] -
                                log(0.25) * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                        if (temp_alpha < ans) {
                            ans = temp_alpha;
                            tmp_path = j;
                        }
                    }
                }
            }

        }
        //from I
        j = (i == M1_STATE) ? I6_STATE :I1_STATE + (i - M1_STATE - 1);
        /// if temp_i < 0 then previous match was on previous edge, else prev. match was at current edge
        char prev_m = (curr_res->temp_i[j - I1_STATE] > 0) ? O[curr_res->temp_i[j - I1_STATE]] : prev_result->O[-curr_res->temp_i[j - I1_STATE]];
        char prev_m2 = (curr_res->temp_i[j - I1_STATE] > 0) ? O[curr_res->temp_i[j - I1_STATE] - 1] : prev_result->O[-curr_res->temp_i[j - I1_STATE] - 1];
        if ((i == M2_STATE || i == M5_STATE) && (t < curr_res->len_seq - 1) && prev_m == 'T' &&
                ( (O[t] == 'A' && O[t + 1] == 'A') ||
                  (O[t] == 'A' && O[t + 1] == 'G') ||
                  (O[t] == 'G' && O[t + 1] == 'A') )) {

        } else if ((i == M3_STATE || i == M6_STATE) && prev_m2 == 'T' &&
                   ((prev_m == 'A' && O[t] == 'A') ||
                    (prev_m == 'A' && O[t] == 'G') ||
                    (prev_m == 'G' && O[t] == 'A')) ) {

        } else {
            temp_alpha = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_IM] - log(0.25);
           if (temp_alpha < ans) {
               ans = temp_alpha;
               tmp_path = j;
           }
        }
        curr_res->alpha[i][t] = ans;
        curr_res->path[i][t] = tmp_path;
        return ans;

    }
    return max_dbl + 1;
}

double match1_state_prob_evaluation(int t, int i, HMM *hmm_ptr, ViterbiResult *curr_res, ViterbiResult *prev_result, int whole_genome,
                                    int from, int from2, int to) {
    if (t == 0 && !prev_result) {
        curr_res->alpha[i][0] = -hmm_ptr->pi[i];
        return curr_res->alpha[i][0];
    }
    double max_dbl = 10000000000.0;
    double ans, temp_alpha, alpha1;
    int j, tmp_path, num_d, prev_seq_len = (prev_result) ? prev_result->len_seq : -1;
    double **alpha = curr_res->alpha, **prev_alpha = (!prev_result) ? NULL : prev_result->alpha;
    char* O = curr_res->O, *prev_O = (!prev_result) ? NULL : prev_result->O;
    int **path = curr_res->path;
    int *temp_i_1 = curr_res->temp_i_1, *prev_temp_i_1 = (!prev_result) ? NULL : prev_result->temp_i_1;

    char prev1, prev2, prev3;
    prev1 = (t >= 1) ? O[t - 1] : (prev_result ? prev_O[prev_seq_len - 1] : -1);
    prev2 = (t >= 2) ? O[t - 2] : (prev_result ? prev_O[prev_seq_len + (t - 2)] : -1);
    prev3 = (t >= 3) ? O[t - 3] : (prev_result ? prev_O[prev_seq_len + (t - 3)] : -1);


    if ((i == M1_STATE_1 || i == M4_STATE_1) && (t >= 3 || prev_result) &&
            ( (prev3 == 'T' && prev2 == 'T' && prev1 == 'A') ||
              (prev3 == 'C' && prev2 == 'T' && prev1 == 'A') ||
              (prev3 == 'T' && prev2 == 'C' && prev1 == 'A') )){
        //from S
        alpha1 = (t > 0) ? alpha[S_STATE_1][t - 1] : prev_alpha[S_STATE_1][prev_seq_len - 1];
        alpha[i][t] = alpha1 - hmm_ptr->e_M_1[i - M1_STATE_1][from2][to];
        path[i][t] = S_STATE_1;
        return alpha[i][t];
    } else {
        if (i == M1_STATE_1) {
            //from M
            j = M6_STATE_1;
            //case t = 0 and prev_result is null considered above, so if t == 0 then prev_result is not null
            alpha1 = (t == 0) ?  prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
            ans = alpha1 - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M_1[0][from2][to];
            tmp_path = j;

            //from D
            if (whole_genome == 0) {
                for (j = M5_STATE_1; j >= M1_STATE_1; --j){
                    num_d = i - j + 6;
                    alpha1 = (t == 0) ?  prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
                    temp_alpha = alpha1 - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M_1[0][from2][to] -
                            log(0.25) * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                    if (temp_alpha < ans) {
                        ans = temp_alpha;
                        tmp_path = j;
                    }
                }
            }
        } else {
            //from M
            j = i - 1;
            alpha1 = (t == 0) ?  prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
            ans = alpha1 - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M_1[i - M1_STATE_1][from2][to];
            tmp_path = j;

            //from D
            if (!whole_genome){
                for (j = M6_STATE_1; j >= M1_STATE_1; --j) {
                    if (j >= i) {
                        num_d = i - j + 6;
                    } else if (j + 1 < i) {
                        num_d = i - j;
                    } else {
                        num_d = -10;
                    }

                    if (num_d > 0){
                        alpha1 = (t == 0) ?  prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
                        temp_alpha = alpha1 - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M_1[i - M1_STATE_1][from2][to] -
                                log(0.25) * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                        if (temp_alpha < ans) {
                            ans = temp_alpha;
                            tmp_path = j;
                        }
                    }
                }
            }
        }
        //from I
        j = (i == M1_STATE_1) ? I6_STATE_1 : (I1_STATE_1 + (i - M1_STATE_1 - 1));
        if (t < 2 && !prev_result){

        } else {
            prev1 = (!prev_result) ? O[temp_i_1[j - I1_STATE_1]] :
                    ( (t == 0) ? prev_O[prev_temp_i_1[j - I1_STATE_1]] :
                        (temp_i_1[j - I1_STATE_1] >= 0 ? O[temp_i_1[j - I1_STATE_1]] : prev_O[-temp_i_1[j - I1_STATE_1]]));
            if ((i == M2_STATE_1 || i == M5_STATE_1) && (t + 1 < curr_res->len_seq) && O[t + 1] == 'A' &&
                    ( (prev1 == 'T' && O[t] == 'T') ||
                      (prev1 == 'C' && O[t] == 'T') ||
                      (prev1 == 'T' && O[t] == 'C') )) {
                //nothing
            } else {
                prev1 = (!prev_result) ? O[temp_i_1[j - I1_STATE_1]] :
                        ((t == 0) ? prev_O[prev_temp_i_1[j - I1_STATE_1]] :
                            (temp_i_1[j - I1_STATE_1] >= 0 ? O[temp_i_1[j - I1_STATE_1]] : prev_O[-temp_i_1[j - I1_STATE_1]]));
                prev2 =(!prev_result) ? O[temp_i_1[j - I1_STATE_1] - 1] :
                            ((t == 0) ? prev_O[prev_temp_i_1[j - I1_STATE_1] - 1] :
                                (temp_i_1[j - I1_STATE_1] > 0 ? O[temp_i_1[j - I1_STATE_1] - 1] :
                                    (temp_i_1[j - I1_STATE_1] < 0 ? prev_O[prev_temp_i_1[j - I1_STATE_1] - 1] :
                                        prev_O[prev_seq_len - 1]) ));
                if ((i == M3_STATE_1 || i == M6_STATE_1) && (temp_i_1[j - I1_STATE_1] > 0 || prev_result) && O[t] == 'A' &&
                        ( (prev2 == 'T' && prev1 == 'T') ||
                          (prev2 == 'C' && prev1 == 'T') ||
                          (prev2 == 'T' && prev1 == 'C') )) {
                    //nothing
                } else {
                    alpha1 = (t == 0) ?  prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
                    temp_alpha = alpha1 - hmm_ptr->tr[TR_IM] - log(0.25);
                    if (temp_alpha < ans) {
                        ans = temp_alpha;
                        tmp_path = j;
                    }
                }
            }
        }
        curr_res->alpha[i][t] = ans;
        curr_res->path[i][t] = tmp_path;
        return ans;
    }
    return max_dbl + 1;
}

double insertion_state_prob_evaluation(int t, int i, HMM *hmm_ptr, ViterbiResult *curr_res, ViterbiResult *prev_result,
                                       int from, int from2, int to) {
    double max_dbl = 10000000000.0;
    if (!prev_result && t == 0){
        curr_res->alpha[i][0] = -hmm_ptr->pi[i];
        return curr_res->alpha[i][0];
    }
    double ans, temp_alpha;
    int j, temp_path;
    if (/*!prev_result &&*/ t > 0) {
        //from I
        j = i;
        ans = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
        temp_path = j;
        //from M
        j = i - I1_STATE + M1_STATE;
        if (i == I6_STATE) {
            temp_alpha = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
        } else {
            temp_alpha = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
        }
        if (temp_alpha < ans) {
            ans = temp_alpha;
            temp_path = j;
            curr_res->temp_i[i - I1_STATE] = t - 1;
        }
        curr_res->alpha[i][t] = ans;
        curr_res->path[i][t] = temp_path;
        return ans;
    }
    /// there is previous edge and t == 0
    if (t == 0 && prev_result) {
        //from I
        j = i;
        ans = prev_result->alpha[j][prev_result->len_seq - 1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
        temp_path = j;
        curr_res->temp_i[i - I1_STATE] = -(prev_result->temp_i[i - I1_STATE]);
        //from M
        j = i - I1_STATE + M1_STATE;
        if (i == I6_STATE) {
            temp_alpha = prev_result->alpha[j][prev_result->len_seq - 1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
        } else {
            temp_alpha = prev_result->alpha[j][prev_result->len_seq - 1] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
        }
        if (temp_alpha < ans) {
            ans = temp_alpha;
            temp_path = j;

            curr_res->temp_i[i - I1_STATE] = -(prev_result->len_seq - 1);
        }
        curr_res->alpha[i][0] = ans;
        curr_res->path[i][0] = temp_path;
        return ans;
    }
    return max_dbl + 1;
}

double insertion1_state_prob_evaluation(int t, int i, HMM *hmm_ptr, ViterbiResult *curr_res, ViterbiResult *prev_result, int from, int from2, int to){
    double max_dbl = 10000000000.0;
    if (t == 0 && !prev_result){
        curr_res->alpha[i][0] = -hmm_ptr->pi[i];
        return curr_res->alpha[i][0];
    }
    double **alpha = curr_res->alpha, ans, temp_alpha, alpha1;
    int j, temp_path, prev_seq_len = (prev_result) ? prev_result->len_seq: -10;
    int ** path = curr_res->path;

    if (prev_result && t < 5) {
        //if (t == 0) {
        //    curr_res->temp_i_1[i - I1_STATE_1] = -prev_result->temp_i_1[i - I1_STATE_1];
        //}
        //from I
        j = i;
        alpha1 = (t > 0) ? alpha[j][t - 1] : prev_result->alpha[j][prev_seq_len - 1];
        ans = alpha1 - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
        temp_path = j;

        //from M
        int path3, path4, path5;
        path3 = (t >= 3) ? path[S_STATE_1][t - 3] : prev_result->path[S_STATE_1][prev_seq_len + (t - 3)];
        path4 = (t >= 4) ? path[S_STATE_1][t - 4] : prev_result->path[S_STATE_1][prev_seq_len + (t - 4)];
        path5 = (t >= 5) ? path[S_STATE_1][t - 5] : prev_result->path[S_STATE_1][prev_seq_len + (t - 5)];
        if (path3 != R_STATE && path4 != R_STATE && path5 != R_STATE) {
            j = i - I1_STATE_1 + M1_STATE_1;
            alpha1 = (t > 0) ? alpha[j][t - 1] : prev_result->alpha[j][prev_seq_len - 1];
            if (i == I6_STATE_1) {
                temp_alpha = alpha1 - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
            } else {
                temp_alpha = alpha1 - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
            }
            if (temp_alpha < ans) {
                ans = temp_alpha;
                temp_path = j;

                curr_res->temp_i_1[i - I1_STATE_1] = (t > 0) ? (t - 1) : -(prev_seq_len - 1);
            }
        }
        alpha[i][t] = ans;
        path[i][t] = temp_path;
        return ans;
    } else if (t >= 5) {
        //from I
        j = i;
        ans = alpha[j][t - 1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
        temp_path = j;

        //from M
        if (path[S_STATE_1][t - 3] != R_STATE && path[S_STATE_1][t - 4] != R_STATE && path[S_STATE_1][t - 5] != R_STATE){
            j = i - I1_STATE_1 + M1_STATE_1;
            if (i == I6_STATE_1) {
                temp_alpha = alpha[j][t - 1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
            } else {
                temp_alpha = alpha[j][t - 1] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
            }
            if (temp_alpha < ans) {
                ans = temp_alpha;
                temp_path  = j;

                curr_res->temp_i_1[i - I1_STATE_1] = t - 1;
            }
        }
        path[i][t] = temp_path;
        alpha[i][t] = ans;
        return ans;
    } else if (!prev_result && t < 5){
        path[i][t] = i;
        alpha[i][t] = alpha[i][t - 1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
        return alpha[i][t];
    }
    return max_dbl + 1;
}

double non_coding_state_prob_evaluation(int t, HMM *hmm_ptr, ViterbiResult *curr_res, ViterbiResult *prev_result, int from, int to) {
    if (!prev_result && t == 0){
        curr_res->alpha[R_STATE][0] = -hmm_ptr->pi[R_STATE];
        return curr_res->alpha[R_STATE][0];
    }
    double ans, temp_alpha, max_dbl = 10000000000.0;
    int tmp_path;
    if (prev_result && t == 0) {
        //from R
        ans = prev_result->alpha[R_STATE][prev_result->len_seq - 1] - hmm_ptr->tr[TR_RR] - hmm_ptr->tr_R_R[from][to];
        tmp_path = R_STATE;

        //from E
        temp_alpha = prev_result->alpha[E_STATE][prev_result->len_seq - 1] - hmm_ptr->tr[TR_ER];
        if (temp_alpha < ans) {
            ans = temp_alpha;
            tmp_path = E_STATE;
        }

        //from E'
        temp_alpha = prev_result->alpha[E_STATE_1][prev_result->len_seq - 1] - hmm_ptr->tr[TR_ER];
        if (temp_alpha < ans) {
            ans = temp_alpha;
            tmp_path = E_STATE_1;
        }

        ans -= log(0.95);
        curr_res->alpha[R_STATE][0] = ans;
        curr_res->path[R_STATE][0] = tmp_path;
        return curr_res->alpha[R_STATE][0];
    }
    double **alpha = curr_res->alpha;
    if (t > 0) {

        //from R
        ans = alpha[R_STATE][t - 1] - hmm_ptr->tr_R_R[from][to] - hmm_ptr->tr[TR_RR];
        tmp_path = R_STATE;

        //from E
        temp_alpha = alpha[E_STATE][t - 1] - hmm_ptr->tr[TR_ER];
        if (temp_alpha < ans) {
            ans = temp_alpha;
            tmp_path = E_STATE;
        }

        //from E'
        temp_alpha = alpha[E_STATE_1][t - 1] - hmm_ptr->tr[TR_ER];
        if (temp_alpha < ans) {
            ans = temp_alpha;
            tmp_path = E_STATE_1;
        }

        ans -= log(0.95);
        curr_res->alpha[R_STATE][t] = ans;
        curr_res->path[R_STATE][t] = tmp_path;
        return curr_res->alpha[R_STATE][t];
    }
    return max_dbl + 1;
}

double end_state_prob_evaluation(int t, HMM *hmm_ptr, ViterbiResult *curr_res, ViterbiResult *prev_result) {
    double **alpha = curr_res->alpha, ans, temp_alpha, tmp_path, max_dbl = 10000000000.0;
    int **path = curr_res->path, seq_len = curr_res->len_seq;
    int prev_seq_len = (prev_result) ? prev_result->len_seq : 0;
    char *O = curr_res->O;

    if (!prev_result && t == 0) {
        curr_res->alpha[E_STATE][0] = -hmm_ptr->pi[E_STATE];
        return curr_res->alpha[E_STATE][0];
    } else { //t > 0 || (t == 0 && prev_result)
        if (alpha[E_STATE][t] == 0){

            alpha[E_STATE][t] = max_dbl;
            path[E_STATE][t] = NOSTATE;

            if (t < seq_len - 2 && O[t] == 'T' &&
                ((O[t + 1] == 'A' && O[t + 2] == 'A') ||
                 (O[t + 1] == 'A' && O[t + 2] == 'G') ||
                 (O[t + 1] == 'G' && O[t + 2] == 'A'))){

                alpha[E_STATE][t + 2] = max_dbl;
                //transition from M6
                if (t == 0) {
                    temp_alpha = prev_result->alpha[M6_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_GE];
                } else {
                    temp_alpha = alpha[M6_STATE][t - 1] - hmm_ptr->tr[TR_GE];
                }

                if (temp_alpha < alpha[E_STATE][t + 2]) {
                    alpha[E_STATE][t + 2] = temp_alpha;
                    path[E_STATE][t] = M6_STATE;
                }

                //transition from M6
                if (t == 0){
                    temp_alpha = prev_result->alpha[M3_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_GE];
                } else {
                    temp_alpha = alpha[M3_STATE][t - 1] - hmm_ptr->tr[TR_GE];
                }

                if (temp_alpha < alpha[E_STATE][t + 2]) {
                    alpha[E_STATE][t + 2] = temp_alpha;
                    path[E_STATE][t] = M3_STATE;
                }

                alpha[E_STATE][t] = max_dbl;
                alpha[E_STATE][t+1] = max_dbl;
                path[E_STATE][t+1] = E_STATE;
                path[E_STATE][t+2] = E_STATE;

                alpha[M6_STATE][t+2] = max_dbl;
                alpha[M5_STATE][t+1] = max_dbl;
                alpha[M4_STATE][t] = max_dbl;
                alpha[M3_STATE][t+2] = max_dbl;
                alpha[M2_STATE][t+1] = max_dbl;
                alpha[M1_STATE][t] = max_dbl;

                if (O[t + 1] == 'A' && O[t + 2] == 'A') {
                    alpha[E_STATE][t + 2] -= log(0.54);
                } else if (O[t + 1] == 'A' && O[t + 2] == 'G') {
                    alpha[E_STATE][t + 2] -= log(0.16);
                } else if (O[t + 1] == 'G' && O[t + 2] == 'A') {
                    alpha[E_STATE][t + 2] -= log(0.30);
                }

                //adjustment based on probanility distribution
                double start_freq = 0, h_kd, p_kd, r_kd;
                int lbound, i;
                if (t >= 60 || !prev_result){
                    lbound = min(60, t);
                    for (i = -lbound; i <= -3; ++i) {
                        start_freq -= hmm_ptr->tr_E[i + 60][trinucleotide(O[t + i], O[t + i + 1], O[t + i + 2])];
                    }
                } else { //t < 60 && prev_result
                    lbound = min(60, t + prev_seq_len);
                    char nt1, nt2, nt3;
                    for (i = -lbound; i <= -3; ++i) {
                        nt1 = (t + i >= 0) ? O[t + i] : prev_result->O[prev_seq_len + t + i];
                        nt2 = (t + i + 1 >= 0) ? O[t + i + 1] : prev_result->O[prev_seq_len + t + i + 1];
                        nt3 = (t + i + 2 >= 0) ? O[t + i + 2] : prev_result->O[prev_seq_len + t + i + 2];
                        start_freq -= hmm_ptr->tr_E[i + 60][trinucleotide(nt1, nt2, nt3)];
                    }
                }
                start_freq *= 58.0 / (lbound - 2);

                h_kd = hmm_ptr->E_dist[2] * exp(-1*pow(start_freq-hmm_ptr->E_dist[1],2)/(2*pow(hmm_ptr->E_dist[0],2)));
                r_kd = hmm_ptr->E_dist[5] * exp(-1*pow(start_freq-hmm_ptr->E_dist[4],2)/(2*pow(hmm_ptr->E_dist[3],2)));
                p_kd = h_kd / (h_kd + r_kd);
                if (p_kd<0.01){
                  p_kd=0.01;
                }else if (p_kd>0.99){
                  p_kd=0.99;
                }
                alpha[E_STATE][t + 2] -= log(p_kd);
                return alpha[E_STATE][t];
            }
            return alpha[E_STATE][t]; //no stop codone
        }
        return alpha[E_STATE][t];//alpha[E_STATE][t] != 0
    }
    //also case when stop codon is on the edge connection
    return max_dbl + 1;
}


Edge *create_raw_edge(Graph* g){
    Edge *ans = (Edge*)malloc(sizeof(Edge));
    ans->g = g;
    return ans;
}




double any_state_prob(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int *prev_ind, int n_prev, int whole_genome){
    TmpResult ans_res, temp;
    double max_dbl = 10000000000.0;
    int to = nt2int(curr_res->O[t]), group = state2group(i);

    if (n_prev == 0) {
        switch (group){
            case M_GROUP:
                ans_res = match_state_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, whole_genome, to);
            break;
            case I_GROUP:
                ans_res = insertion_state_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, to);
            break;
            case R_GROUP:
                ans_res = non_coding_state_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, to);
            break;
            case E_GROUP:
                ans_res = end_state_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev);
            break;
            case S_GROUP:
                ans_res = start_state_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev);
                break;
            case M_GROUP_1:
                ans_res = match_state1_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, whole_genome, to);
            break;
            case I_GROUP_1:
                ans_res = insertion_state1_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev, to);
            break;
            case E_GROUP_1:
                ans_res = end_state1_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev);
            break;
            case S_GROUP_1:
                ans_res = start_state1_prob_eval(hmm_ptr, t, i, curr_res, NULL, -2, n_prev);
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
                        temp = match_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, whole_genome, to);
                    break;
                    case I_GROUP:
                        temp = insertion_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, to);
                    break;
                    case R_GROUP:
                        temp = non_coding_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, to);
                    break;
                    case E_GROUP:
                        temp = end_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev);

                    break;
                    case S_GROUP:
                        temp = start_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev);
                    break;
                    case M_GROUP_1:
                        temp = match_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, whole_genome, to);
                    break;
                    case I_GROUP_1:
                        temp = insertion_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, to);
                    break;
                    case E_GROUP_1:
                        temp = end_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev);
                    break;
                    case S_GROUP_1:
                        temp = start_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev);
                    break;
                    default:
                        return max_dbl + 1;
                }

                if (j == 0 || temp.alpha < ans_res.alpha) {
                    ans_res = temp;
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
                    ans_res = match_state_prob_eval(hmm_ptr, t, i, curr_res,prev_res, curr_prev_ind, n_prev, whole_genome, to);
                break;
                case I_GROUP:
                    ans_res = insertion_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, to);
                break;
                case R_GROUP:
                    ans_res = non_coding_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev, to);
                break;
                case E_GROUP:
                    ans_res = end_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev);

                break;
                case S_GROUP:
                    ans_res = start_state_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev);
                break;
                case E_GROUP_1:
                    ans_res = end_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev);
                break;
                case S_GROUP_1:
                    ans_res = start_state1_prob_eval(hmm_ptr, t, i, curr_res, prev_res, curr_prev_ind, n_prev);
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

    curr_res->curr_column_prev[i] = ans_res.prev_ind;
    curr_res->path[i][t] = ans_res.path;
    if (t == 0) {
        curr_res->first_column_prev[i] = ans_res.prev_ind;
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

TmpResult match_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_ind, int n_prev, int whole_genome, int to) {

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
    } else { //t > 0 || prev_res
        int from2;
        char prev1, prev2;
        //prev1 = (t >= 1) ? O[t - 1] : (prev_res ? prev_O[prev_seq_len - 1] : -1);
        //prev2 = (t >= 2) ? O[t - 2] : (prev_res ? prev_O[prev_seq_len + (t - 2)] : -1);
        //from2 = (prev2 == -1 ? 2 : nt2int(prev2)) * 4 + nt2int(prev1);

        if (i == M1_STATE) {
            //from M
            j = M6_STATE;
            if (t == 0) {
                prev_alpha = prev_res[prev_ind].alpha;
                prev_seq_len = prev_res[prev_ind].len_seq;
                prev_O = prev_res[prev_ind].O;
            } else {
                prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
                prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
            }
            from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);
            alpha1 = (t == 0) ? prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1] ;
            ans_res.alpha = alpha1 - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[0][from2][to];
            ans_res.path = j;
            ans_res.prev_ind = (t == 0) ? prev_ind : curr_res->curr_column_prev[j + 1];

            //from D
            if (whole_genome == 0) {
                for (j = M5_STATE; j>= M1_STATE; --j){
                    num_d = i - j + 6;

                    if (t == 0) {
                        prev_alpha = prev_res[prev_ind].alpha;
                        prev_seq_len = prev_res[prev_ind].len_seq;
                        prev_O = prev_res[prev_ind].O;
                    } else {
                        prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                        prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
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
                prev_seq_len = prev_res[prev_ind].len_seq;
                prev_O = prev_res[prev_ind].O;
            } else {
                prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
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
                prev_seq_len = prev_res[prev_ind].len_seq;
                prev_O = prev_res[prev_ind].O;
            } else {
                prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
                prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
            }
            from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);

            alpha1 = (t == 0) ? prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
            ans_res.alpha = alpha1 - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[i - M1_STATE][from2][to];
            ans_res.path = j;
            ans_res.prev_ind = (t == 0) ? prev_ind : curr_res->curr_column_prev[j + 1];
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
                            prev_seq_len = prev_res[prev_ind].len_seq;
                            prev_O = prev_res[prev_ind].O;
                        } else {
                            prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                            prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
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
            prev_seq_len = prev_res[prev_ind].len_seq;
            prev_O = prev_res[prev_ind].O;
        } else {
            prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
            prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
            prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
        }
        from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);

        ///
        ///Обязательно проверить этот блок!!!!
        ///
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
                if ((i == M3_STATE && i == M6_STATE) && prev2 == 'T' &&
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
    }
    return ans_res;
}


TmpResult match_state1_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, int whole_genome, int to) {
    double max_dbl = 10000000000.0, temp_alpha;
    TmpResult ans_res;
    int j, num_d;
    double **alpha = curr_res->alpha, alpha1, **prev_alpha; //**prev_alpha = prev_res ? prev_res[prev_index].alpha : NULL, alpha1;
    char* O = curr_res->O, *prev_O; //= prev_res ? prev_res[prev_index].O : NULL;
    int *temp_i_1 = curr_res->temp_i_1, *prev_temp_i_1;// = prev_res ? prev_res[prev_index].temp_i_1 : NULL;
    int prev_seq_len;// = prev_res ? prev_res[prev_index].len_seq : -1;

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
            prev_seq_len = prev_res[prev_index].len_seq;
            prev_O = prev_res[prev_index].O;
        } else {
            prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
            prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
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
                    prev_seq_len = prev_res[prev_index].len_seq;
                    prev_O = prev_res[prev_index].O;
                } else {
                    prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                    prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
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
                            prev_seq_len = prev_res[prev_index].len_seq;
                            prev_O = prev_res[prev_index].O;
                        } else {
                            prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                            prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
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
                    prev_seq_len = prev_res[prev_index].len_seq;
                    prev_O = prev_res[prev_index].O;
                } else {
                    prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                    prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
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
                            prev_seq_len = prev_res[prev_index].len_seq;
                            prev_O = prev_res[prev_index].O;
                        } else {
                            prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                            prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
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
                prev_seq_len = prev_res[prev_index].len_seq;
                prev_O = prev_res[prev_index].O;
            } else {
                int x = curr_res->curr_column_prev[j + 1];
                prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
                prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
                prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
            }
            prev_temp_i_1 = (t == 0 && curr_res->curr_column_prev[j + 1] >= 0) ? prev_res[curr_res->curr_column_prev[j + 1]].temp_i_1 : NULL;
            from2 = count_from2(t, O, curr_res->len_seq, prev_O, prev_seq_len);
            if (t < 2 && n_prev > 0) {

            } else {
                prev1 = (curr_res->curr_column_prev[j + 1] < 0) ? O[temp_i_1[j - I1_STATE_1]] :
                        ( (t == 0) ? prev_O[prev_temp_i_1[j - I1_STATE_1]] :
                            (temp_i_1[j - I1_STATE_1] >= 0 ? O[temp_i_1[j - I1_STATE_1]] : prev_O[-temp_i_1[j - I1_STATE_1]]));
                if ((i == M2_STATE_1 || i == M5_STATE_1) && (t + 1 < curr_res->len_seq) && O[t + 1] == 'A' &&
                       ( ( prev1 == 'T' && O[t] == 'T') ||
                         ( prev1 == 'C' && O[t] == 'T') ||
                         ( prev1 == 'T' && O[t] == 'C'))) {

                } else {
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

TmpResult insertion_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, int to){
    double max_dbl = 10000000000.0;
    TmpResult ans_res, temp;
    int j, from;
    if (n_prev == 0 && t == 0){
        ans_res.alpha = -hmm_ptr->pi[i];
        ans_res.path = 0;
        ans_res.prev_ind = -1;
    } else if (t == 0 && n_prev > 0) {
        ans_res.prev_ind = prev_index;
        int prev_seq_len = prev_res[prev_index].len_seq;
        //from I
        j = i;
        from = nt2int(prev_res[prev_index].O[prev_res[prev_index].len_seq - 1]);
        ans_res.alpha = prev_res[prev_index].alpha[j][prev_res[prev_index].len_seq - 1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
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
            curr_res->temp_i[j - I1_STATE] = t - 1;
        }
    } else {
        ans_res.alpha = max_dbl + 1;
    }
    return ans_res;
}

TmpResult insertion_state1_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, int to){
    double max_dbl = 10000000000.0, alpha1, temp_alpha;
    double **alpha = curr_res->alpha, **prev_alpha ;//= prev_res ? prev_res[prev_index].alpha : NULL;
    TmpResult ans_res;
    char *prev_O;// = prev_res ? prev_res[prev_index].O : NULL;
    int j, from, prev_seq_len;// = prev_res ? prev_res[prev_index].len_seq : -1;
    int **path = curr_res->path, **prev_path;// = prev_res ? prev_res[prev_index].path : NULL;
    if (n_prev == 0 && t ==0){
        ans_res.alpha = -hmm_ptr->pi[i];
        ans_res.path = 0;
        ans_res.alpha2 = 0;
        ans_res.prev_ind = -1;
    } else if (t >= 5 ||
               (n_prev > 0 &&
                curr_res->curr_column_prev[i] >= 0 &&
                curr_res->curr_column_prev[i - I1_STATE_1 + M1_STATE_1 + 1] >= 0 &&
                prev_res[curr_res->curr_column_prev[i] + 1].len_seq + t >= 5
                && prev_res[curr_res->curr_column_prev[i - I1_STATE_1 + M1_STATE_1 + 1]].len_seq + t >= 5)) {
        //from = (t == 0) ? nt2int(prev_O[prev_seq_len - 1]) : nt2int(curr_res->O[t - 1]);
        //from I
        j = i;
        if (t == 0) {
            prev_alpha = prev_res[prev_index].alpha;
            prev_seq_len = prev_res[prev_index].len_seq;
            prev_path = prev_res[prev_index].path;
            prev_O = prev_res[prev_index].O;
        } else {
            prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
            prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
            prev_path = (curr_res->curr_column_prev[j + 1] >= 0) ? prev_res[curr_res->curr_column_prev[j + 1]].path : NULL;
            prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
        }

        from = (t == 0) ? nt2int(prev_O[prev_seq_len - 1]) : nt2int(curr_res->O[t - 1]);

        alpha1 = (t == 0) ?  prev_alpha[j][prev_seq_len - 1] : alpha[j][t - 1];
        ans_res.alpha = alpha1 - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
        ans_res.path = j;
        ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[j + 1];

        //from M
        j = i - I1_STATE_1 + M1_STATE_1;
        if (t == 0) {
            prev_alpha = prev_res[prev_index].alpha;
            prev_seq_len = prev_res[prev_index].len_seq;
            prev_path = prev_res[prev_index].path;
            prev_O = prev_res[prev_index].O;
        } else {
            prev_alpha = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].alpha : NULL;
            prev_seq_len = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].len_seq : -1;
            prev_path = (curr_res->curr_column_prev[j + 1] >= 0) ? prev_res[curr_res->curr_column_prev[j + 1]].path : NULL;
            prev_O = (curr_res->curr_column_prev[j + 1] >=0) ? prev_res[curr_res->curr_column_prev[j + 1]].O : NULL;
        }
        from = (t == 0) ? nt2int(prev_O[prev_seq_len - 1]) : nt2int(curr_res->O[t - 1]);

        int path3, path4, path5;
        path3 = (t >= 3) ? path[S_STATE_1][t - 3] : prev_path[S_STATE_1][prev_seq_len + (t - 3)]; //prev path is not null because we considered a condition
        path4 = (t >= 4) ? path[S_STATE_1][t - 4] : prev_path[S_STATE_1][prev_seq_len + (t - 4)];
        path5 = (t >= 5) ? path[S_STATE_1][t - 5] : prev_path[S_STATE_1][prev_seq_len + (t - 5)];
        if (path3 != R_STATE && path4 != R_STATE && path5 != R_STATE) {

            alpha1 = (t >= 0) ? alpha[j][t - 1] : prev_alpha[j][prev_seq_len - 1];
            if (i == I6_STATE_1) {
                temp_alpha = alpha1 - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_M_I[from][to];
            } else {
                temp_alpha = alpha1 - hmm_ptr->tr[TR_II] - hmm_ptr->tr_M_I[from][to];
            }
            if (temp_alpha < ans_res.alpha){
                ans_res.alpha = temp_alpha;
                ans_res.alpha = j;
                ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[j + 1];

                if (t > 0) {
                    curr_res->temp_i[i - I1_STATE_1] = t - 1;
                } //else we take argmax over all previos edges
            }
        }
    } else  {
        ans_res.alpha = max_dbl;
        ans_res.alpha2 = 0;
        ans_res.path = 0;
        ans_res.prev_ind = -1;
    }
    return ans_res;
}

TmpResult non_coding_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev, int to) {
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
        prev_seq_len = prev_res[prev_index].len_seq;
        prev_O = prev_res[prev_index].O;

        from = nt2int(prev_O[prev_seq_len - 1]);

        ans_res.alpha = prev_alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_RR] - hmm_ptr->tr_R_R[from][to];
        ans_res.path = j;
        ans_res.prev_ind = prev_index;

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
    } else if (t > 0){
        from = nt2int(curr_res->O[t - 1]);
        //from R
        j = R_STATE;
        ans_res.alpha = curr_res->alpha[j][t - 1] - hmm_ptr->tr[TR_RR] - hmm_ptr->tr_R_R[from][to];
        ans_res.path = j;
        ans_res.prev_ind = curr_res->curr_column_prev[j + 1];

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
    } else {
        ans_res.alpha = max_dbl + 1;
        ans_res.path = -2;
        ans_res.prev_ind = -2;
    }

    return ans_res;
}

TmpResult end_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev) {
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
        double **prev_alpha;// = (prev_res) ? prev_res[prev_index].alpha : NULL;
        int prev_seq_len;// = (prev_res) ? prev_res[prev_index].len_seq : 0;
        char *prev_O;// = (prev_res) ? prev_res[prev_index].O : NULL;
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

                ans_res.alpha2 = max_dbl;
                //transition from M6
                if(t == 0) {
                    //since t == 0 we now the edge which we consider to be previos
                    prev_alpha = prev_res[prev_index].alpha;
                    prev_seq_len = prev_res[prev_index].len_seq;
                    temp_alpha = prev_alpha[M6_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_GE];
                } else {
                    temp_alpha = alpha[M6_STATE][t - 1] - hmm_ptr->tr[TR_GE];
                }
                if (temp_alpha < ans_res.alpha2) {
                    ans_res.alpha2 = temp_alpha;
                    ans_res.path = M6_STATE;
                    ans_res.prev_ind = curr_res->curr_column_prev[M6_STATE + 1];
                }
                //transition from M3
                if (t == 0){
                    prev_alpha = prev_res[prev_index].alpha;
                    prev_seq_len = prev_res[prev_index].len_seq;
                    temp_alpha = prev_alpha[M3_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_GE];
                } else {
                    temp_alpha = alpha[M3_STATE][t - 1] - hmm_ptr->tr[TR_GE];
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
                    prev_seq_len = (t == 0) ? prev_res[prev_index].len_seq : prev_res[curr_res->curr_column_prev[ans_res.path + 1]].len_seq;
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
            ans_res.path = curr_res->alpha[E_STATE][t];
            ans_res.prev_ind = curr_res->curr_column_prev[E_STATE + 1];
            ans_res.alpha2 = (t < curr_res->len_seq - 2) ? curr_res->alpha[E_STATE][t + 2] : 0;
        }
    }
    return ans_res;
}

TmpResult start_state1_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev){
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

                alpha1 = (t > 0) ? alpha[R_STATE][t - 1] : prev_res[prev_index].alpha[R_STATE][prev_res[prev_index].len_seq - 1];
                ans_res.alpha2 = alpha1 - hmm_ptr->tr[TR_RS];
                ans_res.path = R_STATE;
                ans_res.prev_ind = curr_res->curr_column_prev[R_STATE + 1];
                //не забыть в main переназначить alpha и path для кодона

                //from E'
                alpha1 = (t > 0) ? alpha[E_STATE_1][t - 1] : prev_res[prev_index].alpha[E_STATE_1][prev_res[prev_index].len_seq - 1];
                temp_alpha = alpha1 - hmm_ptr->tr[TR_ES];
                if (temp_alpha < ans_res.alpha2) {
                    ans_res.alpha2 = temp_alpha;
                    ans_res.path = E_STATE_1;
                    ans_res.prev_ind = curr_res->curr_column_prev[E_STATE_1 + 1];
                }

                //from E
                alpha1 = (t > 0) ? alpha[E_STATE][t - 1] : prev_res[prev_index].alpha[E_STATE][prev_res[prev_index].len_seq - 1];
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


TmpResult end_state1_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev){
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
        if (alpha[E_STATE_1][t] == 0) {
            ans_res.alpha = max_dbl;
            ans_res.path = NOSTATE;
            ans_res.alpha2 = 0;
            ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[i + 1];

            if (t < len_seq - 2 && O[t] == 'C' && O[t + 1] == 'A' &&
                    ( O[t + 2] == 'T' || O[t + 2] == 'A' )) {
                //transition from frame 6
                alpha1 = (t > 0) ? alpha[M6_STATE_1][t - 1] : prev_res[prev_index].alpha[M6_STATE_1][prev_res[prev_index].len_seq - 1]; //here only if t >0 || t == 0 && prev_res => no null
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
                    prev_seq_len = (t == 0) ? prev_res[prev_index].len_seq : prev_res[ans_res.prev_ind].len_seq;
                    lbound = min(t, 30);
                    char nt1, nt2, nt3;
                    for (i = -lbound; i <= 30; ++i){
                        nt1 = (t + i >= 0) ? O[t + i] : prev_O[prev_seq_len + t + i];
                        nt2 = (t + i + 1 >= 0) ? O[t + i + 1] : prev_O[prev_seq_len + t + i + 1];
                        nt3 = (t + i + 2 >= 0) ? O[t + i + 2] : prev_O[prev_seq_len + t + i + 2];
                        start_freq -= hmm_ptr->tr_E_1[i + 30][trinucleotide(nt1, nt2, nt3)];
                    }
                }
                start_freq /= 61.0 * (31 + lbound);

                double h_kd = hmm_ptr->E1_dist[2] * exp(-1*pow(start_freq-hmm_ptr->E1_dist[1],2)/(2*pow(hmm_ptr->E1_dist[0],2)));
                double r_kd = hmm_ptr->E1_dist[5] * exp(-1*pow(start_freq-hmm_ptr->E1_dist[4],2)/(2*pow(hmm_ptr->E1_dist[3],2)));
                double p_kd = h_kd / (h_kd + r_kd);

                if (p_kd < 0.01) {
                    p_kd = 0.01;
                } else if (p_kd > 0.99) {
                    p_kd = 0.99;
                }
                ans_res.alpha2 -= log(p_kd);
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

TmpResult start_state_prob_eval(HMM *hmm_ptr, int t, int i, ViterbiResult *curr_res, ViterbiResult *prev_res, int prev_index, int n_prev) {
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
                    ans_res.alpha2 = prev_res[prev_index].alpha[R_STATE][prev_res[prev_index].len_seq - 1] - hmm_ptr->tr[TR_RS];
                } else {
                    ans_res.alpha2 = alpha[R_STATE][t - 1] - hmm_ptr->tr[TR_RS];
                }
                ans_res.path = R_STATE;
                ans_res.prev_ind = (t == 0) ? prev_index : curr_res->curr_column_prev[R_STATE + 1];

                //from E state
                if (t == 0) {
                    temp_alpha = prev_res[prev_index].alpha[E_STATE][prev_res[prev_index].len_seq - 1] - hmm_ptr->tr[TR_ES];
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
                    temp_alpha = prev_res[prev_index].alpha[E_STATE_1][prev_res[prev_index].len_seq - 1] - hmm_ptr->tr[TR_ES1];
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
                    prev_seq_len = (t == 0) ? prev_res[prev_index].len_seq : prev_res[ans_res.prev_ind].len_seq;
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
