#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <ctype.h>
#include "hmm.h"
#include "util_lib.h"
#define viterbi_out_flg
#define I_state_debug

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

    int temp_i[6] = {0,0,0,0,0,0};
    int temp_i_1[6] = {0,0,0,0,0,0};

    int num_N=0;
    ViterbiResult res;

    char *prev_O;


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
        prev_seq_len = (!prev_result) ? 0 : strlen(prev_result->O);
        prev_O = prev_result->O;

        from = prev_O[prev_seq_len - 1];
        from0 = prev_O[prev_seq_len - 2];
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

        /*
         * M state
         */
        for (i = M1_STATE; i <= M6_STATE; ++i){
            if (i == M1_STATE){
                /*from M state*/
                j = M6_STATE;
                alpha[i][0] = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[0][from2][to];
                path[i][0] = j;

                /*from D state*/
                if (whole_genome == 0){
                    for (j = M5_STATE; j >= M1_STATE; --j){
                        if (j >= i){
                            num_d = i - j + 6;
                        } else if (j + 1 < i){
                            num_d = i - j;
                        } else {
                            num_d = -10;
                        }
                        if ( num_d > 0 ){
                            temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[0][from2][to] -
                                    log25 * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                            if ( temp_alpha < alpha[i][0] ) {
                                alpha[i][0] = temp_alpha;
                                path[i][0] = j;
                            }
                        }
                    }
                }
                // From Start state
                temp_alpha = prev_result->alpha[S_STATE][prev_seq_len - 1] - hmm_ptr->e_M[0][from2][to];
                if ( temp_alpha < alpha[i][0] ) {
                    alpha[i][0] = temp_alpha;
                    path[i][0] = S_STATE;
                }
            } else { // i = M2, ..., M6

                //from M state
                j = i - 1;
                alpha[i][0] = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[i-M1_STATE][from2][to];
                path[i][0] = j;

                //from D state
                if ( whole_genome == 0){
                    for (j = M6_STATE; j >= M1_STATE; --j){
                        if (j >= i){
                            num_d = i - j + 6;
                        } else if (j + 1 < i) {
                            num_d = i - j;
                        } else {
                            num_d = -10;
                        }

                        if (num_d > 0){
                            temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[i - M1_STATE][from2][to]
                                        - log25 * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                            if ( temp_alpha < alpha[i][0] ){
                                alpha[i][0] = temp_alpha;
                                path[i][0] = j;
                            }
                        }
                    }
                }
            }

            //from I state
            if (i == M1_STATE) {
                j = I6_STATE;
            } else {
                j = I1_STATE + (i - M1_STATE - 1);
            }

            //to avoid stop-codon
            if (prev_seq_len < 1) {

            } else if ((i == M2_STATE || i == M5_STATE) &&
                       (prev_O[prev_result->temp_i[j - I1_STATE]] == 'T') &&
                       (( O[0] == 'A'  && O[1] == 'A') ||
                        ( O[0] == 'A' && O[1] == 'G') ||
                        ( O[0] == 'G' && O[1] == 'A')) ) {

            } else if ((i == M3_STATE || i == M6_STATE) &&
                       prev_result->temp_i[j - I1_STATE] > 0 &&
                       prev_O[prev_result->temp_i[j - I1_STATE] - 1] == 'T' &&
                       ( (prev_O[prev_result->temp_i[j - I1_STATE]] == 'A' && O[0] == 'A') ||
                         (prev_O[prev_result->temp_i[j - I1_STATE]] == 'A' && O[0] == 'G') ||
                         (prev_O[prev_result->temp_i[j - I1_STATE]] == 'G' && O[0] == 'A'))){

            } else {
                temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_IM] - log25;
                if (temp_alpha < alpha[i][0]) {
                    alpha[i][0] = temp_alpha;
                    path[i][0] = j;
                }
            }

        }

        ///
        /// I state
       ///
        for (i = I1_STATE; i <= I6_STATE; ++i){
            //from I state
            j = i;
            alpha[i][0] = prev_result->alpha[j][prev_seq_len] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
            path[i][0] = j;

            //from M state
            j = i - I1_STATE + M1_STATE;
            if (i == I6_STATE){
                temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
            } else {
                temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
            }
            if (temp_alpha < alpha[i][0]){
                alpha[i][0] = temp_alpha;
                path[i][0] = j;

                temp_i[i - I1_STATE] = -1 * (prev_seq_len - 1);
            }

        }

        ///
        /// M' state
        ///

        for (i = M1_STATE_1; i <= M6_STATE_1; ++i) {
            if ((i == M1_STATE_1 || i == M4_STATE_1) && (prev_seq_len >= 3) &&
                ( (prev_O[prev_seq_len - 3] == 'T' && prev_O[prev_seq_len - 2] == 'T' && prev_O[prev_seq_len - 1] == 'A') ||
                   (prev_O[prev_seq_len - 3] == 'C' && prev_O[prev_seq_len - 2] == 'T' && prev_O[prev_seq_len - 1] == 'A') ||
                   (prev_O[prev_seq_len - 3] == 'T' && prev_O[prev_seq_len - 2] == 'C' && prev_O[prev_seq_len - 1] == 'A'))) {
                /* from Start state  since this is actually stop codon in minus strand */
                alpha[i][0] = prev_result->alpha[i][prev_seq_len - 1] - hmm_ptr->e_M_1[i-M1_STATE_1][from2][to];
                path[i][0] = S_STATE_1;
            } else {

                if (i == M1_STATE_1){
                    //from M state
                    j = M6_STATE_1;
                    alpha[i][0] = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M_1[0][from2][to];
                    path[i][0] = j;

                ///
                ///from D state
                ///
                    if (whole_genome == 0){
                        for (j = M5_STATE_1; j >= M1_STATE_1; --j){
                            if (j >= i){
                                num_d = i - j + 6;
                            } else if (j + 1 < i){
                                num_d = i - j;
                            } else {
                                num_d = -10;
                            }
                            if (num_d > 0) {
                                temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M_1[0][from2][to] -
                                        log25 * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                                if (temp_alpha < alpha[i][0]) {
                                    alpha[i][0] = temp_alpha;
                                    path[i][0] = j;
                                }
                            }
                        }
                    }
                } else {
                    //from M state
                    j = i - 1;
                    alpha[i][0] = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M_1[0][from2][to];
                    path[i][0] = j;

                    ///
                    /// from D state
                    ///
                    if (whole_genome == 0){
                        for (j = M6_STATE_1; j <= M1_STATE_1; --j){
                            if (j >= i){
                                num_d = i - j + 6;
                            } else if (j + 1 < i) {
                                num_d = i - j;
                            } else {
                                num_d = -10;
                            }
                            if (num_d > 0) {
                                temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M_1[i - M1_STATE_1][from2][to] -
                                        log25 * (num_d - 1) - hmm_ptr->tr[TR_DD] * (num_d - 2) - hmm_ptr->tr[TR_DM];
                                if (temp_alpha < alpha[i][0]) {
                                    alpha[i][0] = temp_alpha;
                                    path[i][0] = j;
                                }
                            }
                        }
                    }
                }
                ///
                ///from I state
                ///
                if (i == M1_STATE_1) {
                    j = I6_STATE_1;
                } else {
                    j = I1_STATE_1 + (i - M1_STATE_1 - 1);
                }

                //to avoid stop-codon
                if (prev_seq_len < 1) {

                } else if ((i == M2_STATE_1 || i == M5_STATE_1) && O[1] == 'A' &&
                           ( (prev_O[prev_result->temp_i_1[j - I1_STATE_1]] == 'T' && O[0] == 'T') ||
                             (prev_O[prev_result->temp_i_1[j - I1_STATE_1]] == 'C' && O[0] == 'T') ||
                              (prev_O[prev_result->temp_i_1[j - I1_STATE_1]] == 'T' && O[0] == 'C') ) ) {

                } else if ((i == M3_STATE_1 || i == M6_STATE_1) &&
                           prev_result->temp_i[j - I1_STATE_1] > 0 && O[0] == 'A' &&
                           ( ( prev_O[prev_result->temp_i_1[j - I1_STATE_1] - 1] == 'T' && prev_O[prev_result->temp_i_1[j - I1_STATE_1]] == 'T') ||
                             ( prev_O[prev_result->temp_i_1[j - I1_STATE_1] - 1] == 'C' && prev_O[prev_result->temp_i_1[j - I1_STATE_1]] == 'T') ||
                             ( prev_O[prev_result->temp_i_1[j - I1_STATE_1] - 1] == 'T' && prev_O[prev_result->temp_i_1[j - I1_STATE_1]] == 'C') ) ){

                } else {
                    temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_IM] - log25;
                    if (temp_alpha < alpha[i][0]) {
                        alpha[i][0] = temp_alpha;
                        path[i][0] = j;
                    }
                }
            }
        }

        ///
        /// I' state
        ///
        for (i = I1_STATE_1; i <= I6_STATE_1; ++i){
            ///
            ///From I state
            ///
            j = i;
            alpha[i][0] = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
            path[i][0] = j;
            ///
            ///From M state
            ///
            if (prev_seq_len > 5 &&
                prev_result->path[S_STATE_1][prev_seq_len - 3] != R_STATE &&
                prev_result->path[S_STATE_1][prev_seq_len - 4] != R_STATE &&
                prev_result->path[S_STATE_1][prev_seq_len - 5] != R_STATE) {

                j = i - I1_STATE_1 + M1_STATE_1;
                if (i == I6_STATE_1) {
                    temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
                } else {
                    temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
                }
                if (temp_alpha < alpha[i][0]) {
                    alpha[i][0] = temp_alpha;
                    path[i][0] = j;

                    temp_i_1[i - I1_STATE_1] = -(prev_seq_len - 1);
                }
            }
        }

        ///
        ///Non-coding state
        ///
        //from R_state
        alpha[R_STATE][0] = prev_result->alpha[R_STATE][prev_seq_len - 1] - hmm_ptr->tr_R_R[from][to] - hmm_ptr->tr[TR_RR];
        path[R_STATE][0] = R_STATE;

        //from E state
        temp_alpha = prev_result->alpha[E_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_ER];
        if (temp_alpha < alpha[R_STATE][0]){
            alpha[R_STATE][0] = temp_alpha;
            path[R_STATE][0] = E_STATE;
        }

        //from E' state
        temp_alpha = prev_result->alpha[R_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_ER];
        if (temp_alpha < alpha[R_STATE][0]){
            alpha[R_STATE][0] = temp_alpha;
            path[R_STATE][0] = E_STATE_1;
        }
        alpha[R_STATE][0] -= log95;


        ///
        ///END state
        ///
        if (alpha[E_STATE][0] == 0){
            alpha[E_STATE][0] = max_dbl;
            path[E_STATE][0] = NOSTATE;

            if (len_seq > 2 && O[0] == 'T' &&
                ( (O[1] == 'A' && O[2] == 'A') ||
                  (O[1] == 'A' && O[2] == 'G') ||
                  (O[1] == 'G' && O[2] == 'A'))) {

                alpha[E_STATE][2] = max_dbl;
                //transition from frame 4, frame 5, frame 6;
                temp_alpha = prev_result->alpha[M6_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_GE];
                if (temp_alpha < alpha[E_STATE][2]) {
                    alpha[E_STATE][2] = temp_alpha;
                    path[E_STATE][0] = M6_STATE;
                }

                //transition from frame 1, frame 2, frame 3
                temp_alpha = prev_result->alpha[M3_STATE][prev_seq_len - 1] - hmm_ptr->tr[TR_GE];
                if (temp_alpha < alpha[E_STATE][2]) {
                    alpha[E_STATE][2] = temp_alpha;
                    path[E_STATE][0] = M3_STATE;
                }
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

                if (O[1] == 'A' && O[2] == 'A'){
                    alpha[E_STATE][2] -= log54;
                } else if (O[1] == 'A' && O[2] == 'G') {
                    alpha[E_STATE][2] -=log16;
                } else if (O[1] == 'G' && O[2] == 'A') {
                    alpha[E_STATE][2] -= log30;
                }

                //adjustment based on probability distribution
                start_freq = 0;
                freq_id = 0;

                double sub_sum = 0;
                int sub_count = 0;

                int lower = min(60, prev_seq_len);
                for (i = -lower; i <= -3; ++i){
                    start_freq -= hmm_ptr->tr_E[i + 60][trinucleotide(prev_O[prev_seq_len + i], prev_O[prev_seq_len + i + 1], prev_O[prev_seq_len + i + 1])];
                }
                start_freq *= 58.0 /  (lower - 2);

                    h_kd = hmm_ptr->E_dist[2] * exp(-1 * pow(start_freq - hmm_ptr->E_dist[1], 2) / (2 * pow(hmm_ptr->E_dist[0], 2)));
                    r_kd = hmm_ptr->E_dist[5] * exp(-1 * pow(start_freq - hmm_ptr->E_dist[4], 2) / (2 * pow(hmm_ptr->E_dist[3], 2)));
                    p_kd = h_kd / (h_kd + r_kd);

                    if (p_kd < 0.01){
                        p_kd = 0.01;
                    } else if (p_kd > 0.99) {
                        p_kd = 0.99;
                    }
                    alpha[E_STATE][2] = alpha[E_STATE][2] - log(p_kd);
            }
        }

        ///
        /// START' state
        ///
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

                //int upper = min(len_seq, 60);
                //for (i = 3; i <= upper; ++i) {
                //    if (i + 2 < len_seq) {
                //        start_freq -= hmm_ptr->tr_S_1[i - 3][trinucleotide(O[i], O[i + 1], O[i + 2])];
                //    }
                //}
                //start_freq *= 58.0 / (upper - 2);
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

        /******************/
        /* M state        */
        /******************/

        for (i = M1_STATE; i <= M6_STATE; i++) {

            if (alpha[i][t]<max_dbl){

				if (t == 0){
                    //watch higher
				}else{
					if (i == M1_STATE){ //i == 5

						/* from M state */
						j = M6_STATE;
						alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[0][from2][to]; // 0 = M1_STATE - M1_STATE
						path[i][t] = j;

						/* from D state */
						if (whole_genome==0){ // if !whole genome, then transition to M1 is possible from any state through Deletions.  
							for (j = M5_STATE; j >= M1_STATE; j--){ 
								if (j >= i ){
									num_d = i - j + 6; // Если добавлять 6, то получается число пропусков + 1. Пусть j = M5_state = 9, тогда i - j + 6 = 2, хотя у нас всего один пропуск: M5 -> D6 -> M1. Далее это корректируется вычитанием. 
								}else if (j + 1 < i){ // Условие --- бред. Здесь j >= i, по условию в цикле for, т.к. i = M1_STATE. Зачем этот забор из if? 
									num_d = i - j;
								}else{
                                    num_d = -10; //Что это за костыль вообще? Почему не 15, не 25?
								}
								if(num_d > 0){ // num_d > 0 for any j in this loop. 
									temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[0][from2][to]    // 0 = M1_STATE - M1_STATE
										- log25*(num_d-1) - hmm_ptr->tr[TR_DD]*(num_d-2) - hmm_ptr->tr[TR_DM]; //What does log25 mean? | 
									if ( temp_alpha < alpha[i][t]){
										alpha[i][t] = temp_alpha;
										path[i][t] = j;
									}
								}
							}
						}

						/* from Start state */
						temp_alpha = alpha[S_STATE][t-1] - hmm_ptr->e_M[0][from2][to];
						if ( temp_alpha < alpha[i][t] ){
							alpha[i][t] = temp_alpha;
							path[i][t] = S_STATE;
						}

					} else {   /*i == M2-M6*/

						/* from M state */
						j = i - 1;
						alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M[i-M1_STATE][from2][to];
						path[i][t] = j;


						/* from D state */
						if (whole_genome == 0){
							for (j = M6_STATE; j >= M1_STATE; j--){
								if (j >= i ){
									num_d = i-j+6;
								}else if (j+1 < i){ // j <= i - 1 
									num_d = i-j;
								}else{
									num_d = -10;
								}
								if (num_d > 0){

									temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M[i-M1_STATE][from2][to]
											- log25*(num_d-1) - hmm_ptr->tr[TR_DD]*(num_d-2) - hmm_ptr->tr[TR_DM];
									if ( temp_alpha < alpha[i][t]){
										alpha[i][t] = temp_alpha;
										path[i][t] = j;
									}
								}
							}
						}
					}

					/* from I state */
					if (i==M1_STATE) { 
						j = I6_STATE;
					}else{ 
							j = I1_STATE + (i - M1_STATE -1); 
					}

                    /* to aviod stop codon */

					if (t<2){
                        //t = 1
                        if (prev_result) {
                            if ((i == M2_STATE || i == M5_STATE) &&
                                prev_O[-temp_i[j - I1_STATE]] == 'T' &&
                                ( (O[1] == 'A' && O[2] == 'A') ||
                                  (O[1] == 'A' && O[2] == 'G') ||
                                  (O[1] == 'G' && O[2] == 'A') )) {

                            } else if ((i == M3_STATE || i == M6_STATE) &&
                                       prev_O[-temp_i[j - I1_STATE] + 1] == 'T' &&
                                       ( (prev_O[-temp_i[j - I1_STATE]] == 'A' && O[1] == 'A') ||
                                         (prev_O[-temp_i[j - I1_STATE]] == 'A' && O[1] == 'G') ||
                                         (prev_O[-temp_i[j - I1_STATE]] == 'G' && O[1] == 'A') )) {

                            } else {
                                temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_IM] - log25;
                                if (temp_alpha < alpha[i][t]) {
                                    alpha[i][t] = temp_alpha;
                                    path[i][t] = j;
                                }
                            }
                        }

					}else if((i==M2_STATE || i==M5_STATE) && (O[temp_i[j-I1_STATE]] == 'T'||O[temp_i[j-I1_STATE]] =='t') &&
                         (( O[t] == 'A' && O[t+1] =='A') ||
                            (O[t] == 'A' && O[t+1] =='G') ||
                            (O[t] == 'G' && O[t+1] =='A'))){

                    }else if ( temp_i[j - I1_STATE] < 1 ||
                     ((i==M3_STATE || i==M6_STATE) && (O[temp_i[j-I1_STATE]-1] == 'T'||O[temp_i[j-I1_STATE]-1] =='t') &&
                     (((O[temp_i[j-I1_STATE]] == 'A'||O[temp_i[j-I1_STATE]] == 'a') && (O[t] =='A'||O[t] == 'a')) ||
                        ((O[temp_i[j-I1_STATE]] == 'A'||O[temp_i[j-I1_STATE]] == 'a') && (O[t] =='G'||O[t] == 'g')) ||
                        ((O[temp_i[j-I1_STATE]] == 'G'||O[temp_i[j-I1_STATE]] == 'g') && (O[t] =='A'||O[t] == 'a'))))){

					}else{
						temp_alpha = alpha[j][t-1]  - hmm_ptr->tr[TR_IM] - log25;
						if ( temp_alpha < alpha[i][t]){
							alpha[i][t] = temp_alpha;
							path[i][t] = j;
						}
                    }
				}
            }
        }

        /******************/
        /* I state        */
        /******************/
        for (i=I1_STATE; i<=I6_STATE; i++) {

            if (t==0){
                if (prev_result) {
                    //from I state
                    j = i;
                    alpha[i][t] = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
                    path[i][t] = j;

                    //from M state
                    j = i - I1_STATE + M1_STATE;
                    if (i == I6_STATE) {
                        temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
                    } else {
                        temp_alpha = prev_result->alpha[j][prev_seq_len - 1] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
                    }
                    if (temp_alpha < alpha[i][t]){
                        alpha[i][t] = temp_alpha;
                        path[i][t] = j;

                        temp_i[i-I1_STATE] = t-1;
                    }

                }
            }else{
#ifdef I_state_debug
                FILE *param_f;
                if (t == 1 && head[2] == '1' && head[3] == '.' && i == I1_STATE) {
                    printf("%s\n", head);
                    param_f = fopen("../run_result/with_graph/single_edge/I_tr_parsmeters.txt", "w");
                    if (!param_f) {
                        printf("Can not open file for writing transition to Insertion state parameters\n");
                        exit(0);
                    }
                    fprintf(param_f,"t = %d, i = %d\nobs. symb = %d\nfrom = %d\nto = %d\n", t, i, O[t], from, to);
                    printf("from = %d\nto = %d\n", from, to);
                    fprintf(param_f, "Model parameters:\n");
                    fprintf(param_f, "TR_II = %lf\nTR_GG = %lf\nTR_MI = %lf\n", hmm_ptr->tr[TR_II], hmm_ptr->tr[TR_GG], hmm_ptr->tr[TR_MI]);
                    fprintf(param_f, "tr_I_I[from][to] = %lf\ntr_M_I[from][to] = %lf\n", hmm_ptr->tr_I_I[from][to], hmm_ptr->tr_M_I[from][to]);
                    //fclose(param_f);
                }

#endif
				/* from I state */
				j = i;
				alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
				path[i][t] = j;

				/* from M state */
				j = i - I1_STATE + M1_STATE ; // 5 <= j <= 10, so j is match state to every insertion; I1_STATE = 17, I6_STATE = 22, M1_STATE = 5;
				if (i==I6_STATE){
					temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
				}else{
					temp_alpha = alpha[j][t-1]  - hmm_ptr->tr[TR_MI] - hmm_ptr->tr_M_I[from][to];
				}
#ifdef I_state_debug
                if (t == 1 && head[2] == '1' && head[3] == '.' && i == I1_STATE) {
                    fprintf(param_f, "i = %d, j = %d, temp_alpha = %lf\n", i, i, alpha[i][t]);
                    fprintf(param_f, "i = %d, j = %d, temp_alpha = %lf\n", i, j, temp_alpha);
                    //fclose(param_f);
                }
#endif
				if (temp_alpha < alpha[i][t]){
					alpha[i][t] = temp_alpha;
					path[i][t] = j;

					temp_i[i-I1_STATE] = t-1;
#ifdef I_state_debug
                if (t == 1 && head[2] == '1' && head[3] == '.' && i == I1_STATE) {
                    fprintf(param_f, "alpha = %lf, path = %d", alpha[i][1], path[i][1]);
                    fclose(param_f);
                }
#endif
				}
            }
        }

        /******************/
        /* M' state        */
        /******************/

        for (i=M1_STATE_1; i<=M6_STATE_1; i++)   {
            if  ((i==M1_STATE_1 || i==M4_STATE_1) &&
                ( ( t>=3 &&
                    (((O[t-3] == 'T'||O[t-3] == 't') && (O[t-2] == 'T'||O[t-2] == 't') && (O[t-1] == 'A'||O[t-1] =='a')) ||
                    ((O[t-3] == 'C'||O[t-3] == 'c') && (O[t-2] == 'T'||O[t-2] == 't') && (O[t-1] == 'A'||O[t-1] =='a')) ||
                    ((O[t-3] == 'T'||O[t-3] == 't') && (O[t-2] == 'C'||O[t-2] == 'c') && (O[t-1] == 'A'||O[t-1] =='a'))))
                || (t == 2 && prev_O &&
                    ( (prev_O[prev_seq_len - 1] == 'T' && O[0] == 'T' && O[1] == 'A') ||
                    (prev_O[prev_seq_len - 1] == 'C' && O[0] == 'T' && O[1] == 'A') ||
                    (prev_O[prev_seq_len - 1] == 'T' && O[0] == 'C' && O[1] == 'A') ) )
                || (t == 1 && prev_O &&
                    ( (prev_O[prev_seq_len - 2] == 'T' && prev_O[prev_seq_len - 1] == 'T' && O[0] == 'A') ||
                    (prev_O[prev_seq_len - 2] == 'C' && prev_O[prev_seq_len - 1] == 'T' && O[0] == 'A') ||
                    (prev_O[prev_seq_len - 2] == 'T' && prev_O[prev_seq_len - 1] == 'C' && O[0] == 'A') ) ) ) )  {

				/* from Start state  since this is actually stop codon in minus strand */
				alpha[i][t] = alpha[S_STATE_1][t-1] - hmm_ptr->e_M_1[i-M1_STATE_1][from2][to];
				path[i][t] = S_STATE_1;

            }else{

                if (t==0){
				
            }else{

                if (i==M1_STATE_1 ){

                    /* from M state */
                    j = M6_STATE_1;
                    alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M_1[0][from2][to];
                    path[i][t] = j;

                    /* from D state */
                    if (whole_genome==0){
                        for (j=M5_STATE_1; j>=M1_STATE_1; j--){
                            if (j >= i){
                                num_d = i-j+6;
                            }else if (j+1 <i){
                                num_d = i-j;
                            } else {
                                num_d = -10;
                            }
                            if (num_d > 0){
                                temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M_1[0][from2][to]
                                                            - log25*(num_d-1) - hmm_ptr->tr[TR_DD]*(num_d-2) - hmm_ptr->tr[TR_DM];
                                if ( temp_alpha < alpha[i][t]){
                                    alpha[i][t] = temp_alpha;
                                    path[i][t] = j;
                                }
                            }
                        }
                    }

                }else{

                    /* from M state */
                    j = i - 1;
                    alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_MM] - hmm_ptr->e_M_1[i-M1_STATE_1][from2][to];
                    path[i][t] = j;

                    /* from D state */
                    if (whole_genome==0){
                        for (j=M6_STATE_1; j>=M1_STATE_1; j--){
                            if (j >= i ){
                                num_d = i-j+6;
                            }else if (j+1 < i){
                                num_d = i-j;
                            }else{
                                num_d = -10;
                            }
                            if (num_d>0){
                                temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_MD] - hmm_ptr->e_M_1[i-M1_STATE_1][from2][to]
                                    - log25*(num_d-1) - hmm_ptr->tr[TR_DD]*(num_d-2) - hmm_ptr->tr[TR_DM];

                                if ( temp_alpha < alpha[i][t]){
                                    alpha[i][t] = temp_alpha;
                                    path[i][t] = j;
                                }
                            }
                        }
                    }
                }

                /* from I state */
                if (i==M1_STATE_1) {
                    j = I6_STATE_1;

                }else{
                    j = I1_STATE_1 + (i - M1_STATE_1 -1);
                }


                /* to aviod stop codon */
                if (t<2){
                    //t = 1
                    if (prev_result) {
                        if ((i == M2_STATE_1 || i == M5_STATE_1) && O[2] == 'A' &&
                                ( (prev_O[-temp_i[j - I1_STATE_1]] == 'T' && O[1] == 'T') ||
                                  (prev_O[-temp_i[j - I1_STATE_1]] == 'C' && O[1] == 'T') ||
                                  (prev_O[-temp_i[j - I1_STATE_1]] == 'T' && O[1] == 'C') )) {
                        } else if ((i == M3_STATE_1 || i == M6_STATE_1) && O[1] == 'A' &&
                                       ( ( prev_O[-temp_i[j - I1_STATE_1] - 1] == 'T' && prev_O[-temp_i[j - I1_STATE_1]] == 'T') ||
                                         ( prev_O[-temp_i[j - I1_STATE_1] - 1] == 'C' && prev_O[-temp_i[j - I1_STATE_1]] == 'T') ||
                                         ( prev_O[-temp_i[j - I1_STATE_1] - 1] == 'T' && prev_O[-temp_i[j - I1_STATE_1]] == 'C') ) ){

                            } else {
                                temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_IM] - log25;
                                if (temp_alpha < alpha[i][t]) {
                                    alpha[i][t] = temp_alpha;
                                    path[i][t] = j;
                                }
                        }

                    }
                } else  if ((i==M2_STATE_1 || i==M5_STATE_1) && (O[t+1] == 'A'||O[t+1] == 'a') &&
					   (((O[temp_i_1[j-I1_STATE_1]] == 'T'|| O[temp_i_1[j-I1_STATE_1]] == 't') && (O[t] =='T'|| O[t] =='t')) ||
							((O[temp_i_1[j-I1_STATE_1]] == 'C'|| O[temp_i_1[j-I1_STATE_1]] == 'c') && (O[t] =='T'|| O[t] =='t')) ||
							((O[temp_i_1[j-I1_STATE_1]] == 'T'|| O[temp_i_1[j-I1_STATE_1]] == 't') && (O[t] =='C'|| O[t] =='c')))){

                    }else if ( temp_i_1[j - I1_STATE_1] > 0 &&
                            (i==M3_STATE_1 || i==M6_STATE_1) && (O[t] == 'A'||O[t] == 'a') &&
                            (((O[temp_i_1[j-I1_STATE_1]-1] == 'T'|| O[temp_i_1[j-I1_STATE_1]-1]=='t') &&
							 (O[temp_i_1[j-I1_STATE_1]] =='T'|| O[temp_i_1[j-I1_STATE_1]] =='t')) || 
							((O[temp_i_1[j-I1_STATE_1]-1] == 'C'|| O[temp_i_1[j-I1_STATE_1]-1]=='c') &&
							 (O[temp_i_1[j-I1_STATE_1]] =='T'|| O[temp_i_1[j-I1_STATE_1]] =='t')) || 
							((O[temp_i_1[j-I1_STATE_1]-1] == 'T'|| O[temp_i_1[j-I1_STATE_1]-1]=='t') &&
							 (O[temp_i_1[j-I1_STATE_1]] =='C'|| O[temp_i_1[j-I1_STATE_1]] =='c')))){
					}else {

						temp_alpha = alpha[j][t-1]  - hmm_ptr->tr[TR_IM] - log25;
						if ( temp_alpha < alpha[i][t]){
							alpha[i][t] = temp_alpha;
							path[i][t] = j;
						}
					}
				}
            }
        }

    /******************/
    /* I' state        */
    /******************/
        for (i=I1_STATE_1; i<=I6_STATE_1; i++) {

            if (t==0){

            }else{
				/* from I state */
                j = i;
				alpha[i][t] = alpha[j][t-1] - hmm_ptr->tr[TR_II] - hmm_ptr->tr_I_I[from][to];
				path[i][t] = j;

				/* from M state */
                if (t > 4) {
                    if (path[S_STATE_1][t-3] != R_STATE && path[S_STATE_1][t-4] != R_STATE && path[S_STATE_1][t-5] != R_STATE){
                        j = i - I1_STATE_1 + M1_STATE_1;
                        if (i==I6_STATE_1){
                            temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] -hmm_ptr->tr_M_I[from][to];
                        }else{
                            temp_alpha = alpha[j][t-1]  - hmm_ptr->tr[TR_MI] -hmm_ptr->tr_M_I[from][to];
                        }
                        if (temp_alpha < alpha[i][t]){
                            alpha[i][t] = temp_alpha;
                            path[i][t] = j;

                            temp_i_1[i-I1_STATE_1] = t-1;
                        }
                    }
                } else if (prev_result) {
                    int flag = 0;
                    switch (t){
                    case 4:
                        if (path[S_STATE_1][1] != R_STATE &&
                                path[S_STATE_1][0] != R_STATE &&
                                prev_result->path[S_STATE_1][prev_seq_len - 1] != R_STATE) {

                            flag = 1;
                        }
                        break;
                    case 3:
                        if (prev_seq_len > 1 && path[S_STATE_1][0] != R_STATE &&
                                prev_result->path[S_STATE_1][prev_seq_len - 1] != R_STATE &&
                                prev_result->path[S_STATE_1][prev_seq_len - 2] != R_STATE) {
                            flag = 1;
                        }
                        break;
                    case 2:
                        if (prev_seq_len > 2 && prev_result->path[S_STATE_1][prev_seq_len - 1] != R_STATE &&
                                prev_result->path[S_STATE_1][prev_seq_len - 2] != R_STATE &&
                                prev_result->path[S_STATE_1][prev_seq_len - 3] != R_STATE) {
                            flag = 1;
                        }
                        break;
                    case 1:
                        if (prev_seq_len > 3 && prev_result->path[S_STATE_1][prev_seq_len - 2] != R_STATE &&
                                prev_result->path[S_STATE_1][prev_seq_len - 3] != R_STATE &&
                                prev_result->path[S_STATE_1][prev_seq_len - 4] != R_STATE) {
                            flag = 1;
                        }

                        break;
                    }
                    if (flag == 1){
                        j = i - I1_STATE_1 + M1_STATE_1;
                        if (i==I6_STATE_1){
                            temp_alpha = alpha[j][t-1] - hmm_ptr->tr[TR_GG] - hmm_ptr->tr[TR_MI] -hmm_ptr->tr_M_I[from][to];
                        }else{
                            temp_alpha = alpha[j][t-1]  - hmm_ptr->tr[TR_MI] -hmm_ptr->tr_M_I[from][to];
                        }
                        if (temp_alpha < alpha[i][t]){
                            alpha[i][t] = temp_alpha;
                            path[i][t] = j;

                            temp_i_1[i-I1_STATE_1] = t-1;
                        }
                    }
                }
            }
        }

    /***********************/
    /* Non_coding state    */
    /***********************/

        if (t==0){
    
		}else{
            alpha[R_STATE][t] = alpha[R_STATE][t-1] - hmm_ptr->tr_R_R[from][to] - hmm_ptr->tr[TR_RR];
            path[R_STATE][t] = R_STATE;

            temp_alpha = alpha[E_STATE][t-1]  - hmm_ptr->tr[TR_ER];
            if (temp_alpha < alpha[R_STATE][t] ){
				alpha[R_STATE][t] = temp_alpha;
				path[R_STATE][t] = E_STATE;
            }

            temp_alpha = alpha[E_STATE_1][t-1] - hmm_ptr->tr[TR_ER] ;
            if (temp_alpha < alpha[R_STATE][t] ){
                alpha[R_STATE][t] = temp_alpha;
                path[R_STATE][t] = E_STATE_1;
            }
            alpha[R_STATE][t] -= log95;
        }

    /******************/
    /* END state      */
    /******************/
    if (alpha[E_STATE][t] == 0){ // Check if previously we haven't considered state number t to be a part of stop codon? 

        alpha[E_STATE][t] = max_dbl;
        path[E_STATE][t] = NOSTATE;

        if (t < len_seq -2 && (O[t] == 'T'||O[t] == 't')  &&
				 (((O[t+1] == 'A'||O[t+1] == 'a') && (O[t+2] == 'A'||O[t+2] =='a')) ||
					((O[t+1] == 'A'||O[t+1] == 'a') && (O[t+2] == 'G'||O[t+2] =='g')) ||
					((O[t+1] == 'G'||O[t+1] == 'g') && (O[t+2] == 'A'||O[t+2] =='a')))) {

            alpha[E_STATE][t+2] = max_dbl;
            /* transition from frame4,frame5,and frame6 */             //Frame numerations here and in the article are different.
            temp_alpha = alpha[M6_STATE][t-1] - hmm_ptr->tr[TR_GE];
            if (temp_alpha < alpha[E_STATE][t+2]){
                alpha[E_STATE][t+2] = temp_alpha;
                path[E_STATE][t] = M6_STATE;
            }

            /* transition from frame1,frame2,and frame3 */
            temp_alpha  = alpha[M3_STATE][t-1] - hmm_ptr->tr[TR_GE];
            if (temp_alpha < alpha[E_STATE][t+2]){
                alpha[E_STATE][t+2] = temp_alpha;
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

            if ((O[t+1] == 'A'||O[t+1] =='a') && (O[t+2] == 'A'||O[t+2] =='a')){
                alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log54;
            }else if ((O[t+1] == 'A'||O[t+1] =='a') && (O[t+2] == 'G'||O[t+2] =='g')){
                alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log16;
            }else if((O[t+1] == 'G'||O[t+1] == 'g') && (O[t+2] == 'A'||O[t+2] =='a')) {
                alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log30;
            }

            /* adjustment based on probability distribution */
            start_freq=0;
            freq_id = 0;

            double sub_sum = 0;
            int sub_count = 0;

            if (!prev_result || t >= 60){ // !prev_result => no prev edge; t >= 60 => we don't need prev edge
                int lbound = min(t, 60);

                for (i = -lbound; i <= -3; ++i){
                    if (t + i + 2 < len_seq) {
                        sub_sum += hmm_ptr->tr_E[i+60][trinucleotide(O[t+i], O[t+i+1], O[t+i+2])];
                    }
                }
                sub_sum  *=  58.0 / (lbound - 2);
                start_freq -= sub_sum;
            } else  if (prev_result){
                int lbound = min(60, prev_seq_len + t);
                for (i = -lbound; i <= -3; ++i) {
                    int cd1, cd2, cd3;
                    cd1 = (t + i > 0) ? O[t + i] : prev_O[prev_seq_len + t + i];
                    cd2 = (t + i + 1 > 0) ? O[t + i + 1] : prev_O[prev_seq_len + t + i + 1];
                    cd3 = (t + i + 2 > 0) ? O[t + i + 2] : prev_O[prev_seq_len + t + 2];
                    start_freq -= hmm_ptr->tr_E[i + 60][trinucleotide(cd1, cd2, cd3)];
                }
                start_freq *= 58.0 / (lbound - 2);
            }
            h_kd = hmm_ptr->E_dist[2] * exp(-1*pow(start_freq-hmm_ptr->E_dist[1],2)/(2*pow(hmm_ptr->E_dist[0],2)));
            r_kd = hmm_ptr->E_dist[5] * exp(-1*pow(start_freq-hmm_ptr->E_dist[4],2)/(2*pow(hmm_ptr->E_dist[3],2)));
            p_kd = h_kd / (h_kd + r_kd);
            if (p_kd<0.01){
                p_kd=0.01;
            }else if (p_kd>0.99){
                p_kd=0.99;
            }
            alpha[E_STATE][t+2] = alpha[E_STATE][t+2] - log(p_kd);
        }
    }

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
    sprintf(fname, "../run_result/with_graph/single_edge/%s-matrix.csv", head);
    FILE *f = fopen(fname, "w");
    if (!f) {
        printf("The file was not opened\n");
    }
    print_viterbi(alpha, len_seq, NUM_STATE, f);
    fclose(f);
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
    res.adjacency_matrix = (int**)malloc(res.n_edge * sizeof(int*));
    for (i = 0; i < res.n_edge; ++i){
        res.adjacency_matrix[i] = (int*)malloc(res.n_edge * sizeof(int));
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
            res.seq_len[i] += sizeof(tmp_str);
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

    return;
}
