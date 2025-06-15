/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
All rights reserved.

    This file is part of primer3 software suite.

    This software suite is is free software;
    you can redistribute it and/or modify it under the terms
    of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at
    your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this software (file gpl-2.0.txt in the source
    distribution); if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "dpal.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* Valid command line arguments for ntdpal */
static const char* valid_args[] = {
    "-g", "-l", "-m", "-a", "-e", "-h", "-p", "-s", "-f1", "-f2", "-f3", NULL
};

/* Forward declaration */
int parse_ntdpal_args(int argc, const char **argv, const char *usage,
                      dpal_args *a, int *print_align_end, int *use_ambiguity_codes,
                      int *use_h_matrix, int *dpal_debug, int *dpal_only, 
                      int *seq_start_index);

int
main(int argc, const char**argv)
{
    dpal_args a;
    dpal_results r;
    const unsigned char *s1, *s2;
    int tmp_ret;
    int i;
    int print_align_end = 0; /* 
                              * Print align_end_1 and align_end_2 from
                              * dpal_results.
                              */
    int use_ambiguity_codes = 0;
    int use_h_matrix = 0;
    char mode;
    dpal_mode dp_mode = DPM_GENERAL;
    int dpal_debug = 0;
    int dpal_only = 0;
    char *endptr;
    const char *msg = 
      "\nUsage: %s [-g <gval>] [-l <lval>] [-m <mval>]\n"
      "                 [-f2] [-p] [-s] [-e] <seq1> <seq2> <mode>"
      "\n\nwhere\n\n"
      "<gval> and <lval> are (positive) floats (.01 precision)\n"
      "    specifying penalties for creating or lengthening a gap\n"
      "    respectively (the penalties are subtracted from the\n"
      "    output score).\n\n"
      "-a causes the scoring matrix to be modified by dpal_set_ambiguity_codes.\n\n"
      "-e causes the end postion of the alignment in both sequences to\n"
      "   be printed.  Do not confuse with the 'e' <mode>.\n\n"
      "-f1, -f2, f3\n"
      "    force specific implementations.\n"
      "    -f2 forces use an implementation that might provide more\n"
      "    informative error messages, possibly at the expense\n"
      "    of some speed.\n\n"
      "-h use a different scoring matrix: G and C matches = 3, A and T = 2,\n"
      "   and mismatches = -0.5.\n"
      "   (The default scoring matrix assigns 1 to a match,\n"
      "   and -1 to a mismatch.)\n\n"
      "-p causes the alignment to be displayed on stderr.\n\n"
      "-s causes _only_ the score to printed.\n\n"
      "<mval> is the maximum allowed gap (default is 3).\n\n"
      "<seq1> and <seq2> are the sequences to be aligned.\n\n"
      "<mode> is one of g, G, l, or L specifying a global,\n"
      "       global end-anchored, local, or local end-achored\n"
      "       alignment respectively.  For backward compatibility\n"
      "       e is equivalent to G.\n\n";

    int seq_start_index = 0;
    
    if (argc < 4) {
      tmp_ret = fprintf(stderr, msg, argv[0]);
      exit(-1);
      return tmp_ret;
    }
    dpal_set_default_nt_args(&a);
    
    /* Parse command line arguments */
    if (parse_ntdpal_args(argc, argv, msg, &a, &print_align_end, 
                          &use_ambiguity_codes, &use_h_matrix, 
                          &dpal_debug, &dpal_only, &seq_start_index) != 0) {
      return -1;
    }
    
    /* Update i for backward compatibility */
    i = seq_start_index;
    if (use_h_matrix) dpal_set_h_nt_matrix(&a);
    if (use_ambiguity_codes) dpal_set_ambiguity_code_matrix(&a);
    if (dpal_debug && dpal_only) {
        fprintf(stderr, msg, argv[0]);
        exit(-1);
    } else {
      if (dpal_debug) {
        dp_mode = DPM_DEBUG;
      } else if (dpal_only) {
        dp_mode = DPM_FAST;
      } else {
        dp_mode = DPM_GENERAL;
      } 
    }
    a.flag = -1;
    
    if ((argc - i) != 3){
      /* Not enough remaining arguments. */
      fprintf(stderr, msg, argv[0]);
      exit(-1);
    }
    s1 = (unsigned char *) argv[i];
    s2 = (unsigned char *) argv[i+1];
    mode = *argv[i+2];
    if ('l' == mode)
        a.flag = DPAL_LOCAL;
    else if ('e' == mode || 'G' == mode)
        a.flag = DPAL_GLOBAL_END;
    else if ('g' == mode)
        a.flag = DPAL_GLOBAL;
    else if ('L' == mode)
        a.flag = DPAL_LOCAL_END;
    else {
        fprintf(stderr, msg, argv[0]);
        exit(-1);
    }
    if(print_align_end == 1) a.force_long_generic = 1;

    dpal(s1, s2, &a, dp_mode, &r);
    if (r.score == DPAL_ERROR_SCORE) {
      tmp_ret = fprintf(stderr, "Error: %s\n", r.msg);
      exit(-1);
      return tmp_ret;
    }
    if (dp_mode == DPM_FAST) {
      printf("%.2f\n", 0.01 * r.score);
      if (print_align_end) {
        if(r.align_end_1 >= 0) printf("align_end_1=%d ",r.align_end_1);
        if(r.align_end_2 >= 0) printf("align_end_2=%d\n ",r.align_end_2);
      }
    } else {
      printf("|%s|  |%s| %c ", s1, s2, mode); 
      printf("score=%.2f len=%d ", (0.01 * r.score), r.path_length);
      if (print_align_end) {
        if(r.align_end_1 >= 0) printf("align_end_1=%d ",r.align_end_1);
        if(r.align_end_2 >= 0) printf("align_end_2=%d ",r.align_end_2);
      }
      for (i=0; i<r.path_length; i++)
        printf("|%d,%d", r.path[i][0],r.path[i][1]);
      printf("|\n");
    }
    return 0;
}

/* Parse command line arguments with hybrid parser that supports both exact matches and unique abbreviations */
int parse_ntdpal_args(int argc, const char **argv, const char *usage,
                      dpal_args *a, int *print_align_end, int *use_ambiguity_codes,
                      int *use_h_matrix, int *dpal_debug, int *dpal_only, 
                      int *seq_start_index) {
    char *endptr;
    int i, j;
    
    for (i = 1; i < argc; ++i) {
        const char *current_arg = argv[i];
        const char *matched_arg = NULL;
        int match_count = 0;
        
        /* Skip if not an option (doesn't start with '-') */
        if (current_arg[0] != '-') {
            *seq_start_index = i;
            return 0; /* Successfully parsed all arguments */
        }
        
        /* First check for exact match */
        for (j = 0; valid_args[j] != NULL; j++) {
            if (strcmp(current_arg, valid_args[j]) == 0) {
                matched_arg = valid_args[j];
                match_count = 1;
                break; /* Exact match found */
            }
        }
        
        /* If no exact match, check for unique abbreviation */
        if (match_count == 0) {
            for (j = 0; valid_args[j] != NULL; j++) {
                if (strncmp(current_arg, valid_args[j], strlen(current_arg)) == 0) {
                    if (match_count == 0) {
                        matched_arg = valid_args[j];
                    }
                    match_count++;
                }
            }
        }
        
        /* Handle the match result */
        if (match_count == 0) {
            /* Unknown argument */
            fprintf(stderr, "error: unknown argument '%s'\n", current_arg);
            fprintf(stderr, usage, argv[0]);
            return -1;
        } else if (match_count > 1) {
            /* Ambiguous argument */
            fprintf(stderr, "error: ambiguous argument '%s'\n", current_arg);
            fprintf(stderr, usage, argv[0]);
            return -1;
        }
        
        /* Process the matched argument */
        if (strcmp(matched_arg, "-p") == 0) {
            *dpal_debug = 1;
        } else if (strcmp(matched_arg, "-l") == 0) {
            if (i+1 >= argc) {
                fprintf(stderr, usage, argv[0]);
                return -1;
            }
            a->gapl = (int) (strtod(argv[i+1],(char **)NULL) * -100);
            i++;
        } else if (strcmp(matched_arg, "-e") == 0) {
            *print_align_end = 1;
        } else if (strcmp(matched_arg, "-a") == 0) {
            *use_ambiguity_codes = 1;
        } else if (strcmp(matched_arg, "-h") == 0) {
            *use_h_matrix = 1;
        } else if (strcmp(matched_arg, "-g") == 0) {
            if (i+1 >= argc) {
                fprintf(stderr, usage, argv[0]);
                return -1;
            }
            a->gap = (int) (strtod(argv[i+1],(char **)NULL) * -100);
            i++;
        } else if (strcmp(matched_arg, "-m") == 0) {
            if (i+1 >= argc) {
                fprintf(stderr, usage, argv[0]);
                return -1;
            }
            a->max_gap = strtol(argv[i+1], &endptr, 10);
            if ('\0' != *endptr) {
                fprintf(stderr, usage, argv[0]);
                return -1;
            }
            i++;
        } else if (strcmp(matched_arg, "-s") == 0) {
            *dpal_only = 1;
        } else if (strcmp(matched_arg, "-f1") == 0) {
            a->force_generic = 1;
        } else if (strcmp(matched_arg, "-f2") == 0) {
            a->force_long_generic = 1;
        } else if (strcmp(matched_arg, "-f3") == 0) {
            a->force_long_maxgap1 = 1;
        }
    }
    
    /* All arguments processed but no sequence found */
    *seq_start_index = argc;
    return 0;
}
