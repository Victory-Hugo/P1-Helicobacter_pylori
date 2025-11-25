#ifndef CHROMOPAINTEREM_H
#define CHROMOPAINTEREM_H


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include "zlib.h"

static const char *helpfilestring="to run: use './ChromoPainter' with following options:\n \
       -g <geno.filein>  (REQUIRED; no default)\n\
       -r <recommap.filein>  (REQUIRED; no default)\n\
       -f <donorlist.filein>  file listing breakdown of donor haps by population (REQUIRED; no default -- unless using -a switch)\n\
       -i <int>  number of EM iterations for estimating parameters (default=0)\n\
       -in  maximize over recombination scaling constant (N_e) using E-M\n\
       -ip  maximize over copying proportions using E-M\n\
       -im  maximize over mutation (emission) probabilities using E-M\n\
       -iM  maximize over global mutation (emission) probability using E-M\n\
       -s <int>  number of samples per recipient haplotype (default=10)\n\
       -n <double>  recombination scaling constant start-value (N_e; default=400000 divided by total number of haplotypes in <geno.filein>)\n\n\
       -p  specify to use prior copying probabilities in donor list file\n\
       -m <double>  specify to use mutation (emission) probabilities in donor list file (and provide self-copying mutation rate -- if -c switch is NOT specified this value will be ignored)\n\
       -M <double>  global mutation (emission) probability (default=Li & Stephen's (2003) fixed estimate)\n\
       -k <double>  specify number of expected chunks to define a 'region' (default=100)\n\
       -c  condition on own population's individuals (default is to condition only on donor haps)\n\
       -j  specify that individuals are haploid\n\
       -u  specify that data are unlinked\n\
       -a <a_1> <a_2>  condition individuals a_1 through a_2 on every other individual (use '-a 0 0' to do all inds)\n\
       -b  print-out zipped file with suffix '.copyprobsperlocus.out' containing prob each recipient copies each donor at every SNP (note: file can be quite large)\n\
       -y  do NOT print individual count numbers next to population labels in output files (only relevant if '-a' switch is used)\n\
       -o <outfile-prefix>  (default = 'geno.filein')\n\
       -J jitter SNP locations if they are invalid.  SNPs that have the same location are placed 1 basepair after the previous SNP.\n";
extern int jitter_locations;

struct data_t {
  int nhaps_startpop;
  double nind;
  int nsnps;
  double *positions;
  double *lambda;
  int **cond_chromosomes;
  int **ind_chromosomes;
};

void stop_on_error(int val);

int reading(  char **st, char *format,  void *res);

struct data_t *ReadData(FILE *fd, int ind_val);
struct data_t *ReadData2(FILE *fd, int ind_val);
struct data_t *ReadDataHap(FILE *fd, int ind_val);
struct data_t *ReadDataHap2(FILE *fd, int ind_val);
void DestroyData(struct data_t *dat);

void DestroyData2(struct data_t *dat);

void DestroyDataHap(struct data_t *dat);

void DestroyDataHap2(struct data_t *dat);

double standard_normal_samp();

double ** sampler(int * newh, int ** existing_h, int *p_Nloci, int *p_Nhaps, int *p_nchr,  double p_rhobar, double * MutProb_vec, int * allelic_type_count_vec, double * lambda, double * pos, int nsampTOT, double * copy_prob, double * copy_probSTART, int * pop_vec, int ndonorpops, double region_size, int run_num, int run_thres, int all_versus_all_ind, int haploid_ind, int unlinked_ind, int ind_val, int print_file9_ind, FILE *fout, FILE *fout3, FILE *fout9);
int loglik(int nhaps_startpop, int *p_nloci, int p_nhaps, double N_e_start, double * recom_map, double * MutProb_vec, int nsampTOT, int ndonorpops, int * ndonorhaps, double * copy_prob, double * copy_probSTART, int * pop_vec, double region_size, int EMruns, int estimate_copyprob_ind, int estimate_recom_ind, int estimate_mutation_ind, int estimate_mutationALL_ind, int recipient_cond_ind, int all_versus_all_ind, int start_val, int end_val, char *filename, char *filenameDL, int donorlist_ind, int haploid_ind, int unlinked_ind, int print_file9_ind, int indcount_suppress_ind, FILE *fout, FILE *fout2, FILE *fout3, FILE *fout4, FILE *fout5, FILE *fout6, FILE *fout7, FILE *fout8, FILE *fout9);

void usage();
int runprogram(int argc, char *argv[]);


#endif
