#include "ChromoPainterMutEM.h"

#include <time.h>

#define PI 3.141593

/*****************************************************************
// TAKES A PHASE-STYLE SET OF n PHASED HAPLOTYPES, SUCH THAT THE TOP x HAPLOTYPES ARE POTENTIAL DONORS TO THE BOTTOM n-x HAPLOTYPES (THESE n-x HAPLOTYPES ARE ASSUMED TO BE DIPLOID INDIVIDUALS, WITH EACH INDIVIDUAL'S TWO HAPLOTYPES ON CONSECUTIVE LINES). THEN PAINTS EACH OF THE n-x BOTTOM HAPLOTYPES AS A MOSAIC OF THE x DONOR HAPLOTYPES, WITH OR WITHOUT MAXIMIXING OVER N_e (DEFAULT) OR COPYING PROPORTIONS (IF USER SELECTS) USING E-M, OUTPUTTING AS MANY SAMPLES AS DESIRED (AFTER MAXIMIZATION, IF CHOSEN -- MAXIMIZATION IS WHERE THE DIPLOID ASSUMPTION COMES IN, AS N_E ESTIMATES OR COPYING PROPORTIONS ARE AVERAGED OVER AN INDIVIDUAL'S TWO HAPS AT EACH E-M STEP)

// to compile:  gcc -Wall -o ChromoPainter ChromoPainter.c -lm -lz

// to run: "./ChromoPainter" with following options:
//                 -g geno.filein (PHASEish-style file; no default, required)
//                 -r recommap.filein (no default, required)
//                 -f file listing breakdown of donor haps by population (no default, required)
//                 -i number of EM iterations for estimating parameters (default=0)
//                 -in maximize over recombination scaling constant (N_e) using E-M
//                 -ip maximize over copying proportions using E-M
//                 -im maximize over mutation (emission) probabilities using E-M
//                 -iM maximize over global mutation (emission) probability using E-M
//                 -s number of samples per recipient haplotype (default=10)
//                 -n recombination scaling constant (N_e; default=400000 divided by total number of haplotypes in 'geno.filein')
//                 -p specify to use prior copying probabilities in donor list file
//                 -m specify to use mutation (emission) probabilities in donor list file (and provide self-copying mutation rate -- if -c switch is NOT specified this value will be ignored)
//                 -M global mutation (emission) probability (default=Li & Stephen's (2003) fixed estimate)
//                 -k specify number of expected chunks to define a 'region' (default=100)
//                 -c condition on own population's individuals (default is to condition only on donor haps)
//                 -j specify that individuals are haploid
//                 -u specify that data are unlinked
//                 -a a_1 a_2 condition individuals a_1 through a_2 on every other individual (use '-a 0 0' to do all inds)
//                 -b print-out zipped file  with suffix '.copyprobsperlocus.out' containing prob each recipient copies each donor at every SNP (note: file can be quite large)
//                 -y do NOT print individual count numbers next to population labels in output files (only relevant if '-a' switch is used)
//                 -o outfile-prefix (default = "geno.filein")
//                 --help print help menu

// example:
//        ./ChromoPainter -g example.genotypes -r example.recommap -f example.donorlist -o example.out

*******************************************************************/
int jitter_locations;

int errormode;
void stop_on_error(int val) {
    if(errormode==0) {
        exit(val);
    }else{
        printf("Unrecoverable error, cancel execution via the GUI and reconfigure!\n");
       //# Sleep(10000);
        getchar();
        exit(val);
    }
}

int reading(char **st, char *format, void *res)
{
    int i;
    char *rs;
    rs = *st;
    for(i = 0; isspace(rs[i]); i++) ;
    if (!rs[i]) return 0;
    for(; !isspace(rs[i]); i++) ;
    if (rs[i]) rs[i++] = 0;
    if (!sscanf(*st, format, res)) return 0;
    *st += i;
    return 1;
}


/*struct data_t *ReadData(char *filename) {*/
struct data_t *ReadData(FILE *fd, int ind_val){
  struct data_t *dat;
  char * line = malloc(100000000 * sizeof(char));
  char *step;
  char waste[10];
  int i,j;
  int nhaps;

  dat=malloc(sizeof(struct data_t));
  /* if dat==NULL etc ... */

  /* Number of chromosomes from start population (to condition on)*/
  fgets(line,2047,fd);
  if (line==NULL) { printf("ReadData: error with PHASE-style input file; first (haplotype donor) line empty\n"); stop_on_error(1);}
  sscanf(line,"%d",&dat->nhaps_startpop);
  if (dat->nhaps_startpop <= 0) { printf("Number of donor haplotypes must be > 0 (or >= 0 if you choose the self-copying '-c' or all-versus-all '-a' options). Found %d. Exiting...\n",dat->nhaps_startpop); stop_on_error(1);}

  /* Number of individuals */
  fgets(line,2047,fd);
  if (line==NULL) { printf("ReadData: error with PHASE-style input file; second (individual) line empty\n"); stop_on_error(1);}
  sscanf(line,"%lf",&dat->nind);
  nhaps = (int) 2*dat->nind;
  if (dat->nind < 0) { printf("Number of total individuals must be > 0. Exiting...\n"); stop_on_error(1);}
  if ((nhaps-2) < dat->nhaps_startpop) { printf("Number of total haps must be greater than or equal to number of donor haplotypes plus one extra individual.\n"); stop_on_error(1);}

  /* Number of SNPs */
  fgets(line,2047,fd);
  if (line==NULL) { printf("ReadData: error with PHASE-style input file; third (snps() line empty\n"); stop_on_error(1);}
  sscanf(line,"%d",&dat->nsnps);
  if (dat->nsnps <= 0) { printf("Number of sites must be > 0. Exiting...\n"); stop_on_error(1);}

  dat->positions=malloc(dat->nsnps*sizeof(double));
  dat->lambda=malloc((dat->nsnps-1)*sizeof(double));
  dat->cond_chromosomes=malloc(dat->nhaps_startpop*sizeof(int *));
  dat->ind_chromosomes=malloc(2*sizeof(int *));
  for (i=0; i<dat->nhaps_startpop; i++)
    dat->cond_chromosomes[i]=malloc(dat->nsnps*sizeof(int));
  for (i=0; i<2; i++)
    dat->ind_chromosomes[i]=malloc(dat->nsnps*sizeof(int));

  /* Positions */
  fgets(line,100000000,fd);
  if (line==NULL) { printf("ReadData: error with PHASE-style input file; fourth (SNP location) line empty\n"); stop_on_error(1);}
  step=line;
  reading(&step,"%s",waste);
  for (i=0; i<dat->nsnps; i++)
    {
      reading(&step,"%lf",&dat->positions[i]);
      if (dat->positions[i] < 0) { printf("Basepair positions must be >= 0. Exiting...\n"); stop_on_error(1);}
      if(i>0) if (dat->positions[i] == dat->positions[i-1])
	{
	  if(jitter_locations==0) { 
	    printf("Basepair positions must be increasing (at SNPs %i-%i with positions %f-%f). Rerun with option \"-J\" to ignore. Exiting...\n",i-1,i,dat->positions[i-1],dat->positions[i]); stop_on_error(1);
	  }else{
	    double newloc=dat->positions[i-1]+1;
	    dat->positions[i]=newloc;
	  }
	}
      if (i < (dat->nsnps-1)) dat->lambda[i] = 1.0;
    }

  /* Position type */
  fgets(line,100000000,fd);
  if (line==NULL) { printf("ReadData: error with PHASE-style input file; fifth (mutation type) line empty\n"); stop_on_error(1);}

  for(i = 0; i < dat->nhaps_startpop; i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("error with PHASE-style input file, line length is %i, expected %i\n",(int)strlen(line),dat->nsnps+1); stop_on_error(1);}
      for (j=0; j<dat->nsnps;j++) {
	if(line[j]=='0')
	  dat->cond_chromosomes[i][j]=0;
	if(line[j]=='1')
	  dat->cond_chromosomes[i][j]=1;
	if(line[j]=='A')
	  dat->cond_chromosomes[i][j]=2;
	if(line[j]=='C')
	  dat->cond_chromosomes[i][j]=3;
	if(line[j]=='G')
	  dat->cond_chromosomes[i][j]=4;
	if(line[j]=='T')
	  dat->cond_chromosomes[i][j]=5;
	if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T'))
	  {
	    printf("Allele-type invalid for hap%d, snp%d. Exiting...\n",i+1,j+1);
	    stop_on_error(1);
	  }
      }
    }
  for (i= 0; i < (2*ind_val); i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("ReadData: error with PHASE-style input file at haplotype %i of %i\n",i, 2*ind_val); stop_on_error(1);}
   }
  for(i = 0; i < 2; i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("ReadData: error with PHASE-style input file at haplotype %i of %i (reading)\n",i, 2*ind_val); stop_on_error(1);}
      for (j=0; j<dat->nsnps;j++) {
	if(line[j]=='0')
	  dat->ind_chromosomes[i][j]=0;
	if(line[j]=='1')
	  dat->ind_chromosomes[i][j]=1;
	if(line[j]=='A')
	  dat->ind_chromosomes[i][j]=2;
	if(line[j]=='C')
	  dat->ind_chromosomes[i][j]=3;
	if(line[j]=='G')
	  dat->ind_chromosomes[i][j]=4;
	if(line[j]=='T')
	  dat->ind_chromosomes[i][j]=5;
	if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T'))
	  {
	    printf("Allele type invalid for hap%d, snp%d. Exiting...\n",dat->nhaps_startpop+2*ind_val+i+1,j+1);
	    stop_on_error(1);
	  }
      }
    }

  free(line);
  return dat;

}

struct data_t *ReadData2(FILE *fd, int ind_val){
  struct data_t *dat;
  char * line = malloc(100000000 * sizeof(char));
  char *step;
  char waste[10];
  int i,j;
  int nhaps;

  dat=malloc(sizeof(struct data_t));
  /* if dat==NULL etc ... */

  /* Number of chromosomes from start population (to condition on)*/
  fgets(line,2047,fd);
  if (line==NULL) { printf("ReadData2: error with PHASE-style input file\n"); stop_on_error(1);}
  sscanf(line,"%d",&dat->nhaps_startpop);
  if (dat->nhaps_startpop < 0) { printf("Number of donor haplotypes must be >= 0. Exiting...\n"); stop_on_error(1);}

  /* Number of individuals */
  fgets(line,2047,fd);
  if (line==NULL) { printf("ReadData2: error with PHASE-style input file\n"); stop_on_error(1);}
  sscanf(line,"%lf",&dat->nind);
  nhaps = (int) 2*dat->nind;
  if (dat->nind <= 0) { printf("Number of total individuals must be > 0. Exiting...\n"); stop_on_error(1);}
  if ((nhaps-2) < dat->nhaps_startpop) { printf("Number of total haps must be greater than or equal to number of donor haplotypes plus one extra individual.\n"); stop_on_error(1);}

  /* Number of SNPs */
  fgets(line,2047,fd);
  if (line==NULL) { printf("ReadData2: error with PHASE-style input file\n"); stop_on_error(1);}
  sscanf(line,"%d",&dat->nsnps);
  if (dat->nsnps <= 0) { printf("Number of sites must be > 0. Exiting...\n"); stop_on_error(1);}

  dat->positions=malloc(dat->nsnps*sizeof(double));
  dat->lambda=malloc((dat->nsnps-1)*sizeof(double));
  dat->cond_chromosomes=malloc((nhaps-2)*sizeof(int *));
  dat->ind_chromosomes=malloc(2*sizeof(int *));
  for (i=0; i<(nhaps-2); i++)
    dat->cond_chromosomes[i]=malloc(dat->nsnps*sizeof(int));
  for (i=0; i<2; i++)
    dat->ind_chromosomes[i]=malloc(dat->nsnps*sizeof(int));

  /* Positions */
  fgets(line,100000000,fd);
  if (line==NULL) { printf("ReadData2: error with PHASE-style input file\n"); stop_on_error(1);}
  step=line;
  reading(&step,"%s",waste);
  for (i=0; i<dat->nsnps; i++)
    {
      reading(&step,"%lf",&dat->positions[i]);
      if (dat->positions[i] < 0) { printf("Basepair positions must be >= 0. Exiting...\n"); stop_on_error(1);}
      if(i>0) if (dat->positions[i] == dat->positions[i-1])
	{
	  if(jitter_locations==0) { 
	    printf("Basepair positions must be increasing (at SNPs %i-%i with positions %f-%f). Rerun with option \"-J\" to ignore. Exiting...\n",i-1,i,dat->positions[i-1],dat->positions[i]); stop_on_error(1);
	  }else{
	    double newloc=dat->positions[i-1]+1;
	    dat->positions[i]=newloc;
	  }
	}
      if (i < (dat->nsnps-1)) dat->lambda[i] = 1.0;
    }

  /* Position type */
  fgets(line,100000000,fd);
  if (line==NULL) { printf("ReadData2: error with PHASE-style input file\n"); stop_on_error(1);}

  for(i = 0; i < dat->nhaps_startpop; i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("ReadData2: error with PHASE-style input file\n"); stop_on_error(1);}
      for (j=0; j<dat->nsnps;j++) {
	if(line[j]=='0')
	  dat->cond_chromosomes[i][j]=0;
	if(line[j]=='1')
	  dat->cond_chromosomes[i][j]=1;
	if(line[j]=='A')
	  dat->cond_chromosomes[i][j]=2;
	if(line[j]=='C')
	  dat->cond_chromosomes[i][j]=3;
	if(line[j]=='G')
	  dat->cond_chromosomes[i][j]=4;
	if(line[j]=='T')
	  dat->cond_chromosomes[i][j]=5;
	if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T'))
	  {
	    printf("Allele-type invalid for hap%d, snp%d. Exiting...\n",i+1,j+1);
	    stop_on_error(1);
	  }
      }
    }
  for (i= 0; i < (2*ind_val); i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("ReadData2: error with PHASE-style input file\n"); stop_on_error(1);}
      for (j=0; j<dat->nsnps;j++) {
	if(line[j]=='0')
	  dat->cond_chromosomes[(dat->nhaps_startpop+i)][j]=0;
	if(line[j]=='1')
	  dat->cond_chromosomes[(dat->nhaps_startpop+i)][j]=1;
	if(line[j]=='A')
	  dat->cond_chromosomes[(dat->nhaps_startpop+i)][j]=2;
	if(line[j]=='C')
	  dat->cond_chromosomes[(dat->nhaps_startpop+i)][j]=3;
	if(line[j]=='G')
	  dat->cond_chromosomes[(dat->nhaps_startpop+i)][j]=4;
	if(line[j]=='T')
	  dat->cond_chromosomes[(dat->nhaps_startpop+i)][j]=5;
	if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T'))
	  {
	    printf("Allele-type invalid for hap%d, snp%d. Exiting...\n",dat->nhaps_startpop+i+1,j+1);
	    stop_on_error(1);
	  }
      }
    }
  for(i = 0; i < 2; i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("ReadData2: error with PHASE-style input file\n"); stop_on_error(1);}
      for (j=0; j<dat->nsnps;j++) {
	if(line[j]=='0')
	  dat->ind_chromosomes[i][j]=0;
	if(line[j]=='1')
	  dat->ind_chromosomes[i][j]=1;
	if(line[j]=='A')
	  dat->ind_chromosomes[i][j]=2;
	if(line[j]=='C')
	  dat->ind_chromosomes[i][j]=3;
	if(line[j]=='G')
	  dat->ind_chromosomes[i][j]=4;
	if(line[j]=='T')
	  dat->ind_chromosomes[i][j]=5;
	if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T'))
	  {
	    printf("Allele type invalid for hap%d, snp%d. Exiting...\n",dat->nhaps_startpop+2*ind_val+i+1,j+1);
	    stop_on_error(1);
	  }
      }
    }

  for (i= 0; i < (nhaps-2*ind_val-dat->nhaps_startpop-2); i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
      for (j=0; j<dat->nsnps;j++) {
	if(line[j]=='0')
	  dat->cond_chromosomes[(dat->nhaps_startpop+2*ind_val+i)][j]=0;
	if(line[j]=='1')
	  dat->cond_chromosomes[(dat->nhaps_startpop+2*ind_val+i)][j]=1;
	if(line[j]=='A')
	  dat->cond_chromosomes[(dat->nhaps_startpop+2*ind_val+i)][j]=2;
	if(line[j]=='C')
	  dat->cond_chromosomes[(dat->nhaps_startpop+2*ind_val+i)][j]=3;
	if(line[j]=='G')
	  dat->cond_chromosomes[(dat->nhaps_startpop+2*ind_val+i)][j]=4;
	if(line[j]=='T')
	  dat->cond_chromosomes[(dat->nhaps_startpop+2*ind_val+i)][j]=5;
	if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T'))
	  {
	    printf("Allele-type invalid for hap%d, snp%d. Exiting...\n",dat->nhaps_startpop+2*ind_val+i+1,j+1);
	    stop_on_error(1);
	  }
      }
    }

  free(line);
  return dat;

}

struct data_t *ReadDataHap(FILE *fd, int ind_val){
  struct data_t *dat;
  char * line = malloc(100000000 * sizeof(char));
  char *step;
  char waste[10];
  int i,j;
  int nhaps;

  dat=malloc(sizeof(struct data_t));
  /* if dat==NULL etc ... */

  /* Number of chromosomes from start population (to condition on)*/
  fgets(line,2047,fd);
  if (line==NULL) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
  sscanf(line,"%d",&dat->nhaps_startpop);
  if (dat->nhaps_startpop <= 0) { printf("Number of donor haplotypes must be > 0 (or >= 0 if you choose the self-copying '-c' or all-versus-all '-a' options). Exiting...\n"); stop_on_error(1);}

  /* Number of individuals */
  fgets(line,2047,fd);
  if (line==NULL) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
  sscanf(line,"%lf",&dat->nind);
  nhaps = (int) dat->nind;
  if (dat->nind < 0) { printf("Number of total individuals must be > 0. Exiting...\n"); stop_on_error(1);}
  if ((nhaps-1) < dat->nhaps_startpop) { printf("Number of total haps must be greater than or equal to number of donor haplotypes plus one extra individual.\n"); stop_on_error(1);}

  /* Number of SNPs */
  fgets(line,2047,fd);
  if (line==NULL) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
  sscanf(line,"%d",&dat->nsnps);
  if (dat->nsnps <= 0) { printf("Number of sites must be > 0. Exiting...\n"); stop_on_error(1);}

  dat->positions=malloc(dat->nsnps*sizeof(double));
  dat->lambda=malloc((dat->nsnps-1)*sizeof(double));
  dat->cond_chromosomes=malloc(dat->nhaps_startpop*sizeof(int *));
  dat->ind_chromosomes=malloc(1*sizeof(int *));
  for (i=0; i<dat->nhaps_startpop; i++)
    dat->cond_chromosomes[i]=malloc(dat->nsnps*sizeof(int));
  for (i=0; i<1; i++)
    dat->ind_chromosomes[i]=malloc(dat->nsnps*sizeof(int));

  /* Positions */
  fgets(line,100000000,fd);
  if (line==NULL) { printf("ReadDataHap: error with PHASE-style input file\n"); stop_on_error(1);}
  step=line;
  reading(&step,"%s",waste);
  for (i=0; i<dat->nsnps; i++)
    {
      reading(&step,"%lf",&dat->positions[i]);
      if (dat->positions[i] < 0) { printf("Basepair positions must be >= 0. Exiting...\n"); stop_on_error(1);}
      if(i>0) if (dat->positions[i] == dat->positions[i-1])
	{
	  if(jitter_locations==0) { 
	    printf("Basepair positions must be increasing (at SNPs %i-%i with positions %f-%f). Rerun with option \"-J\" to ignore. Exiting...\n",i-1,i,dat->positions[i-1],dat->positions[i]); stop_on_error(1);
	  }else{
	    double newloc=dat->positions[i-1]+1;
	    dat->positions[i]=newloc;
	  }
	}
      if (i < (dat->nsnps-1)) dat->lambda[i] = 1.0;
    }

  /* Position type */
  fgets(line,100000000,fd);
  if (line==NULL) { printf("error with PHASE-style input file\n"); stop_on_error(1);}

  for(i = 0; i < dat->nhaps_startpop; i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
      for (j=0; j<dat->nsnps;j++) {
	if(line[j]=='0')
	  dat->cond_chromosomes[i][j]=0;
	if(line[j]=='1')
	  dat->cond_chromosomes[i][j]=1;
	if(line[j]=='A')
	  dat->cond_chromosomes[i][j]=2;
	if(line[j]=='C')
	  dat->cond_chromosomes[i][j]=3;
	if(line[j]=='G')
	  dat->cond_chromosomes[i][j]=4;
	if(line[j]=='T')
	  dat->cond_chromosomes[i][j]=5;
	if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T'))
	  {
	    printf("Allele-type invalid for hap%d, snp%d. Exiting...\n",i+1,j+1);
	    stop_on_error(1);
	  }
      }
    }
  for (i= 0; i < ind_val; i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
   }
  for(i = 0; i < 1; i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
      for (j=0; j<dat->nsnps;j++) {
	if(line[j]=='0')
	  dat->ind_chromosomes[i][j]=0;
	if(line[j]=='1')
	  dat->ind_chromosomes[i][j]=1;
	if(line[j]=='A')
	  dat->ind_chromosomes[i][j]=2;
	if(line[j]=='C')
	  dat->ind_chromosomes[i][j]=3;
	if(line[j]=='G')
	  dat->ind_chromosomes[i][j]=4;
	if(line[j]=='T')
	  dat->ind_chromosomes[i][j]=5;
	if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T'))
	  {
	    printf("Allele type invalid for hap%d, snp%d. Exiting...\n",dat->nhaps_startpop+ind_val+i+1,j+1);
	    stop_on_error(1);
	  }
      }
    }

  free(line);
  return dat;

}

struct data_t *ReadDataHap2(FILE *fd, int ind_val){
  struct data_t *dat;
  char * line = malloc(100000000 * sizeof(char));
  char *step;
  char waste[10];
  int i,j;
  int nhaps;

  dat=malloc(sizeof(struct data_t));
  /* if dat==NULL etc ... */

  /* Number of chromosomes from start population (to condition on)*/
  fgets(line,2047,fd);
  if (line==NULL) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
  sscanf(line,"%d",&dat->nhaps_startpop);
  if (dat->nhaps_startpop < 0) { printf("Number of donor haplotypes must be >= 0. Exiting...\n"); stop_on_error(1);}

  /* Number of individuals */
  fgets(line,2047,fd);
  if (line==NULL) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
  sscanf(line,"%lf",&dat->nind);
  nhaps = (int) dat->nind;
  if (dat->nind <= 0) { printf("Number of total individuals must be > 0. Exiting...\n"); stop_on_error(1);}
  if ((nhaps-1) < dat->nhaps_startpop) { printf("Number of total haps must be greater than or equal to number of donor haplotypes plus one extra individual.\n"); stop_on_error(1);}

  /* Number of SNPs */
  fgets(line,2047,fd);
  if (line==NULL) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
  sscanf(line,"%d",&dat->nsnps);
  if (dat->nsnps <= 0) { printf("Number of sites must be > 0. Exiting...\n"); stop_on_error(1);}

  dat->positions=malloc(dat->nsnps*sizeof(double));
  dat->lambda=malloc((dat->nsnps-1)*sizeof(double));
  dat->cond_chromosomes=malloc((nhaps-1)*sizeof(int *));
  dat->ind_chromosomes=malloc(1*sizeof(int *));
  for (i=0; i<(nhaps-1); i++)
    dat->cond_chromosomes[i]=malloc(dat->nsnps*sizeof(int));
  for (i=0; i<1; i++)
    dat->ind_chromosomes[i]=malloc(dat->nsnps*sizeof(int));

  /* Positions */
  fgets(line,100000000,fd);
  if (line==NULL) { printf("ReadDataHap2: error with PHASE-style input file\n"); stop_on_error(1);}
  step=line;
  reading(&step,"%s",waste);
  for (i=0; i<dat->nsnps; i++)
    {
      reading(&step,"%lf",&dat->positions[i]);
      if (dat->positions[i] < 0) { printf("Basepair positions must be >= 0. Exiting...\n"); stop_on_error(1);}
      if(i>0) if (dat->positions[i] == dat->positions[i-1])
	{
	  if(jitter_locations==0) { 
	    printf("Basepair positions must be increasing (at SNPs %i-%i with positions %f-%f). Rerun with option \"-J\" to ignore. Exiting...\n",i-1,i,dat->positions[i-1],dat->positions[i]); stop_on_error(1);
	  }else{
	    double newloc=dat->positions[i-1]+1;
	    dat->positions[i]=newloc;
	  }
	}
      if (i < (dat->nsnps-1)) dat->lambda[i] = 1.0;
    }

  /* Position type */
  fgets(line,100000000,fd);
  if (line==NULL) { printf("error with PHASE-style input file\n"); stop_on_error(1);}

  for(i = 0; i < dat->nhaps_startpop; i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
      for (j=0; j<dat->nsnps;j++) {
	if(line[j]=='0')
	  dat->cond_chromosomes[i][j]=0;
	if(line[j]=='1')
	  dat->cond_chromosomes[i][j]=1;
	if(line[j]=='A')
	  dat->cond_chromosomes[i][j]=2;
	if(line[j]=='C')
	  dat->cond_chromosomes[i][j]=3;
	if(line[j]=='G')
	  dat->cond_chromosomes[i][j]=4;
	if(line[j]=='T')
	  dat->cond_chromosomes[i][j]=5;
	if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T'))
	  {
	    printf("Allele-type invalid for hap%d, snp%d. Exiting...\n",i+1,j+1);
	    stop_on_error(1);
	  }
      }
    }
  for (i=0; i < ind_val; i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
      for (j=0; j<dat->nsnps;j++) {
	if(line[j]=='0')
	  dat->cond_chromosomes[(dat->nhaps_startpop+i)][j]=0;
	if(line[j]=='1')
	  dat->cond_chromosomes[(dat->nhaps_startpop+i)][j]=1;
	if(line[j]=='A')
	  dat->cond_chromosomes[(dat->nhaps_startpop+i)][j]=2;
	if(line[j]=='C')
	  dat->cond_chromosomes[(dat->nhaps_startpop+i)][j]=3;
	if(line[j]=='G')
	  dat->cond_chromosomes[(dat->nhaps_startpop+i)][j]=4;
	if(line[j]=='T')
	  dat->cond_chromosomes[(dat->nhaps_startpop+i)][j]=5;
	if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T'))
	  {
	    printf("Allele-type invalid for hap%d, snp%d. Exiting...\n",dat->nhaps_startpop+i+1,j+1);
	    stop_on_error(1);
	  }
      }
    }
  for(i = 0; i < 1; i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
      for (j=0; j<dat->nsnps;j++) {
	if(line[j]=='0')
	  dat->ind_chromosomes[i][j]=0;
	if(line[j]=='1')
	  dat->ind_chromosomes[i][j]=1;
	if(line[j]=='A')
	  dat->ind_chromosomes[i][j]=2;
	if(line[j]=='C')
	  dat->ind_chromosomes[i][j]=3;
	if(line[j]=='G')
	  dat->ind_chromosomes[i][j]=4;
	if(line[j]=='T')
	  dat->ind_chromosomes[i][j]=5;
	if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T'))
	  {
	    printf("Allele type invalid for hap%d, snp%d. Exiting...\n",dat->nhaps_startpop+ind_val+i+1,j+1);
	    stop_on_error(1);
	  }
      }
    }

  for (i= 0; i < (nhaps-ind_val-dat->nhaps_startpop-1); i++)
    {
      fgets(line,100000000,fd);
      if ((line==NULL)||(strlen(line)!=(dat->nsnps+1))) { printf("error with PHASE-style input file\n"); stop_on_error(1);}
      for (j=0; j<dat->nsnps;j++) {
	if(line[j]=='0')
	  dat->cond_chromosomes[(dat->nhaps_startpop+ind_val+i)][j]=0;
	if(line[j]=='1')
	  dat->cond_chromosomes[(dat->nhaps_startpop+ind_val+i)][j]=1;
	if(line[j]=='A')
	  dat->cond_chromosomes[(dat->nhaps_startpop+ind_val+i)][j]=2;
	if(line[j]=='C')
	  dat->cond_chromosomes[(dat->nhaps_startpop+ind_val+i)][j]=3;
	if(line[j]=='G')
	  dat->cond_chromosomes[(dat->nhaps_startpop+ind_val+i)][j]=4;
	if(line[j]=='T')
	  dat->cond_chromosomes[(dat->nhaps_startpop+ind_val+i)][j]=5;
	if ((line[j]!='0')&&(line[j]!='1')&&(line[j]!='A')&&(line[j]!='G')&&(line[j]!='C')&&(line[j]!='T'))
	  {
	    printf("Allele-type invalid for hap%d, snp%d. Exiting...\n",dat->nhaps_startpop+ind_val+i+1,j+1);
	    stop_on_error(1);
	  }
      }
    }

  free(line);
  return dat;

}

void DestroyData(struct data_t *dat)
{
  int i;
  for (i=0; i<dat->nhaps_startpop; i++)
    free(dat->cond_chromosomes[i]);
  for (i=0; i<2; i++)
    free(dat->ind_chromosomes[i]);
  free(dat->cond_chromosomes);
  free(dat->ind_chromosomes);
  free(dat->positions);
  free(dat->lambda);
}

void DestroyData2(struct data_t *dat)
{
  int i, nhaps;
  nhaps = (int) 2*dat->nind;
  for (i=0; i<(nhaps-2); i++)
    free(dat->cond_chromosomes[i]);
  for (i=0; i<2; i++)
    free(dat->ind_chromosomes[i]);
  free(dat->cond_chromosomes);
  free(dat->ind_chromosomes);
  free(dat->positions);
  free(dat->lambda);
}

void DestroyDataHap(struct data_t *dat)
{
  int i;
  for (i=0; i<dat->nhaps_startpop; i++)
    free(dat->cond_chromosomes[i]);
  for (i=0; i<1; i++)
    free(dat->ind_chromosomes[i]);
  free(dat->ind_chromosomes);
  free(dat->cond_chromosomes);
  free(dat->positions);
  free(dat->lambda);
}

void DestroyDataHap2(struct data_t *dat)
{
  int i, nhaps;
  nhaps = (int) dat->nind;
  for (i=0; i<(nhaps-1); i++)
    free(dat->cond_chromosomes[i]);
  for (i=0; i<1; i++)
    free(dat->ind_chromosomes[i]);
  free(dat->cond_chromosomes);
  free(dat->ind_chromosomes);
  free(dat->positions);
  free(dat->lambda);
}

double standard_normal_samp()
{
      /* SAMPLE Z1,Z2 FROM A NORMAL(0,1) USING X1,X2 ~ UNIF(0,1) */
               // (THIS IS THE BOX-MULLER TRANSFORMATION)
  double z1, z2;
  double x1,x2;

  x1 = (double)rand()/RAND_MAX;
  x2 = (double)rand()/RAND_MAX;

  z1 = pow((-2.0 * log(x1)),0.5) * cos(2 * PI * x2);
  z2 = pow((-2.0 * log(x1)),0.5) * sin(2 * PI * x2);

  return(z1);
}

double ** sampler(int * newh, int ** existing_h, int *p_Nloci, int *p_Nhaps, int *p_nchr,  double p_rhobar, double * MutProb_vec, int * allelic_type_count_vec, double * lambda, double * pos, int nsampTOT, double * copy_prob, double * copy_probSTART, int * pop_vec, int ndonorpops, double region_size, int run_num, int run_thres, int all_versus_all_ind, int haploid_ind, int unlinked_ind, int ind_val, int print_file9_ind, FILE *fout, FILE *fout3, FILE *fout9)
{

  double rounding_val = 1.0/10000000.0;  // for regional_counts; c is a bit lame

  int i, j, locus;
  double sum, Theta;
  double prob, total_prob, total_gen_dist;
  double total_prob_from_i_to_i,total_prob_to_i_exclude_i,total_prob_from_i_exclude_i,constant_exclude_i,constant_from_i_to_i,constant_exclude_i_both_sides,total_prob_from_any_to_any_exclude_i;
  double total_regional_chunk_count, total_ind_sum;
  int num_regions;
  double * TransProb = malloc( ((*p_Nloci)-1) * sizeof(double));
  double N_e_new;
  int * sample_state = malloc(*p_Nloci * sizeof(int));
  double ObsStateProb, ObsStateProbPREV;
  double delta;
             //correction to PAC-A rho_est
  double random_unif, random_unifSWITCH;
  double no_switch_prob;
  double ** Alphamat = malloc(*p_Nhaps * sizeof(double *));
  double * BetavecPREV = malloc(*p_Nhaps * sizeof(double));
  double * BetavecCURRENT = malloc(*p_Nhaps * sizeof(double));
  double Alphasum, Alphasumnew, Betasum, Betasumnew;
  double large_num;
  double * copy_prob_new = malloc(*p_Nhaps * sizeof(double));
  double * copy_prob_newSTART = malloc(*p_Nhaps * sizeof(double));
  double * Alphasumvec = malloc(*p_Nloci * sizeof(double));
  double * expected_transition_prob = malloc((*p_Nloci-1)*sizeof(double));
  double * corrected_chunk_count = malloc(*p_Nhaps * sizeof(double));
  double * regional_chunk_count = malloc(*p_Nhaps * sizeof(double));
  double * expected_chunk_length = malloc(*p_Nhaps * sizeof(double));
  double * expected_differences = malloc(*p_Nhaps * sizeof(double));
  double * regional_chunk_count_sum = malloc(ndonorpops * sizeof(double));
  double * regional_chunk_count_sum_final = malloc(ndonorpops * sizeof(double));
  double * regional_chunk_count_sum_squared_final = malloc(ndonorpops * sizeof(double));
  double * ind_snp_sum_vec = malloc(ndonorpops * sizeof(double));
  double * snp_info_measure=malloc(ndonorpops * sizeof(double));
  double * exp_copy_pop=malloc(ndonorpops * sizeof(double));
  double expected_chunk_length_sum, sum_prob;

  double ** copy_prob_new_mat = malloc(8 * sizeof(double *));
  for (i=0; i < 8; i++)
    copy_prob_new_mat[i] = malloc((*p_Nhaps+1) * sizeof(double));

  for(i=0 ; i< *p_Nhaps ; i++)
    {
      Alphamat[i] = malloc(*p_Nloci * sizeof(double));
    }
  for (i=0; i < ndonorpops; i++)
    {
      regional_chunk_count_sum[i] = 0.0;
      regional_chunk_count_sum_final[i] = 0.0;
      regional_chunk_count_sum_squared_final[i] = 0.0;
      ind_snp_sum_vec[i]=0.0;
      snp_info_measure[i]=0.0;
    }


				// Theta as given in Li and Stephens
  sum = 0;
  for(i = 1; i < *p_nchr; i++){
    sum = sum + 1.0/i;
  }
  Theta = 1.0 / sum;

  for (i=0; i < *p_Nhaps; i++)
    {
      if (MutProb_vec[i]<0) MutProb_vec[i]=0.5 * Theta/(*p_Nhaps + Theta);
      //if (MutProb_vec[i]<rounding_val) MutProb_vec[i]=rounding_val;
    }

    // TransProb[i] is probability of copying mechanism "jumping" between
  //   loci i and i+1
  delta = 1.0;

  if (unlinked_ind==0 && lambda[0] >= 0) TransProb[0] = 1 - exp(-1 * (pos[1]-pos[0]) * delta * p_rhobar*lambda[0]);
  if (unlinked_ind==1 || lambda[0] < 0) TransProb[0] = 1.0;
  /*
  if (TransProb[0]<rounding_val)
	{
	  printf("Transition prob is too low; will likely cause rounding errors. Exiting...\n");
	  stop_on_error(1);
	}
  */

  for(locus = 1; locus < *p_Nloci - 1; locus++)
    {
      delta = 1.0;
      if (unlinked_ind==0 && lambda[locus] >= 0) TransProb[locus] = 1 - exp(-1 * (pos[locus+1]-pos[locus]) * delta * p_rhobar*lambda[locus]);
      if (unlinked_ind==1 || lambda[locus] < 0) TransProb[locus] = 1.0;
      /*
      if (TransProb[locus]<rounding_val)
	{
	  printf("Transition prob is too low; will likely cause rounding errors. Exiting...\n");
	  stop_on_error(1);
	}
      */
    }

      /* FORWARDS ALGORITHM: (Rabiner 1989, p.262) */
      /* INITIALIZATION: */
  Alphasum = 0.0;
  for (i=0; i < *p_Nhaps; i++)
    {
      ObsStateProb = (1-MutProb_vec[i]) * (newh[0] == existing_h[i][0]) + MutProb_vec[i] * (newh[0] != existing_h[i][0]);
      /*
      if (ObsStateProb<rounding_val)
	{
	  printf("Mutation (emission) rate is too low; will likely cause rounding errors. Exiting...\n");
	  stop_on_error(1);
	}
      */

      Alphamat[i][0] = log(copy_probSTART[i]*ObsStateProb);
      Alphasum = Alphasum + exp(Alphamat[i][0])*TransProb[0];
      //if (copy_probSTART[i] < 0) printf("%d %lf %lf\n",i,Alphamat[i][0],10000000000000*copy_probSTART[i]);
    }

       /* INDUCTION: */
  Alphasum = log(Alphasum);
  for (locus=1; locus < *p_Nloci; locus++)
    {
      Alphasumnew = 0.0;
      large_num = -1.0*Alphasum;
      for (i=0; i < *p_Nhaps; i++)
	{
	  ObsStateProb = (1-MutProb_vec[i]) * (newh[locus] == existing_h[i][locus]) + MutProb_vec[i] * (newh[locus] != existing_h[i][locus]);
	  /*
	  if (ObsStateProb<rounding_val)
	    {
	      printf("Mutation (emission) rate is too low; will likely cause rounding errors. Exiting...\n");
	      stop_on_error(1);
	    }
	  */

	  Alphamat[i][locus] = log(ObsStateProb*copy_prob[i]*exp(Alphasum+large_num) + ObsStateProb*(1-TransProb[(locus-1)])*exp(Alphamat[i][(locus-1)]+large_num)) - large_num;
	  if (locus < (*p_Nloci - 1)) Alphasumnew = Alphasumnew + exp(Alphamat[i][locus]+large_num)*TransProb[locus];
	  if (locus == (*p_Nloci - 1)) Alphasumnew = Alphasumnew + exp(Alphamat[i][locus]+large_num);
	}
      Alphasum = log(Alphasumnew)-large_num;
    }
  //if (Alphasum == (Alphasum*5))
  if (isnan(Alphasum))
    {
      printf("Negative or NaN likelihood. Could be because emission or transition probabilities are too low??...Exiting...\n");
      //for (i=0; i < *p_Nhaps; i++) printf("%d %lf %lf\n",i,copy_prob[i],log(copy_prob[i]));
      stop_on_error(1);
    }
  fprintf(fout3," %.10lf",Alphasum);

  for (i = 0; i < *p_Nhaps; i++)
    {
      copy_prob_new[i] = 0.0;
      corrected_chunk_count[i] = 0.0;
      expected_chunk_length[i] = 0.0;
      expected_differences[i] = 0.0;
      regional_chunk_count[i] = 0.0;
    }
  total_regional_chunk_count=0.0;
  num_regions=0;
  if (run_num <= (run_thres-1))
    {
          /* BACKWARDS ALGORITHM: (Rabiner 1989, p.263) */
          /* INITIALIZATION: */
      Betasum = 0.0;
      if (run_num == (run_thres-1))
	{
	  for (i=0; i < ndonorpops; i++)
	    exp_copy_pop[i]=0.0;
	}
      for(i=0; i < *p_Nhaps; i++)
	{
	  ObsStateProb = (1-MutProb_vec[i]) * (newh[(*p_Nloci-1)] == existing_h[i][(*p_Nloci-1)]) + MutProb_vec[i] * (newh[(*p_Nloci-1)] != existing_h[i][(*p_Nloci-1)]);
	  /*
	  if (ObsStateProb<rounding_val)
	    {
	      printf("Mutation (emission) rate is too low; will likely cause rounding errors. Exiting...\n");
	      stop_on_error(1);
	    }
	  */
	  BetavecPREV[i] = 0.0;
	  Betasum = Betasum + TransProb[(*p_Nloci-2)]*copy_prob[i]*ObsStateProb*exp(BetavecPREV[i]);
	  if (run_num == (run_thres-1)) exp_copy_pop[pop_vec[i]]=exp_copy_pop[pop_vec[i]]+exp(BetavecPREV[i]+Alphamat[i][(*p_Nloci-1)]-Alphasum);
	                    // for estimating new mutation rates:
	  expected_differences[i]=expected_differences[i]+exp(Alphamat[i][(*p_Nloci-1)]-Alphasum)*(newh[(*p_Nloci-1)] != existing_h[i][(*p_Nloci-1)]);
	}
      if ((run_num == (run_thres-1)) && (print_file9_ind==1))
 	{
	  gzprintf(fout9,"%.0lf",pos[(*p_Nloci-1)]);
	  if (all_versus_all_ind==0)
	    {
	      for (i=0; i < ndonorpops; i++)
		gzprintf(fout9," %lf",exp_copy_pop[i]);
	    }
	  if (all_versus_all_ind==1)
	    {
	      for (i=0; i < ndonorpops; i++)
		{
		  if (i == ind_val) gzprintf(fout9," 0.0");
		  gzprintf(fout9," %lf",exp_copy_pop[i]);
		}
	    }
	  gzprintf(fout9,"\n");
	}

           /* INDUCTION: */
      Betasum = log(Betasum);
          /* CALCULATE EXPECTED NUMBER OF TIMES OF COPYING TO EACH DONOR POP (Rabiner 1989, p.263,265 or Scheet/Stephens 2006 Appendix C): */
      for (locus = (*p_Nloci-2); locus >= 0; locus--)
	{
	  Betasumnew = 0.0;
	  large_num = -1.0*Betasum;
	  total_prob=0.0;

	  constant_exclude_i = 0.5;
	  constant_from_i_to_i = 1.0;
	  constant_exclude_i_both_sides = 0.0;
	  expected_chunk_length_sum=0.0;
	  sum_prob=0.0;
	  if (run_num == (run_thres-1))
	    {
	      for (i=0; i < ndonorpops; i++)
		exp_copy_pop[i]=0.0;
	    }
 	  for (i = 0; i < *p_Nhaps; i++)
	    {
	      ObsStateProb = (1-MutProb_vec[i]) * (newh[locus] == existing_h[i][locus]) + MutProb_vec[i] * (newh[locus] != existing_h[i][locus]);
	      ObsStateProbPREV = (1-MutProb_vec[i]) * (newh[(locus+1)] == existing_h[i][(locus+1)]) + MutProb_vec[i] * (newh[(locus+1)] != existing_h[i][(locus+1)]);
	      /*
	      if ((ObsStateProb<rounding_val) || (ObsStateProbPREV<rounding_val))
		{
		  printf("Mutation (emission) rate is too low; will likely cause rounding errors. Exiting...\n");
		  stop_on_error(1);
		}
	      */
	      BetavecCURRENT[i] = log(exp(Betasum+large_num) + (1-TransProb[locus]) * ObsStateProbPREV*exp(BetavecPREV[i] + large_num)) - large_num;
	      if (locus > 0) Betasumnew = Betasumnew + TransProb[(locus-1)]*copy_prob[i]*ObsStateProb*exp(BetavecCURRENT[i] + large_num);
	      if (locus == 0) copy_prob_newSTART[i] = exp(Alphamat[i][0] + BetavecCURRENT[i] - Alphasum);
	      total_prob = total_prob + exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]);

	      copy_prob_new[i] = copy_prob_new[i] + exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]);

	      total_prob_from_i_to_i = exp(Alphamat[i][locus]+BetavecPREV[i]-Alphasum)*ObsStateProbPREV*(1-TransProb[locus]+TransProb[locus]*copy_prob[i]);
	      total_prob_to_i_exclude_i = exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]+TransProb[locus]*copy_prob[i]);
	      total_prob_from_i_exclude_i = exp(Alphamat[i][locus]+BetavecCURRENT[i]-Alphasum) - exp(Alphamat[i][locus]+BetavecPREV[i]-Alphasum)*ObsStateProbPREV*(1-TransProb[locus]+TransProb[locus]*copy_prob[i]);
	      total_prob_from_any_to_any_exclude_i = 1.0-exp(Alphamat[i][locus]+BetavecCURRENT[i]-Alphasum)-exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)+exp(Alphamat[i][locus]+BetavecPREV[i]-Alphasum)*ObsStateProbPREV*(1-TransProb[locus]+TransProb[locus]*copy_prob[i]);

	      regional_chunk_count[i]=regional_chunk_count[i]+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]));
	      total_regional_chunk_count=total_regional_chunk_count+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]));
	      ind_snp_sum_vec[pop_vec[i]]=ind_snp_sum_vec[pop_vec[i]]+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]));

	      //corrected_chunk_count[i]=corrected_chunk_count[i]+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]))*(1.0+(1.0/(*p_Nhaps))*((p_rhobar*(pos[(locus+1)]-pos[locus])*delta*lambda[locus]/(*p_Nhaps))/(1.0-exp(-1.0*p_rhobar*(pos[(locus+1)]-pos[locus])*delta*lambda[locus]/(*p_Nhaps)))-1.0));
	      corrected_chunk_count[i]=corrected_chunk_count[i]+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]));
	      if (unlinked_ind==0 && lambda[locus]>=0) expected_chunk_length[i]=expected_chunk_length[i]+100*(pos[locus+1]-pos[locus])*delta*lambda[locus]*(constant_from_i_to_i*total_prob_from_i_to_i+constant_exclude_i*(total_prob_to_i_exclude_i+total_prob_from_i_exclude_i)+constant_exclude_i_both_sides*total_prob_from_any_to_any_exclude_i);  // multiply by 100 to get cM
	      expected_chunk_length_sum=expected_chunk_length_sum+constant_from_i_to_i*total_prob_from_i_to_i+constant_exclude_i*(total_prob_to_i_exclude_i+total_prob_from_i_exclude_i)+constant_exclude_i_both_sides*total_prob_from_any_to_any_exclude_i;
	                    // for estimating new mutation rates:
	      expected_differences[i]=expected_differences[i]+exp(Alphamat[i][locus]+BetavecCURRENT[i]-Alphasum)*(newh[locus] != existing_h[i][locus]);
	      BetavecPREV[i] = BetavecCURRENT[i];

	      if (run_num == (run_thres-1)) exp_copy_pop[pop_vec[i]]=exp_copy_pop[pop_vec[i]]+exp(BetavecCURRENT[i]+Alphamat[i][locus]-Alphasum);

	      sum_prob=sum_prob+total_prob_from_i_to_i+total_prob_to_i_exclude_i+total_prob_from_i_exclude_i;
	    }

	  if ((run_num == (run_thres-1)) && (print_file9_ind==1))
	    {
	      gzprintf(fout9,"%.0lf",pos[locus]);
	      if (all_versus_all_ind==0)
		{
		  for (i=0; i < ndonorpops; i++)
		    gzprintf(fout9," %lf",exp_copy_pop[i]);
		}
	      if (all_versus_all_ind==1)
		{
		  for (i=0; i < ndonorpops; i++)
		    {
		      if (i == ind_val) gzprintf(fout9," 0.0");
		      gzprintf(fout9," %lf",exp_copy_pop[i]);
		    }
		}
	      gzprintf(fout9,"\n");
	    }


	  expected_transition_prob[locus]=total_prob;
	  if (locus > 0) Betasum = log(Betasumnew) - large_num;

	  if ((total_regional_chunk_count+rounding_val) >= region_size)
	    {
	      for (i = 0; i < *p_Nhaps; i++)
		{
		  //printf("%d %d %lf %lf\n",i,pop_vec[i],regional_chunk_count[i],regional_chunk_count_sum[pop_vec[i]]);
		  regional_chunk_count_sum[pop_vec[i]]=regional_chunk_count_sum[pop_vec[i]]+regional_chunk_count[i];
		  regional_chunk_count[i]=0.0;
		}
	      for (i = 0; i < ndonorpops; i++)
		{
		  regional_chunk_count_sum_final[i]=regional_chunk_count_sum_final[i]+regional_chunk_count_sum[i];
		  regional_chunk_count_sum_squared_final[i]=regional_chunk_count_sum_squared_final[i]+pow(regional_chunk_count_sum[i],2.0);
		  regional_chunk_count_sum[i]=0.0;
		}
	      total_regional_chunk_count=0.0;
	      num_regions=num_regions+1;
	    }
	  total_ind_sum=0.0;
	  for (i = 0; i < ndonorpops; i++)
	    total_ind_sum=total_ind_sum+ind_snp_sum_vec[i];
	  for (i = 0; i < ndonorpops; i++)
	    {
	      snp_info_measure[i]=snp_info_measure[i]+pow((ind_snp_sum_vec[i]/total_ind_sum),2.0);
	      ind_snp_sum_vec[i]=0.0;
	    }
	}
  for (i=0; i < ndonorpops; i++)
    snp_info_measure[i]=snp_info_measure[i]/(*p_Nloci);

          /* CALCULATE EXPECTED NUMBER OF TOTAL TRANSITIONS, IN ORDER TO ESTIMATE N_e (Scheet/Stephens 2006 Appendix C (C3)): */
      total_prob=0.0;
      total_gen_dist=0.0;
      for (locus = 0; locus < (*p_Nloci-1); locus++)
	{
	  if (unlinked_ind==0 && lambda[locus] >= 0) total_gen_dist=total_gen_dist+(pos[(locus+1)]-pos[locus])*delta*lambda[locus];
	  if (unlinked_ind==0 && lambda[locus] >= 0) total_prob=total_prob+((p_rhobar*(pos[(locus+1)]-pos[locus])*delta*lambda[locus])/(1.0-exp(-1.0*p_rhobar*(pos[(locus+1)]-pos[locus])*delta*lambda[locus])))*expected_transition_prob[locus];
	}
      if (unlinked_ind==0) N_e_new = total_prob/total_gen_dist;
      if (unlinked_ind==1) N_e_new = 0.0;

              /* CALCULATE SOMETHING ANALAGOUS TO EXPECTED NUMBER OF TIMES EACH HAP i IS VISITED, CONDITIONAL ON THE OBSERVED DATA (I.E  (27) AND PARAGRAPH UNDER (38) IN RABINER 1989, Proceedings of the IEEE 77(2):257-286), BUT -- AS WE'RE ONLY COUNTING CHUNKS -- SUBTRACT OUT TIMES YOU DO NOT SWITCH */
      for (i=0; i < *p_Nhaps; i++) corrected_chunk_count[i]=corrected_chunk_count[i]+copy_prob_newSTART[i];
    }

       /* print-out samples if we've done enough iterations: */
   if (run_num == (run_thres-1))
     {

       N_e_new = p_rhobar;

       for (i=0; i < *p_Nhaps; i++)
	 {
	   copy_prob_new[i] = copy_prob[i];
	   copy_prob_newSTART[i] = copy_probSTART[i];
	 }

           /* calculate Alphasums (for efficient sampling): */
       for (locus=0; locus < *p_Nloci; locus++)
	 {
	   Alphasumvec[locus] = 0.0;
	   large_num = Alphamat[0][locus];
	   for (i = 1; i < *p_Nhaps; i++)
	     {
	       if (Alphamat[i][locus] > large_num)
		 large_num = Alphamat[i][locus];
	     }
	   large_num = -1.0*large_num;
	   for (i = 0; i < *p_Nhaps; i++)
	     Alphasumvec[locus] = Alphasumvec[locus] + exp(Alphamat[i][locus]+large_num);
	   Alphasumvec[locus] = log(Alphasumvec[locus]) - large_num;
	 }

              /* SAMPLING ALGORITHM: (from Falush, Stephens, & Pritchard (2003) Genetics 164:1567-1587) */
       for (j = 0; j < nsampTOT; j++)
	 {
	   //printf("sample %d\n",j);
	      /* sample last position: */
	   total_prob = 0.0;
	   large_num = Alphamat[0][(*p_Nloci-1)];
	   for (i = 1; i < *p_Nhaps; i++)
	     {
	       if (Alphamat[i][(*p_Nloci-1)] > large_num)
		 large_num = Alphamat[i][(*p_Nloci-1)];
	     }
	   large_num = -1.0*large_num;
	   random_unif = (double) rand()/RAND_MAX;
	   total_prob = Alphasumvec[(*p_Nloci-1)];
	   prob = 0.0;
	   for (i = 0; i < *p_Nhaps; i++)
	     {
	       prob = prob + exp(Alphamat[i][(*p_Nloci-1)]+large_num);
	       if (random_unif <= exp(log(prob)-large_num-total_prob))
		 {
		   sample_state[(*p_Nloci-1)] = i;
		   break;
		 }
	     }

              /* sample remaining positions: */
	   for (locus = (*p_Nloci-2); locus >= 0; locus--)
	     {
	        // first sample prob you switch and see if you need to
                   // if you do need to switch, you need to go through the below loop to figure out where to switch to
               large_num = -1.0 * Alphasumvec[locus];
	       total_prob = log(exp(Alphasumvec[locus]+large_num)*TransProb[locus]*copy_prob[sample_state[(locus+1)]] + exp(Alphamat[sample_state[(locus+1)]][locus]+large_num)*(1.0-TransProb[locus]))-large_num;
	       no_switch_prob = exp(log(exp(Alphamat[sample_state[(locus+1)]][locus]+large_num)*(1.0-TransProb[locus])) - large_num - total_prob);
               random_unifSWITCH = (double) rand()/RAND_MAX;
	       if (random_unifSWITCH <= no_switch_prob) sample_state[locus] = sample_state[(locus+1)];

	       //if (j ==0 && locus > 9500) printf("%d %d %lf %lf %lf %lf %lf\n",locus,sample_state[(locus+1)],no_switch_prob,Alphamat[sample_state[(locus+1)]][locus],large_num,total_prob,1.0-TransProb[locus]);
               if (random_unifSWITCH > no_switch_prob)
		 {
		   total_prob = 0.0;
		   large_num = Alphamat[0][locus];
		   for (i = 1; i < *p_Nhaps; i++)
		     {
		       if (Alphamat[i][locus] > large_num)
			 large_num = Alphamat[i][locus];
		     }
		   large_num = -1.0*large_num;

		   random_unif = (double) rand()/RAND_MAX;
		   total_prob = log(exp(Alphasumvec[locus]+large_num)*TransProb[locus]*copy_prob[sample_state[(locus+1)]]) - large_num;
		   prob = 0.0;
		   for (i = 0; i < *p_Nhaps; i++)
		     {
		       prob = prob + exp(Alphamat[i][locus]+large_num)*TransProb[locus]*copy_prob[sample_state[(locus+1)]];
		       if (random_unif <= exp(log(prob)-large_num-total_prob))
			 {
			   sample_state[locus] = i;
			   break;
			 }
		     }
		 }
	     }

	   fprintf(fout,"%d",j+1);
	   for (i = 0; i < *p_Nloci; i++)
	     {
	       if (all_versus_all_ind==0) fprintf(fout," %d",sample_state[i]+1);
	       if (all_versus_all_ind==1)
		 {
		   if (sample_state[i] >= ((2-haploid_ind)*ind_val))
		     fprintf(fout," %d",sample_state[i]+2-haploid_ind+1);
		   if (sample_state[i] < ((2-haploid_ind)*ind_val))
		     fprintf(fout," %d",sample_state[i]+1);
		 }
	     }
	   fprintf(fout,"\n");
	 }
     }

  for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[0][i] = copy_prob_new[i];
   for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[1][i] = copy_prob_newSTART[i];
   for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[2][i] = corrected_chunk_count[i];
   for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[3][i] = expected_chunk_length[i];
   for (i = 0; i < *p_Nhaps; i++)
     copy_prob_new_mat[4][i] = expected_differences[i];
   for (i = 0; i < ndonorpops; i++)
     copy_prob_new_mat[5][i] = regional_chunk_count_sum_final[i];
   for (i = 0; i < ndonorpops; i++)
     copy_prob_new_mat[6][i] = regional_chunk_count_sum_squared_final[i];
   for (i = 0; i < ndonorpops; i++)
     copy_prob_new_mat[7][i] = snp_info_measure[i];
   copy_prob_new_mat[0][(*p_Nhaps)] = N_e_new;
   copy_prob_new_mat[1][(*p_Nhaps)] = num_regions;

   for (i=0; i < *p_Nhaps; i++)
     {
       free(Alphamat[i]);
     }
   free(Alphamat);
   free(BetavecPREV);
   free(BetavecCURRENT);
   free(TransProb);
   free(sample_state);
   free(expected_transition_prob);
   free(Alphasumvec);
   free(copy_prob_new);
   free(copy_prob_newSTART);
   free(corrected_chunk_count);
   free(regional_chunk_count);
   free(expected_chunk_length);
   free(expected_differences);
   free(regional_chunk_count_sum);
   free(regional_chunk_count_sum_final);
   free(regional_chunk_count_sum_squared_final);
   free(ind_snp_sum_vec);
   free(snp_info_measure);
   free(exp_copy_pop);

   return(copy_prob_new_mat);
}

int loglik(int nhaps_startpop, int *p_nloci, int p_nhaps, double N_e_start, double * recom_map, double * MutProb_vec, int nsampTOT, int ndonorpops, int * ndonorhaps, double * copy_prob, double * copy_probSTART, int * pop_vec, double region_size, int EMruns, int estimate_copyprob_ind, int estimate_recom_ind, int estimate_mutation_ind, int estimate_mutationALL_ind, int recipient_cond_ind, int all_versus_all_ind, int start_val, int end_val, char *filename, char *filenameDL, int donorlist_ind, int haploid_ind, int unlinked_ind, int print_file9_ind, int indcount_suppress_ind, FILE *fout, FILE *fout2, FILE *fout3, FILE *fout4, FILE *fout5, FILE *fout6, FILE *fout7, FILE *fout8, FILE *fout9)
{
  double small_copy_val=0.000000000000001; // (!!!) copy props per hap not allowed to go below this value, even if E-M wants to make them lower (!!!)
  int nhaps_condpop, nind_condpop;

  nhaps_condpop = p_nhaps-nhaps_startpop;
  if (((recipient_cond_ind==1) || (all_versus_all_ind==1)) && (haploid_ind==0)) nhaps_startpop=nhaps_startpop+nhaps_condpop-2;
  if (((recipient_cond_ind==1) || (all_versus_all_ind==1)) && (haploid_ind==1)) nhaps_startpop=nhaps_startpop+nhaps_condpop-1;

  FILE *fd, *fd3;
  struct data_t *Data;
  char *step;
  char line[2047];
  char waste[2047];
  int i, j, m, n, count, indcount, r, h;
  int nhaps, num_regions_tot;
  int ndonorval, totalhaps;
  double sum_total_diff;
  double * total_back_prob = malloc((ndonorpops+recipient_cond_ind) * sizeof(double));
  double * total_back_probSTART = malloc((ndonorpops+recipient_cond_ind) * sizeof(double));
  double * total_counts = malloc((ndonorpops+recipient_cond_ind) * sizeof(double));
  double * total_lengths = malloc((ndonorpops+recipient_cond_ind) * sizeof(double));
  double * total_differences = malloc((ndonorpops+recipient_cond_ind) * sizeof(double));
  double * total_region_counts = malloc((ndonorpops+recipient_cond_ind) * sizeof(double));
  double * total_squared_region_counts = malloc((ndonorpops+recipient_cond_ind) * sizeof(double));
  double * snp_info_measure_final = malloc((ndonorpops+recipient_cond_ind) * sizeof(double));
  double N_e_new, N_e;
  double total_prob, total_probSTART;
  double * copy_prob_new = malloc(nhaps_startpop * sizeof(double));
  double * copy_prob_newSTART = malloc(nhaps_startpop * sizeof(double));
  double * MutProb_vec_new = malloc(nhaps_startpop * sizeof(double));
  double ** back_prob = malloc(8 * sizeof(double *));
  double ** copy_prob_pop = malloc(2 * sizeof(double *));
  int * newhap = malloc((*p_nloci) * sizeof(int));
  int * ndonorhaps_vec=malloc((ndonorpops+recipient_cond_ind)*sizeof(int));

  for (i=0; i < 8; i++)
    back_prob[i] = malloc((nhaps_startpop+1) * sizeof(double));
  for (i=0; i < 2; i++)
    copy_prob_pop[i] = malloc((ndonorpops+recipient_cond_ind) * sizeof(double));

  for (i=0; i < (ndonorpops+recipient_cond_ind); i++)
    {
      total_back_prob[i] = 0.0;
      total_back_probSTART[i] = 0.0;
    }
  for (i=0; i < ndonorpops; i++)
    ndonorhaps_vec[i]=ndonorhaps[i];
  if (recipient_cond_ind==1 && haploid_ind==0) ndonorhaps_vec[ndonorpops]=nhaps_condpop-2;
  if (recipient_cond_ind==1 && haploid_ind==1) ndonorhaps_vec[ndonorpops]=nhaps_condpop-1;

  if (haploid_ind==0) nind_condpop=nhaps_condpop/2;
  if (haploid_ind==1) nind_condpop=nhaps_condpop;

  if ((all_versus_all_ind==1) && (donorlist_ind==1))
    {
      fd3 = fopen(filenameDL,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filenameDL); stop_on_error(1);}
      totalhaps=0;
      while(totalhaps < ((2-haploid_ind)*start_val+1))
	{
	  fgets(line,2047,fd3);
	  step=line;
	  reading(&step,"%s",waste);
	  reading(&step,"%d",&ndonorval);
	  totalhaps=totalhaps+ndonorval;
	}
    }

  total_prob = 0.0;
  total_probSTART = 0.0;
  indcount=((2-haploid_ind)*start_val-totalhaps+ndonorval)/(2-haploid_ind)+1;
  for (m = start_val; m < end_val; m ++)
    {
      fprintf(fout3,"IND %d\n",m+1);

      for (i=0; i < nhaps_startpop; i++)
	{
	  copy_prob_new[i] = copy_prob[i];
	  copy_prob_newSTART[i] = copy_probSTART[i];
	  MutProb_vec_new[i] = MutProb_vec[i];
	}

      fd = fopen(filename,"r");
      if (fd == NULL) { printf("error opening %s\n",filename); stop_on_error(1);}
      if (haploid_ind==0)
	{
	  if ((recipient_cond_ind==0) && (all_versus_all_ind==0)) Data = ReadData(fd,m);
	  if ((recipient_cond_ind==1) || (all_versus_all_ind==1)) Data = ReadData2(fd,m);
	}
      if (haploid_ind==1)
	{
	  if ((recipient_cond_ind==0) && (all_versus_all_ind==0)) Data = ReadDataHap(fd,m);
	  if ((recipient_cond_ind==1) || (all_versus_all_ind==1)) Data = ReadDataHap2(fd,m);
	}

                  // find number of alleles per snp (this is NOT every used, but perhaps should be to get default mutation rate correct):
      int * allelic_type_count_vec = malloc(Data->nsnps*sizeof(int));
      int * found_vec = malloc(6*sizeof(int));
      for (j=0; j < Data->nsnps; j++)
	{
	  allelic_type_count_vec[j] = 0;
	  for (i=0; i < 6; i++)
	    found_vec[i]=0;
	  for (i=0; i < nhaps_startpop; i++)
	    {
	      if ((Data->cond_chromosomes[i][j] == 0) && (found_vec[0] == 0))
		{
		  allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		  found_vec[0] = 1;
		}
	      if ((Data->cond_chromosomes[i][j] == 1) && (found_vec[1] == 0))
		{
		  allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		  found_vec[1] = 1;
		}
	      if ((Data->cond_chromosomes[i][j] == 2) && (found_vec[2] == 0))
		{
		  allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		  found_vec[2] = 1;
		}
	      if ((Data->cond_chromosomes[i][j] == 3) && (found_vec[3] == 0))
		{
		  allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		  found_vec[3] = 1;
		}
	      if ((Data->cond_chromosomes[i][j] == 4) && (found_vec[4] == 0))
		{
		  allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		  found_vec[4] = 1;
		}
	      if ((Data->cond_chromosomes[i][j] == 5) && (found_vec[5] == 0))
		{
		  allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
		  found_vec[5] = 1;
		}
	    }
	}

      for (r=0; r < EMruns; r++)
	{
	  fprintf(fout3,"%d",r);

	  total_prob = 0.0;
	  total_probSTART = 0.0;
	  for (i=0; i < (ndonorpops+recipient_cond_ind); i++)
	    {
	      total_back_prob[i] = 0.0;
	      total_back_probSTART[i] = 0.0;
	      total_counts[i]=0.0;
	      total_lengths[i]=0.0;
	      total_differences[i]=0.0;
	      total_region_counts[i]=0.0;
	      total_squared_region_counts[i]=0.0;
	      snp_info_measure_final[i]=0.0;
	    }

	  N_e_new=0.0;
	  num_regions_tot=0;
	  for(h=0; h < (2-haploid_ind); h++)
	    {
	      if (r==0) N_e=N_e_start;

	      for (n = 0; n < Data->nsnps; n ++)
		newhap[n] = Data->ind_chromosomes[h][n];

	      nhaps = nhaps_startpop;

	      if ((r == (EMruns-1)) && (haploid_ind==0))
		{
		  fprintf(fout,"HAP %d\n",2*m+h+1);
		  if (print_file9_ind==1) gzprintf(fout9,"HAP %d\n",2*m+h+1);
		}
	      if ((r == (EMruns-1)) && (haploid_ind==1))
		{
		  fprintf(fout,"HAP %d\n",m+1);
		  if (print_file9_ind==1) gzprintf(fout9,"HAP %d\n",m+1);
		}
                  /* SAMPLE FROM PAC CONDITIONAL ON COPY-PROBS: */
	      back_prob = sampler(newhap, Data->cond_chromosomes, p_nloci, &nhaps, &p_nhaps, N_e, MutProb_vec_new, allelic_type_count_vec, recom_map, Data->positions, nsampTOT, copy_prob_new, copy_prob_newSTART, pop_vec, (ndonorpops+recipient_cond_ind), region_size, r, EMruns, all_versus_all_ind, haploid_ind, unlinked_ind, m, print_file9_ind, fout, fout3, fout9);

	      N_e_new=N_e_new+back_prob[0][nhaps_startpop];
	      num_regions_tot=num_regions_tot+back_prob[1][nhaps_startpop];

                 /* GET NEW COPY-PROBS BASED ON PAC SAMPLES: */
	      count = 0;
	      for (i = 0; i < (ndonorpops+recipient_cond_ind); i++)
		{
		  for (j = 0; j < ndonorhaps_vec[i]; j++)
		    {
		      total_back_prob[i] = total_back_prob[i] + back_prob[0][count];
		      total_prob = total_prob + back_prob[0][count];

		      total_back_probSTART[i] = total_back_probSTART[i] + back_prob[1][count];
		      total_probSTART = total_probSTART + back_prob[1][count];

		      total_counts[i] = total_counts[i] + back_prob[2][count];
		      total_lengths[i] = total_lengths[i] + back_prob[3][count];
		      total_differences[i] = total_differences[i] + back_prob[4][count];
		      //if (i < 2) printf("%d %d %d %d %d %d %lf %lf %lf %lf\n",m,r,h,i,j,count,back_prob[2][count],back_prob[3][count],total_counts[i],total_lengths[i]);

		      count = count + 1;
		    }
		}
	      for (i = 0; i < (ndonorpops+recipient_cond_ind); i++)
		{
		  total_region_counts[i] = total_region_counts[i] + back_prob[5][i];
		  total_squared_region_counts[i] = total_squared_region_counts[i] + back_prob[6][i];
		  snp_info_measure_final[i] = snp_info_measure_final[i] + back_prob[7][i];
		}
	    }
	  fprintf(fout3," %.10lf %.10lf\n",N_e,MutProb_vec_new[0]);
	  //if (estimate_mutationALL_ind==0) fprintf(fout3,"\n");
	  //if (estimate_mutationALL_ind==1) fprintf(fout3," %lf\n",MutProb_vec_new[0]);
	  if (estimate_recom_ind==1) N_e=N_e_new/(2.0-haploid_ind);

	  for (i = 0; i < (ndonorpops+recipient_cond_ind); i++)
	    {
	      copy_prob_pop[0][i] = total_back_prob[i]/total_prob;
	      copy_prob_pop[1][i] = total_back_probSTART[i]/total_probSTART;
	    }
                     /* RESET COPY-PROBS and MUTATION-PROBS: */
	      // (first check for probabilities of 0:)
	  for (i=0; i < (ndonorpops+recipient_cond_ind); i++)
	    {
	      if (copy_prob_pop[0][i] <= 0)
		copy_prob_pop[0][i] = small_copy_val*ndonorhaps_vec[i];

	      if (copy_prob_pop[1][i] <= 0)
		copy_prob_pop[1][i] = small_copy_val*ndonorhaps_vec[i];
	    }
	  total_prob = 0.0;
	  total_probSTART = 0.0;
	  for (j=0; j < (ndonorpops+recipient_cond_ind); j++)
	    {
	      total_prob = total_prob + copy_prob_pop[0][j];
	      total_probSTART = total_probSTART + copy_prob_pop[1][j];
	    }
	  for (j=0; j < (ndonorpops+recipient_cond_ind); j++)
	    {
	      copy_prob_pop[0][j] = copy_prob_pop[0][j]/total_prob;
	      copy_prob_pop[1][j] = copy_prob_pop[1][j]/total_probSTART;
	    }

	  if (estimate_copyprob_ind==1)
	    {
	      count = 0;
	      for (i=0; i < (ndonorpops+recipient_cond_ind); i++)
		{
		  for (j=0; j < ndonorhaps_vec[i]; j++)
		    {
		      copy_prob_new[count] = copy_prob_pop[0][i]/ndonorhaps_vec[i];
		      copy_prob_newSTART[count] = copy_prob_pop[1][i]/ndonorhaps_vec[i];
		      count = count + 1;
		    }
		}
	    }

	  if (estimate_mutation_ind==1)
	    {
	      count = 0;
	      for (i=0; i < (ndonorpops+recipient_cond_ind); i++)
		{
		  for (j=0; j < ndonorhaps_vec[i]; j++)
		    {
		      MutProb_vec_new[count] = total_differences[i]/(*p_nloci*(2-haploid_ind));
		      count = count + 1;
		    }
		}
	    }

	  if (estimate_mutationALL_ind==1)
	    {
	      sum_total_diff=0.0;
	      for (i=0; i < (ndonorpops+recipient_cond_ind); i++)
		sum_total_diff=sum_total_diff+total_differences[i]/(*p_nloci*(2-haploid_ind));
	      count = 0;
	      for (i=0; i < (ndonorpops+recipient_cond_ind); i++)
		{
		  for (j=0; j < ndonorhaps_vec[i]; j++)
		    {
		      MutProb_vec_new[count] = sum_total_diff;
		      count = count + 1;
		    }
		}
	    }

 	  if (r == (EMruns-1))
	    {
                /* print props, lengths, counts, and differences: */
	      if ((all_versus_all_ind==0) || ((all_versus_all_ind==1) && (donorlist_ind==0)))
		{
		  fprintf(fout2,"IND%d",m+1);
		  fprintf(fout4,"IND%d",m+1);
		  fprintf(fout5,"IND%d",m+1);
		  fprintf(fout6,"IND%d",m+1);
		  fprintf(fout7,"IND%d",m+1);
		  fprintf(fout8,"IND%d",m+1);
		}
	      if ((all_versus_all_ind==1) && (donorlist_ind==1))
		{
		  if ((2-haploid_ind)*indcount > ndonorval)
		    {
		      fgets(line,2047,fd3);
		      step=line;
		      reading(&step,"%s",waste);
		      reading(&step,"%d",&ndonorval);
		      indcount=1;
		    }
		  if (indcount_suppress_ind==0)
		    {
		      fprintf(fout2,"%s%d",waste,indcount);
		      fprintf(fout4,"%s%d",waste,indcount);
		      fprintf(fout5,"%s%d",waste,indcount);
		      fprintf(fout6,"%s%d",waste,indcount);
		      fprintf(fout7,"%s%d",waste,indcount);
		      fprintf(fout8,"%s%d",waste,indcount);
		    }
		  if (indcount_suppress_ind==1)
		    {
		      fprintf(fout2,"%s",waste);
		      fprintf(fout4,"%s",waste);
		      fprintf(fout5,"%s",waste);
		      fprintf(fout6,"%s",waste);
		      fprintf(fout7,"%s",waste);
		      fprintf(fout8,"%s",waste);
		    }
		}
	      if (all_versus_all_ind==0)
		{
		  fprintf(fout7," %d",num_regions_tot);
		  fprintf(fout8," %d",num_regions_tot);
		  for (j=0; j < (ndonorpops+recipient_cond_ind); j++)
		    {
		      fprintf(fout2," %lf",copy_prob_pop[0][j]);
		      fprintf(fout4," %lf",total_counts[j]);
		      fprintf(fout5," %lf",total_lengths[j]);
		      fprintf(fout6," %lf",total_differences[j]);
		      fprintf(fout7," %lf",total_region_counts[j]);
		      fprintf(fout8," %lf",total_squared_region_counts[j]);
		    }
		}
	      if (all_versus_all_ind==1)
		{
		  fprintf(fout7," %d",num_regions_tot);
		  fprintf(fout8," %d",num_regions_tot);
		  for (j=0; j < (ndonorpops+1); j++)
		    {
		      if (j < m)
			{
			  fprintf(fout2," %lf",copy_prob_pop[0][j]);
			  fprintf(fout4," %lf",total_counts[j]);
			  fprintf(fout5," %lf",total_lengths[j]);
			  fprintf(fout6," %lf",total_differences[j]);
			  fprintf(fout7," %lf",total_region_counts[j]);
			  fprintf(fout8," %lf",total_squared_region_counts[j]);
			}
		      if (j == m)
			{
			  fprintf(fout2," 0.00");
			  fprintf(fout4," 0.00");
			  fprintf(fout5," 0.00");
			  fprintf(fout6," 0.00");
			  fprintf(fout7," 0.00");
			  fprintf(fout8," 0.00");
			}
		      if (j > m)
			{
			  fprintf(fout2," %lf",copy_prob_pop[0][(j-1)]);
			  fprintf(fout4," %lf",total_counts[(j-1)]);
			  fprintf(fout5," %lf",total_lengths[(j-1)]);
			  fprintf(fout6," %lf",total_differences[(j-1)]);
			  fprintf(fout7," %lf",total_region_counts[(j-1)]);
			  fprintf(fout8," %lf",total_squared_region_counts[(j-1)]);
			}
		    }
		}
	      fprintf(fout2,"\n");
	      fprintf(fout4,"\n");
	      fprintf(fout5,"\n");
	      fprintf(fout6,"\n");
	      fprintf(fout7,"\n");
	      fprintf(fout8,"\n");
	    }

	}
      if (haploid_ind==0)
	{
	  if ((recipient_cond_ind==0) && (all_versus_all_ind==0)) DestroyData(Data);
	  if ((recipient_cond_ind==1) || (all_versus_all_ind==1)) DestroyData2(Data);
	}
      if (haploid_ind==1)
	{
	  if ((recipient_cond_ind==0) && (all_versus_all_ind==0)) DestroyDataHap(Data);
	  if ((recipient_cond_ind==1) || (all_versus_all_ind==1)) DestroyDataHap2(Data);
	}
      free(allelic_type_count_vec);
      free(found_vec);
      fclose(fd);

      indcount=indcount+1;
    }

  free(newhap);
  free(copy_prob_new);
  free(copy_prob_newSTART);
  free(MutProb_vec_new);
  for (i = 0; i < 8; i++)
    free(back_prob[i]);
  free(back_prob);
  free(total_back_prob);
  free(total_back_probSTART);
  free(ndonorhaps_vec);
  free(total_counts);
  free(total_lengths);
  free(total_differences);
  free(total_region_counts);
  free(total_squared_region_counts);
  free(snp_info_measure_final);

  if ((all_versus_all_ind==1) && (donorlist_ind==1))
    fclose(fd3);

  return(1);
}

void usage()
{
  printf("%s\n",helpfilestring);
}

int runprogram(int argc, char *argv[])
{
  unsigned int c1=0;
  while (c1<argc){ printf("%s\n",argv[c1]); c1++;};
  errormode=0;

  struct data_t *Data;
  int i,j;
  double bpval;
  double totaldonorprobs, leftoverprob, mut_rate_self;
  int log_lik_check;
  int ndonors, count, ndonorpops, ndonorpopsTEMP, ndonorval;
  int nind, nhaps, nsites, nhaps_startpop, cond_nhaps, cond_nind;
  int geno_find, recom_find, donorlist_find, EMiter_find, numsamples_find, outfile_find, ne_find, mut_find, num_found, copy_prop_em_find, recom_em_find, mutation_em_find, mutationALL_em_find, condition_recipient_inds_find, all_versus_all_ind, haploid_ind, region_size_find, unlinked_ind, prior_donor_probs_ind, mutation_rate_ind, indcount_suppress_ind, print_file9_ind;
  char *step;
  char line[2047];
  char waste[2047];
  FILE *fd, *fd2, *fd3, *fout, *fout2, *fout3, *fout4, *fout5, *fout6, *fout7, *fout8, *fout9;
  char * filename = malloc(1000 * sizeof(char));
  char * filenameGEN = malloc(1000 * sizeof(char));
  char * filenameDONORLIST = malloc(1000 * sizeof(char));
  char * filenameOUT = malloc(1000 * sizeof(char));
  srand((unsigned)time(NULL));

  /***********************************************************/
  // DEFAULT VALUES:

  int EMruns = 0;        // number of EMruns
  int samplesTOT = 10;   // number of final hidden-state samples desired after E-M is finished
  double N_e = 400000;   // scaling constant for recombination rate
  double GlobalMutRate = -9.0;    // global mutation rate per site
  double region_size = 100;    // number of chunks per region -- used to look at variability in copied chunks across regions in order to estimate "c" in Dan Lawson's fineSTRUCTURE

      /* 'NUISSANCE' PARAMETER DETAILS: */
  //double theta = 0.0001;
  //double theta = -9.0;
  double small_recom_val=0.000000000000001;    // lower limit for small genmap rates

  int start_val=0;
  int end_val=0;

  /************************************************************/
  geno_find=0;
  recom_find=0;
  donorlist_find=0;
  outfile_find=0;
  EMiter_find=0;
  numsamples_find=0;
  ne_find=0;
  mut_find=0;
  region_size_find=0;
  copy_prop_em_find=0;
  recom_em_find=0;
  mutation_em_find=0;
  mutationALL_em_find=0;
  condition_recipient_inds_find=0;
  all_versus_all_ind=0;
  haploid_ind=0;
  unlinked_ind=0;
  prior_donor_probs_ind=0;
  mutation_rate_ind=0;
  print_file9_ind=0;
  indcount_suppress_ind=0;
  num_found=0;
  jitter_locations=0;
  for (i=1; i < argc; i++)
    {
    if (strcmp(argv[i],"--internalerrors")==0)
	{
	    errormode=1;
	}
    if ((strcmp(argv[i],"-help")==0) || (strcmp(argv[i],"--help")==0) || (strcmp(argv[i],"-h")==0))
	{
	  usage();
	  stop_on_error(1);
	}
      if (strcmp(argv[i],"-g")==0)
	{
	  geno_find=1;
	  num_found=num_found+1;
	}
      if (strcmp(argv[i],"-r")==0)
	{
	  recom_find=1;
	  num_found=num_found+1;
	}
      if (strcmp(argv[i],"-f")==0)
	{
	  donorlist_find=1;
	  num_found=num_found+1;
	}
      if (strcmp(argv[i],"-i")==0)
	{
	  EMiter_find=1;
	  num_found=num_found+1;
	}
      if (strcmp(argv[i],"-s")==0)
	{
	  numsamples_find=1;
	  num_found=num_found+1;
	}
       if (strcmp(argv[i],"-n")==0)
	{
	  ne_find=1;
	  num_found=num_found+1;
	}
       if (strcmp(argv[i],"-M")==0)
	{
	  mut_find=1;
	  num_found=num_found+1;
	}
       if (strcmp(argv[i],"-k")==0)
	{
	  region_size_find=1;
	  num_found=num_found+1;
	}
       if (strcmp(argv[i],"-o")==0)
	{
	  outfile_find=1;
	  num_found=num_found+1;
	}
       if (strcmp(argv[i],"-ip")==0)
	 copy_prop_em_find=1;
       if (strcmp(argv[i],"-in")==0)
	 recom_em_find=1;
       if (strcmp(argv[i],"-im")==0)
	 mutation_em_find=1;
       if (strcmp(argv[i],"-iM")==0)
	 mutationALL_em_find=1;
       if (strcmp(argv[i],"-c")==0)
	 condition_recipient_inds_find=1;
       if (strcmp(argv[i],"-a")==0)
	 all_versus_all_ind=1;
       if (strcmp(argv[i],"-j")==0)
	 haploid_ind=1;
       if (strcmp(argv[i],"-u")==0)
	 unlinked_ind=1;
       if (strcmp(argv[i],"-p")==0)
	 prior_donor_probs_ind=1;
       if (strcmp(argv[i],"-b")==0)
	 print_file9_ind=1;
       if (strcmp(argv[i],"-y")==0)
	 indcount_suppress_ind=1;
       if (strcmp(argv[i],"-m")==0)
	 {
	   mutation_rate_ind=1;
	   num_found=num_found+1;
	 }
       if(strcmp(argv[i],"-J")==0)
	 {
	   jitter_locations=1;
	 }
    }
/*  if (argc != (num_found*2+copy_prop_em_find+recom_em_find+mutation_em_find+mutationALL_em_find+condition_recipient_inds_find+3*all_versus_all_ind+haploid_ind+unlinked_ind+prior_donor_probs_ind+print_file9_ind+1))
    {
      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
      printf("Unexpected number of arguments; expected %i but received %i\n",
             (num_found*2+copy_prop_em_find+recom_em_find+mutation_em_find+mutationALL_em_find+condition_recipient_inds_find+3*all_versus_all_ind+haploid_ind+unlinked_ind+prior_donor_probs_ind+print_file9_ind+1),
             argc);
      usage();
      stop_on_error(1);
    }*/
  if (donorlist_find==0)
    strcpy(filenameDONORLIST,"NULL");
  if (recom_find==0)
    strcpy(filename,"NULL");
  if (((geno_find==0) || (recom_find==0)) && (unlinked_ind==0)) { printf("Error with command line (Each of -g and -r MUST be specified if data are linked). Exiting...\n"); stop_on_error(1);}
  if ((geno_find==0) && (unlinked_ind==1)) { printf("Error with command line (-g MUST be specified). Exiting...\n"); stop_on_error(1);}
  if ((recom_find==1) && (unlinked_ind==1)) { printf("Data specified as containing unlinked sites (-u). Ignoring supplied recombination rate file....\n");}
  if ((all_versus_all_ind==1) && (condition_recipient_inds_find==1))
    {
      printf("You have selected both the -c and -a switches; you can choose at most one! If you have 0 donor haplotypes, you probably want '-a'. Exiting...\n");
      stop_on_error(1);
    }
  if ((mutation_em_find==1) && (mutationALL_em_find==1))
    {
      printf("You have specified to estimate a global mutation (emission) rate and population-specific mutation (emission) rates. Please choose only one of the '-im' and '-iM' switches. Exiting...\n");
      stop_on_error(1);
    }
  if ((mutation_rate_ind==1) && (mut_find==1))
    {
      printf("You have provided values for both a global mutation (emission) rate ('-M') and population-specific mutation (emission) rates ('-m'). Please choose only one of the '-m' and '-M' switches. Exiting...\n");
      stop_on_error(1);
    }
  if ((mutation_rate_ind==1) && (mutationALL_em_find==1))
    {
      printf("You have specified to estimate a global mutation (emission) rate; will ignore population-specific mutation (emission) rates in %s. If you wish to use donor-specific mutation rates, use the '-im' switch.\n",filenameDONORLIST);
    }
 for (i=1; i < argc; i++)
    {
      if (strcmp(argv[i],"-g")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      printf("In argument %i (%s)\n",i,argv[(i+1)]);
	      usage();
	      stop_on_error(1);
	    }
 	  strcpy(filenameGEN,argv[(i+1)]);
	  if (outfile_find==0)
	    {
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout = fopen(strcat(filenameOUT,".samples.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout2 = fopen(strcat(filenameOUT,".prop.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout3 = fopen(strcat(filenameOUT,".EMprobs.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout4 = fopen(strcat(filenameOUT,".chunkcounts.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout5 = fopen(strcat(filenameOUT,".chunklengths.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout6 = fopen(strcat(filenameOUT,".mutationprobs.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout7 = fopen(strcat(filenameOUT,".regionchunkcounts.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout8 = fopen(strcat(filenameOUT,".regionsquaredchunkcounts.out"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout9 = gzopen(strcat(filenameOUT,".copyprobsperlocus.out.gz"), "w");
	    }
	}
      if (strcmp(argv[i],"-r")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      stop_on_error(1);
	    }
	  strcpy(filename,argv[(i+1)]);
	}
      if (strcmp(argv[i],"-f")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      stop_on_error(1);
	    }
	  strcpy(filenameDONORLIST,argv[(i+1)]);
	}
     if (strcmp(argv[i],"-i")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      stop_on_error(1);
	    }
	  EMruns = atoi(argv[(i+1)]);
	  if (EMruns < 0)
	    {
	      printf("Number of EM runs must be at least 0. Exiting...\n");
	      stop_on_error(1);
	    }
	  if ((EMruns>0) && (copy_prop_em_find==0) && (recom_em_find==0) && (mutation_em_find==0) && (mutationALL_em_find==0))
	    {
	      printf("You have specified to perform E-M iterations, but have not specified which parameter(s) to maximize. If using '-i' switch, please specify at least one of '-in', '-ip', '-im', and/or '-iM'. Exiting...\n");
	      stop_on_error(1);
	    }
	}
      if (strcmp(argv[i],"-s")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      stop_on_error(1);
	    }
	  samplesTOT = atoi(argv[(i+1)]);
	  if (samplesTOT < 0)
	    {
	      printf("Number of samples must be >= 0. Exiting...\n");
	      stop_on_error(1);
	    }
	}
       if (strcmp(argv[i],"-n")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      stop_on_error(1);
	    }
	  N_e = atof(argv[(i+1)]);
	  if (N_e <= 0)
	    {
	      printf("Recombination scaling parameter N_e must be > 0. Exiting...\n");
	      stop_on_error(1);
	    }
	}
       if (strcmp(argv[i],"-M")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      stop_on_error(1);
	    }
	  GlobalMutRate = atof(argv[(i+1)]);
	  if (GlobalMutRate==0) GlobalMutRate=-9;
	  if (GlobalMutRate < 0)
	    printf("Mutation (emission) parameter must be > 0. Using Li & Stephens (2003) version of Watterson's estimate instead of user-supplied value...\n");
	}
       if (strcmp(argv[i],"-k")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      stop_on_error(1);
	    }
	  region_size = atof(argv[(i+1)]);
	  if (region_size < 1)
	    {
	      printf("Region_size must be >= 1. Exiting...\n");
	      stop_on_error(1);
	    }
	}
       if (strcmp(argv[i],"-m")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      stop_on_error(1);
	    }
	  mut_rate_self = atof(argv[(i+1)]);
	  if ((mut_rate_self > 1) && (condition_recipient_inds_find==1))
	    {
	      printf("Self mutation (emission) probability must be <= 1 (use a negative number to specify default). Exiting...\n");
	      stop_on_error(1);
	    }
	}
        if (strcmp(argv[i],"-a")==0)
	{
	  if ((argv[(i+1)][0] == '-') || (argv[(i+2)][0] == '-'))
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      stop_on_error(1);
	    }
	  start_val = atoi(argv[(i+1)]);
	  end_val = atoi(argv[(i+2)]);
	  if ((end_val < start_val) || (start_val < 0) || (end_val < 0))
	    {
	      printf("Invalid start_ind/stop_ind vals ('-a' switch). If you want to condition each individual on every other individual, use '-a 0 0'. Exiting...\n");
	      stop_on_error(1);
	    }
	  if (start_val > 0) start_val=start_val-1;
	}
       if (strcmp(argv[i],"-o")==0)
	 {
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      stop_on_error(1);
	    }
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout = fopen(strcat(filenameOUT,".samples.out"), "w");
  	   strcpy(filenameOUT,argv[(i+1)]);
	   fout2 = fopen(strcat(filenameOUT,".prop.out"), "w");
  	   strcpy(filenameOUT,argv[(i+1)]);
	   fout3 = fopen(strcat(filenameOUT,".EMprobs.out"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout4 = fopen(strcat(filenameOUT,".chunkcounts.out"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout5 = fopen(strcat(filenameOUT,".chunklengths.out"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout6 = fopen(strcat(filenameOUT,".mutationprobs.out"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout7 = fopen(strcat(filenameOUT,".regionchunkcounts.out"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout8 = fopen(strcat(filenameOUT,".regionsquaredchunkcounts.out"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout9 = gzopen(strcat(filenameOUT,".copyprobsperlocus.out.gz"), "w");
	 }
    }
  if (fout == NULL) {printf("error opening closing file1\n"); stop_on_error(1);}
  if (fout2 == NULL) {printf("error opening closing file2\n"); stop_on_error(1);}
  if (fout3 == NULL) {printf("error opening closing file3\n"); stop_on_error(1);}
  if (fout4 == NULL) {printf("error opening closing file4\n"); stop_on_error(1);}
  if (fout5 == NULL) {printf("error opening closing file5\n"); stop_on_error(1);}
  if (fout6 == NULL) {printf("error opening closing file6\n"); stop_on_error(1);}
  if (fout7 == NULL) {printf("error opening closing file7\n"); stop_on_error(1);}
  if (fout8 == NULL) {printf("error opening closing file8\n"); stop_on_error(1);}
  if (fout9 == NULL) {printf("error opening closing file9\n"); stop_on_error(1);}

      // open first file (to get information on haplotype numbers)
  fd = fopen(filenameGEN,"r");
  if (fd == NULL) { printf("error opening %s\n",filenameGEN); stop_on_error(1);}
  if (haploid_ind==0)
    {
      if ((condition_recipient_inds_find==0) && (all_versus_all_ind==0)) Data = ReadData(fd,0);
      if ((condition_recipient_inds_find==1)  || (all_versus_all_ind==1)) Data = ReadData2(fd,0);
    }
  if (haploid_ind==1)
    {
      if ((condition_recipient_inds_find==0) && (all_versus_all_ind==0)) Data = ReadDataHap(fd,0);
      if ((condition_recipient_inds_find==1)  || (all_versus_all_ind==1)) Data = ReadDataHap2(fd,0);
    }
  nind = (int) Data->nind;
  if (haploid_ind==0) nhaps = 2*nind;
  if (haploid_ind==1) nhaps = nind;
  nsites = (int) Data->nsnps;
  nhaps_startpop = (int) Data->nhaps_startpop;
  cond_nhaps = nhaps-nhaps_startpop;
  if (haploid_ind==0) cond_nind=cond_nhaps/2;
  if (haploid_ind!=0) cond_nind=cond_nhaps;
  if (ne_find==0) N_e=N_e/nhaps;
  if (start_val >= cond_nind)
    {
      printf("Your '-a' switch specifies to start with individual %d but there are only %d individuals in %s. Exiting....\n",start_val+1,cond_nind,filenameGEN);
      stop_on_error(1);
    }
  if ((nhaps_startpop != 0) && (donorlist_find==0))
    {
      printf("You have specified >0 donor haplotypes but have no donor-list (-f) file. Exiting....\n");
      stop_on_error(1);
    }
  if ((nhaps_startpop != 0) && (all_versus_all_ind==1))
    {
      printf("You have specified >0 donor haplotypes but have also specified you want everyone to copy everyone ('-a')? If you want to use the '-a' switch, please make the first row of %s a 0. Exiting....\n",filenameGEN);
      stop_on_error(1);
    }
  if ((nhaps_startpop == 0) && (all_versus_all_ind==0))
    {
      printf("You have specified 0 donor haplotypes, which is only allowed when specifying you want everyone to copy from everyone ('-a' switch). Exiting....\n");
      stop_on_error(1);
    }
  if ((nhaps_startpop == 0) && (donorlist_find==1))
    {
      //printf("You have specified 0 donor haplotypes but have also specified '-f' file. Are you sure you want 0 donor haplotypes?? Exiting....\n");
      //stop_on_error(1);
    }
  if ((nhaps_startpop == 0) && (condition_recipient_inds_find==1))
    {
      printf("You have specified 0 donor haplotypes but have also specified that you want to condition on self ('-c'). Maybe an odd thing to do (since you're only allowed to copy from one donor population) -- are you sure you did not want the '-a' switch, to tabulate copying of each individual from every other individual?\n");
      //stop_on_error(1);
    }
  if (nhaps_startpop >= nhaps)
    {
      printf("You've specified wanting %d donor haps, but you have only %d total haps\n",nhaps_startpop,nhaps);
      stop_on_error(1);
    }
  if (((int) Data->nind) != Data->nind)
    {
      printf("You cannot have fractions of individuals, but you have specified %lf total individuals.\n",Data->nind);
      stop_on_error(1);
    }
  if ((floor(cond_nhaps/2.0) != cond_nhaps/2.0) && (haploid_ind==0))
    {
      printf("You have specified %d recipient haplotypes, but this number must be even for diploid individuals. (Maybe you want to use '-j' to specify haploid individuals?)\n",cond_nhaps);
      stop_on_error(1);
    }
  if (all_versus_all_ind==1) printf("Will condition each individual on every other individual...\n");
  if (nhaps_startpop==0 && donorlist_find==1) printf("Will use %s for population labels in output (even though there are no donor haplotypes)\n",filenameDONORLIST);
  if (indcount_suppress_ind==1 && all_versus_all_ind==1) printf("Excluding counts from individual labels in output files....\n");

  if (haploid_ind==1) printf("Assuming all inds are haploid....\n");
  if (unlinked_ind==1)
    {
      printf("Assuming sites are unlinked....\n");
      //EMruns=0;
    }
  if ((prior_donor_probs_ind==1) && (condition_recipient_inds_find==0))
    printf("Using specified prior donor probs from input file....\n");
  if ((prior_donor_probs_ind==1) && (condition_recipient_inds_find==1))
    printf("Using specified prior donor probs from input file....(leftover probs will be assigned to own pop)\n");
  if ((mutation_rate_ind==1) && (mutationALL_em_find==0))
    printf("Using specified mutation rates from input file....\n");
  if (copy_prop_em_find==1)
    printf("Running E-M to estimate copying proportions....\n");
  if (recom_em_find==1)
    printf("Running E-M to estimate N_e....\n");
  if (mutation_em_find==1)
    printf("Running E-M to estimate mutation (emission) probabilities....\n");
  if (mutationALL_em_find==1)
    printf("Running E-M to estimate global mutation (emission) probability....\n");
  if (condition_recipient_inds_find==1)
    printf("Conditioning on own population's individuals (except self) in copying model....\n");
  if (jitter_locations==1)
    printf("Moving SNP locations where collisions occur....\n");
  printf(" Number of EM-runs = %d\n Number of samples = %d\n N_e value = %lf\n Region size = %lf\n",EMruns,samplesTOT,N_e,region_size);
  printf(" Global mutation value = %lf\n",GlobalMutRate);
  printf(" Number of donor haplotypes = %d\n Number of recipient haplotypes = %d\n",nhaps_startpop,nhaps-nhaps_startpop);

  fprintf(fout, "EM_iter = %d (N_e = %d / copy_prop = %d / mutation = %d / mutationGLOBAL = %d), nsamples = %d, N_e_start = %lf, region_size = %lf, genotype dataset = %s, genmap dataset = %s, donor-list dataset = %s\n", EMruns, recom_em_find, copy_prop_em_find, mutation_em_find, mutationALL_em_find, samplesTOT, N_e, region_size, filenameGEN,filename,filenameDONORLIST);

      // open third file (to get information on donor population hap numbers)
  if (nhaps_startpop != 0)
    {
      fd3 = fopen(filenameDONORLIST,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filenameDONORLIST); stop_on_error(1);}
      ndonorpops=-1;
      while(!feof(fd3))
	{
	  fgets(line,2047,fd3);
	  ndonorpops=ndonorpops+1;
	}
      fclose(fd3);
    }
  if (print_file9_ind==1) gzprintf(fout9,"pos");
  if (nhaps_startpop==0)
    {
      if (donorlist_find==1)
	{
	  fd3 = fopen(filenameDONORLIST,"r");
	  if (fd3 == NULL) { printf("error opening %s\n",filenameDONORLIST); stop_on_error(1);}
	  ndonorpopsTEMP=-1;
	  while(!feof(fd3))
	    {
	      fgets(line,2047,fd3);
	      ndonorpopsTEMP=ndonorpopsTEMP+1;
	    }
	  fclose(fd3);
	  fd3 = fopen(filenameDONORLIST,"r");
	  if (fd3 == NULL) { printf("error opening %s\n",filenameDONORLIST); stop_on_error(1);}
	  ndonors = 0;
	  for (i=0; i < ndonorpopsTEMP; i++)
	    {
	      fgets(line,2047,fd3);
	      step=line;
	      reading(&step,"%s",waste);

	      reading(&step,"%d",&ndonorval);
	      if (print_file9_ind==1)
		{
		  for (j=0; j < (ndonorval/(2-haploid_ind)); j++)
		    gzprintf(fout9," %s%d",waste,j+1);
		}
	      ndonors = ndonors + ndonorval;
	    }
	  fclose(fd3);
	  if (ndonors != nhaps)
	    {
	      printf("Number of donor haplotypes listed in %s must match total number of haplotypes in second line of %s if donor haplotypes = 0 and '-f' switch is specified! Exiting....\n",filenameDONORLIST,filenameGEN);
	      stop_on_error(1);
	    }
	}
      if (haploid_ind==0) ndonors=cond_nhaps-2;
      if (haploid_ind==1) ndonors=cond_nhaps-1;
      if (all_versus_all_ind==1)
	{
	  if (haploid_ind==0) ndonorpops=ndonors/2;
	  if (haploid_ind==1) ndonorpops=ndonors;
	}
      if (all_versus_all_ind==0) ndonorpops=0;
      if (donorlist_find==0)
	{
	  if (print_file9_ind==1)
	    {
	      for (i=0; i < (cond_nhaps/(2-haploid_ind)); i++) gzprintf(fout9," IND%d",i+1);
	    }
	}
    }
  int * ndonorhaps = malloc(ndonorpops * sizeof(int));
  double * ndonorprobs = malloc(ndonorpops * sizeof(double));
  double * ndonormutrates = malloc(ndonorpops * sizeof(double));
  if (nhaps_startpop != 0)
    {
      fd3 = fopen(filenameDONORLIST,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filenameDONORLIST); stop_on_error(1);}
      ndonors = 0;
      totaldonorprobs=0.0;
      for (i=0; i < ndonorpops; i++)
	{
	  fgets(line,2047,fd3);
	  step=line;
	  reading(&step,"%s",waste);
	  if (print_file9_ind==1) gzprintf(fout9," %s",waste);

	  reading(&step,"%d",&ndonorhaps[i]);
	  if (prior_donor_probs_ind==1)
	    {
	      reading(&step,"%lf",&ndonorprobs[i]);
	      if (ndonorprobs[i]<=0.0000000000000001)
		{
		  printf("Donor copying probabilities must be > 0. Exiting...\n");
		  stop_on_error(1);
		}
	      totaldonorprobs=totaldonorprobs+ndonorprobs[i];
	    }
	  if (mutation_rate_ind==1)
	    {
	      if (prior_donor_probs_ind==0) reading(&step,"%s",waste);
	      reading(&step,"%lf",&ndonormutrates[i]);
	      if (ndonormutrates[i]>1)
		{
		  printf("Donor mutation (emission) probabilities must be <= 1 (use a negative number to specify default). Exiting...\n");
		  stop_on_error(1);
		}
	    }
	  ndonors = ndonors + ndonorhaps[i];
	}
      fclose(fd3);
      if (ndonors != nhaps_startpop)
	{
	  printf("Number of donor haplotypes listed in %s does not match first line in %s! Exiting....\n",filenameDONORLIST,filenameGEN);
	  stop_on_error(1);
	}
      if (prior_donor_probs_ind==1)
	{
	  if (((totaldonorprobs > 1.00001) || (totaldonorprobs < 0.99999)) && (condition_recipient_inds_find==0))
	    {
	      printf("Probabilities across all donors in %s does not sum to 1.0 (instead sums to %lf)! Exiting....\n",filenameDONORLIST,totaldonorprobs);
	      stop_on_error(1);
	    }
	  if (((totaldonorprobs >= 1) || (totaldonorprobs <=0)) && (condition_recipient_inds_find==1))
	    {
	      printf("You've specified that there should be some prob of copying from your own pop, but probabilities across all donors in %s do not sum to something >0 and <1.0 (%lf)! Exiting....\n",filenameDONORLIST,totaldonorprobs);
	      stop_on_error(1);
	    }
	  leftoverprob=1.0-totaldonorprobs;
	}
      if (condition_recipient_inds_find==1)
	{
	  if (print_file9_ind==1) gzprintf(fout9," Self");
	  if (haploid_ind==0) ndonors=ndonors+cond_nhaps-2;
	  if (haploid_ind==1) ndonors=ndonors+cond_nhaps-1;
	}
    }
  if (print_file9_ind==1) gzprintf(fout9,"\n");
  if (nhaps_startpop==0)
    {
      if (all_versus_all_ind==1)
	{
	  for (i=0; i < ndonorpops; i++)
	    {
	      if (haploid_ind==0) ndonorhaps[i]=2;
	      if (haploid_ind==1) ndonorhaps[i]=1;
	    }
	}
    }
  printf(" num donor pops = %d\n",ndonorpops);

  int * pop_vec = malloc(ndonors * sizeof(int));
  double * MutProb_vec = malloc(ndonors * sizeof(double));
  double * copy_prob = malloc(ndonors * sizeof(double));
  double * copy_probSTART = malloc(ndonors * sizeof(double));

              // (0) INITIALIZE copy_prob, MutProb_vec, and pop_vec:
  if (prior_donor_probs_ind==0)
  {
    for (i=0; i < ndonors; i++) copy_prob[i] = 1.0/ndonors;
  }
  if (prior_donor_probs_ind==1)
  {
    count=0;
    for (i=0; i < ndonorpops; i++)
      {
	for (j=0; j < ndonorhaps[i]; j++)
	  {
	    copy_prob[count] = ndonorprobs[i]/ndonorhaps[i];
	    count = count + 1;
	  }
      }
    if (condition_recipient_inds_find==1)
      {
	if (haploid_ind==0)
	  {
	    for (j=0; j < (cond_nhaps-2); j++)
	      {
		copy_prob[count] = leftoverprob/(cond_nhaps-2);
		count=count+1;
	      }
	  }
 	if (haploid_ind==1)
	  {
	    for (j=0; j < (cond_nhaps-1); j++)
	      {
		copy_prob[count] = leftoverprob/(cond_nhaps-1);
		count=count+1;
	      }
	  }
      }
  }
 if ((mutation_rate_ind==0) || (mutationALL_em_find==1))
   {
     for (i=0; i < ndonors; i++)
       MutProb_vec[i]=GlobalMutRate;
   }
 if ((mutation_rate_ind==1) && (mutationALL_em_find==0))
   {
    count=0;
    for (i=0; i < ndonorpops; i++)
      {
	for (j=0; j < ndonorhaps[i]; j++)
	  {
	    MutProb_vec[count] = ndonormutrates[i];
	    count = count + 1;
	  }
      }
    if (condition_recipient_inds_find==1)
      {
	if (haploid_ind==0)
	  {
	    for (j=0; j < (cond_nhaps-2); j++)
	      {
		MutProb_vec[count] = mut_rate_self;
		count=count+1;
	      }
	  }
	if (haploid_ind==1)
	  {
	    for (j=0; j < (cond_nhaps-1); j++)
	      {
		MutProb_vec[count] = mut_rate_self;
		count=count+1;
	      }
	  }
      }
   }

  /*
  count = 0;
  for (i=0; i < ndonorpops; i++)
    {
      for (j=0; j < ndonorhaps[i]; j++)
	{
	  copy_prob[count] = (1.0/ndonorpops)/ndonorhaps[i];
	  count = count + 1;
	}
    }
  */

  for (i=0; i < ndonors; i++)
    copy_probSTART[i] = copy_prob[i];


  count = 0;
  for (i=0; i < ndonorpops; i++)
    {
      for (j=0; j < ndonorhaps[i]; j++)
	{
	  pop_vec[count] = i;
	  count = count + 1;
	}
    }
  if (condition_recipient_inds_find==1)
    {
      if (haploid_ind==0)
	{
	  for (j=0; j < (cond_nhaps-2); j++)
	    {
	      pop_vec[count] = ndonorpops;
	      count = count + 1;
	    }
	}
      if (haploid_ind==1)
	{
	  for (j=0; j < (cond_nhaps-1); j++)
	    {
	      pop_vec[count] = ndonorpops;
	      count = count + 1;
	    }
	}
    }


             // (i) GET SAMPLES, RUN E-M:
                  // open recomb map file:
  double * recom_map = malloc((Data->nsnps - 1) * sizeof(double));
  if ((unlinked_ind==1) && (recom_find==0))
    {
      for (j=0; j < (Data->nsnps-1); j++) recom_map[j]=-9.0;
    }
  if (recom_find==1)
    {
      fd2 = fopen(filename,"r");
      if (fd2 == NULL) { printf("error opening recom map input file: %s\n",filename); stop_on_error(1);}
      fgets(line,2047,fd2);   // header
      for (j=0; j < (Data->nsnps-1); j++)
	{
	  fgets(line,2047,fd2);
	  step=line;
	  reading(&step,"%lf",&bpval);    // basepair position
	  if (bpval != Data->positions[j])
	    {
	      if(jitter_locations)
		printf("Warning: genetic map position difference at basepair %d (%lf vs %lf). This will occur if jittering occurred or if using an invalid map.\n",j+1,bpval,Data->positions[j]);
	      else {
		printf("basepair positions do not match between %s and %s at basepair %d (%lf vs %lf). Exiting....\n",filename,filenameGEN,j+1,bpval,Data->positions[j]);
		stop_on_error(1);
	      }
	    }
	  reading(&step,"%lf",&recom_map[j]);
	  if (recom_map[j] >= 0 && recom_map[j] <= small_recom_val)
	    {
	      printf("Warning: recom rate very low at basepair %lf (%lf). Assuming recomb rate between this snp and next one is %lf....\n",Data->positions[j],recom_map[j],small_recom_val);
	      recom_map[j]=small_recom_val;
	    }
	  if (recom_map[j]<0)
	    {
	      printf("recom rate < 0 at basepair %lf. Assuming recomb rate of infinity between this snp and next one....\n",Data->positions[j]);
	      //printf("recom rate must be > 0 (basepair %lf)!! Exiting....\n",Data->positions[j]);
	      //stop_on_error(1);
	    }
	}
      fgets(line,2047,fd2);
      step=line;
      reading(&step,"%lf",&bpval);    // basepair position
      if (bpval != Data->positions[(Data->nsnps-1)])
	{
	  printf("basepair positions do not match between %s and %s at basepair %d. Exiting....\n",filename,filenameGEN,Data->nsnps);
	  stop_on_error(1);
	}
      fclose(fd2);
    }
               // check ordering of snps (only allowed to be less than previous position if recom_map<0 at position -- i.e. suggesting new chromosome):
  for (i=0; i < Data->nsnps; i++)
    {
      if (i > 0)
	{
	  if (Data->positions[i]<Data->positions[(i-1)] && (recom_map[(i-1)]>=0))
	    {
		  if(unlinked_ind==1){
			  printf("WARNING: positions in %s are not listed in increasing order between SNPs %i-%i (basepairs %lf and %lf). Continuing as loci are unlinked....\n",filenameGEN,i-1,i,Data->positions[(i-1)],Data->positions[i]);
		  }else{
			  printf("positions in %s are not listed in increasing order at basepairs %lf and %lf. Exiting....\n",filenameGEN,Data->positions[(i-1)],Data->positions[i]);
			  stop_on_error(1);
		  }
	    }
	}
    }

  if (haploid_ind==0)
    {
      if ((condition_recipient_inds_find==0)  && (all_versus_all_ind==0)) DestroyData(Data);
      if ((condition_recipient_inds_find==1)  || (all_versus_all_ind==1)) DestroyData2(Data);
    }
  if (haploid_ind==1)
    {
      if ((condition_recipient_inds_find==0)  && (all_versus_all_ind==0)) DestroyDataHap(Data);
      if ((condition_recipient_inds_find==1)  || (all_versus_all_ind==1)) DestroyDataHap2(Data);
    }
  fclose(fd);

  EMruns = EMruns + 1;

                 /* print-out headers for copy-props, chunk counts, lengths, and differences: */
  fprintf(fout2,"Recipient");
  fprintf(fout4,"Recipient");
  fprintf(fout5,"Recipient");
  fprintf(fout6,"Recipient");
  fprintf(fout7,"Recipient num.regions");
  fprintf(fout8,"Recipient num.regions");
  if (nhaps_startpop > 0)
    {
      fd3 = fopen(filenameDONORLIST,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filenameDONORLIST); stop_on_error(1);}
      for (j=0; j < ndonorpops; j++)
	{
	  fgets(line,2047,fd3);
	  step=line;
	  reading(&step,"%s",waste);
	  fprintf(fout2," %s",waste);
	  fprintf(fout4," %s",waste);
	  fprintf(fout5," %s",waste);
	  fprintf(fout6," %s",waste);
	  fprintf(fout7," %s",waste);
	  fprintf(fout8," %s",waste);
	}
      fclose(fd3);
    }
  if (condition_recipient_inds_find==1)
    {
      fprintf(fout2," Self");
      fprintf(fout4," Self");
      fprintf(fout5," Self");
      fprintf(fout6," Self");
      fprintf(fout7," Self");
      fprintf(fout8," Self");
    }
  if (all_versus_all_ind==1)
    {
      if (donorlist_find==0)
	{
	  for (j=0; j < (ndonorpops+1); j++)
	    {
	      fprintf(fout2," IND%d",j+1);
	      fprintf(fout4," IND%d",j+1);
	      fprintf(fout5," IND%d",j+1);
	      fprintf(fout6," IND%d",j+1);
	      fprintf(fout7," IND%d",j+1);
	      fprintf(fout8," IND%d",j+1);
	    }
	}
      if (donorlist_find==1)
	{
	  fd3 = fopen(filenameDONORLIST,"r");
	  if (fd3 == NULL) { printf("error opening %s\n",filenameDONORLIST); stop_on_error(1);}
	  ndonorpopsTEMP=-1;
	  while(!feof(fd3))
	    {
	      fgets(line,2047,fd3);
	      ndonorpopsTEMP=ndonorpopsTEMP+1;
	    }
	  fclose(fd3);
	  fd3 = fopen(filenameDONORLIST,"r");
	  if (fd3 == NULL) { printf("error opening %s\n",filenameDONORLIST); stop_on_error(1);}
	  for (j=0; j < ndonorpopsTEMP; j++)
	    {
	      fgets(line,2047,fd3);
	      step=line;
	      reading(&step,"%s",waste);
	      reading(&step,"%d",&ndonorval);
	      for (i=0; i < (ndonorval/(2-haploid_ind)); i++)
		{
		  if (indcount_suppress_ind==0)
		    {
		      fprintf(fout2," %s%d",waste,i+1);
		      fprintf(fout4," %s%d",waste,i+1);
		      fprintf(fout5," %s%d",waste,i+1);
		      fprintf(fout6," %s%d",waste,i+1);
		      fprintf(fout7," %s%d",waste,i+1);
		      fprintf(fout8," %s%d",waste,i+1);
		    }
		  if (indcount_suppress_ind==1)
		    {
		      fprintf(fout2," %s",waste);
		      fprintf(fout4," %s",waste);
		      fprintf(fout5," %s",waste);
		      fprintf(fout6," %s",waste);
		      fprintf(fout7," %s",waste);
		      fprintf(fout8," %s",waste);
		    }
		}
	    }
	  fclose(fd3);
	}
    }
  fprintf(fout2,"\n");
  fprintf(fout4,"\n");
  fprintf(fout5,"\n");
  fprintf(fout6,"\n");
  fprintf(fout7,"\n");
  fprintf(fout8,"\n");
  if (end_val==0 || end_val > cond_nind) end_val=cond_nind;
  log_lik_check = loglik(nhaps_startpop, &nsites, nhaps, N_e, recom_map, MutProb_vec, samplesTOT, ndonorpops, ndonorhaps, copy_prob, copy_probSTART, pop_vec, region_size, EMruns, copy_prop_em_find, recom_em_find, mutation_em_find, mutationALL_em_find, condition_recipient_inds_find, all_versus_all_ind, start_val, end_val, filenameGEN, filenameDONORLIST, donorlist_find, haploid_ind, unlinked_ind, print_file9_ind, indcount_suppress_ind, fout, fout2, fout3, fout4, fout5, fout6, fout7, fout8, fout9);
  if (log_lik_check != 1)
    {
      printf("Algorithm failed. Check input files and parameters. Exiting....\n");
      stop_on_error(1);
    }

  free(recom_map);
  free(filename);
  free(filenameGEN);
  free(filenameDONORLIST);
  free(copy_prob);
  free(copy_probSTART);
  free(pop_vec);
  free(ndonorhaps);
  free(ndonorprobs);
  free(ndonormutrates);
  free(MutProb_vec);

  fclose(fout);
  fclose(fout2);
  fclose(fout3);
  fclose(fout4);
  fclose(fout5);
  fclose(fout6);
  fclose(fout7);
  fclose(fout8);
  gzclose(fout9);

  return 0;
}
