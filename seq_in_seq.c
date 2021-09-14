#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<float.h>
#include<math.h>
#include<fftw3.h>
#include<string.h>
#include"dataset.h"

void file_error(char* path) {
  printf("failed to open file %s\n",path);
  _exit(1);
}


static inline char quick_to_base(size_t index) {
  switch(index){
  case 0:
    return 'A';
    break;
  case 1:
    return 'C';
    break;
  case 2:
    return 'G';
    break;
  case 3:
    return 'T';
    break;
  default:
    return 0;
  }
} 

static inline void sequence_to_fftw_complex(fftw_complex* c,char* sequence,
					    size_t length, size_t fill_length,
					    char base) {

  size_t i;

  for(i=0;i<length;i++) {
    if (sequence[i] == base) {
      c[i][0]=1.;
    } else {
      c[i][0]=0.;
    }
    c[i][1]=0.;
  }
  
  if(fill_length) {
    for(i=length;i<fill_length;i++) {
      c[i][0]=0;
      c[i][1]=0;
    }
  }
}


int main(int argc, char** argv) {

  size_t i,j;

  double a, b, c, d;
  
  FILE* fasta_f;
  FILE* target_f;

  dataset ds, t_ds;

  int sequence_index;

  size_t length;
  size_t flength;

  double d_length;
  double d_target_sequence_length;
  
  fftw_complex * ds_sequence_x, *ds_sequence_k;
  fftw_complex * ds_target_x, *ds_target_k;
  fftw_complex * convolved_x, *convolved_k;

  fftw_plan plan_forward_sequence;
  fftw_plan plan_forward_target;
  fftw_plan plan_reverse;

  double* module;
  
  char * reversed_target;

  double min,max, current_i_value;
  size_t i_plus_one;
  size_t min_index, max_index;

  int cutoff;

  int print_convolution;

  int gen_wisdom;
  int n_threads;

  

  fftw_init_threads();
  
  sscanf(argv[3],"%i",&sequence_index);
  sequence_index = sequence_index -1;

  sscanf(argv[4],"%i",&gen_wisdom);
  sscanf(argv[6],"%i",&cutoff);
  sscanf(argv[8],"%i",&n_threads);
  
  if(argc < 7) {
    printf("Arguments are: \n"
	   "  [fasta] sequence to search sequence in \n"
	   "  [fasta] sequence to be searched for \n"
	   "  [int] number of sequence to search in \n"
	   "        in the first fasta file. Beginning with 1\n"
	   "  [int] 1 created a wisdom file during this run " 
	   "(if you do not have one) \n"
	   "        0 use a wisdom file already available \n"
	   "  [wisdom] path of a wisdom file, always to be specified \n"
	   "  [int] cutoff - number of sequence changes to accept \n"
	   "  [string] chromosome specifier \n"
	   "  [int] number of threads fftw can use \n"
	   );
    return 1;
  }

  
  
  if ( NULL == (fasta_f = fopen(argv[1], "r"))) file_error(argv[1]);
  ds = dataset_from_fasta(fasta_f);
  fclose(fasta_f);

  if ( NULL == (target_f = fopen(argv[2], "r"))) file_error(argv[2]);
  t_ds = dataset_from_fasta(target_f);
  fclose(target_f);

  if(sequence_index > ds.n_values) {
    printf("Sequence Index is out of range\n");
    _exit(1);
  }
  
  reversed_target = (char*)malloc(sizeof(char)*t_ds.sequence_lengths[0]);
  
  j = 0;
  for(i=t_ds.sequence_lengths[0]-1;i>0;i--) {
    reversed_target[j] = t_ds.sequences[0][i];
    j++;
  }
  reversed_target[j] = t_ds.sequences[0][0];

  // reversed_target = t_ds.sequences[0];
  
  length = ds.sequence_lengths[sequence_index];

  d_length = (double)length;
  
  /*flength = 2;
  while(flength < length) {
    flength*=2;
    }*/
  flength=length;
  
  ds_sequence_x=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*flength);
  ds_sequence_k=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*flength);
  ds_target_x=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*flength);
  ds_target_k=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*flength);
  convolved_x=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*flength);
  convolved_k=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*flength);


  module = (double*)malloc(sizeof(double)*length);
  memset(module, 0, sizeof(double)*length);

  fftw_plan_with_nthreads(n_threads);

  if(gen_wisdom) {
  
    plan_forward_sequence = fftw_plan_dft_1d(flength,
					   ds_sequence_x, ds_sequence_k,
					   FFTW_FORWARD, FFTW_MEASURE);
  
    plan_forward_target = fftw_plan_dft_1d(flength,
					   ds_target_x, ds_target_k,
					   FFTW_FORWARD, FFTW_MEASURE);
  
    plan_reverse = fftw_plan_dft_1d(flength,convolved_k, convolved_x,
				    FFTW_BACKWARD, FFTW_MEASURE);
    
    
    if(!fftw_export_wisdom_to_filename(argv[5])) {
      printf("Wisdom file creation failed! Exiting! \n");
      return(1);
    }
  } else {
    if(!fftw_import_wisdom_from_filename(argv[5])) {
      printf("Wisdom file import failed! Exiting \n");
      return(1);
    }

    plan_forward_sequence = fftw_plan_dft_1d(flength,
					   ds_sequence_x, ds_sequence_k,
					     FFTW_FORWARD,
					     FFTW_WISDOM_ONLY | FFTW_MEASURE);
  
    plan_forward_target = fftw_plan_dft_1d(flength,
					   ds_target_x, ds_target_k,
					   FFTW_FORWARD,
					   FFTW_WISDOM_ONLY | FFTW_MEASURE);
    
    plan_reverse = fftw_plan_dft_1d(flength,convolved_k, convolved_x,
				    FFTW_BACKWARD,
				    FFTW_WISDOM_ONLY | FFTW_MEASURE);
  }
    
  for(j=0;j<4;j++) {
    
    sequence_to_fftw_complex(ds_sequence_x,ds.sequences[sequence_index],
			     length,flength,quick_to_base(j));
    
    sequence_to_fftw_complex(ds_target_x,reversed_target,
			     t_ds.sequence_lengths[0],
			     flength,quick_to_base(j));
    
    fftw_execute(plan_forward_sequence);
    fftw_execute(plan_forward_target);
    
    for(i=0;i<length;i++) {
      
      /* complex multiplication */
      
	a = ds_sequence_k[i][0];
	b = ds_sequence_k[i][1];
	c = ds_target_k[i][0];
	d = ds_target_k[i][1];
      
	convolved_k[i][0] = a*c-b*d;
	convolved_k[i][1] = b*c+a*d;
	
    }
    
    fftw_execute(plan_reverse);
    
    for(i=0;i<length;i++) {
      module[i] += sqrt(convolved_x[i][0]*convolved_x[i][0]
			+convolved_x[i][1]*convolved_x[i][1]);
    }
  }

  fftw_destroy_plan(plan_forward_sequence);
  fftw_destroy_plan(plan_forward_target);
  fftw_destroy_plan(plan_reverse);  

  
  fftw_free(ds_target_x);  
  fftw_free(ds_sequence_x);
  fftw_free(ds_target_k);    
  fftw_free(ds_sequence_k);  
  
  printf("track name=seq_pos description="
	 "\"Sequence position with %i mismatches\"\n",cutoff);

  d_target_sequence_length = (double)t_ds.sequence_lengths[0];
  min = (double)cutoff-0.1;
  max = (double)cutoff+0.1;
  for(i=0;i<length;i++) {
    current_i_value = module[i]/d_length;
    i_plus_one = i + 1;
    if(min < (d_target_sequence_length - current_i_value) &&
       max > (d_target_sequence_length - current_i_value)) {

      printf("%s %li %li\n",argv[7],i_plus_one-t_ds.sequence_lengths[0],
	     i_plus_one);
      
    }
    fflush(stdout);
  }
  
  
}	
	
      
      
      
    
  
  
  
  
  
  

 
  
