#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<unistd.h>
#include<math.h>

#include"dataset.h"
#include"binary_array.h"
#include"colors.h"

typedef struct {
  int buffer[2];
} sorthelper;

int comp(const void *a, const void *b) {

  sorthelper* sa = (sorthelper*)a;
  sorthelper* sb = (sorthelper*)b;

  return((sa->buffer[1]-sb->buffer[1]));
}

void sort_unique_sequences(unique_sequences us) {

  int i;
  sorthelper* helper = (sorthelper*)malloc(sizeof(sorthelper)*us.n_seq);

  for(i=0;i<us.n_seq;i++) {
    helper[i].buffer[0] = us.u_seq[i];
    helper[i].buffer[1] = us.multiplicities[i];
  }

  qsort(helper, us.n_seq, sizeof(sorthelper), comp);

  for(i=0;i<us.n_seq;i++) {
    us.u_seq[i] = helper[i].buffer[0];        
    us.multiplicities[i] = helper[i].buffer[1];
  }
}

unique_sequences get_sequence_multiplicities(dataset ds) {

  unique_sequences us;
  
  char* visited_sites = alloc_and_set_zero_binary_array(ds.n_values);

  int i,j, u_seq_counter;

  char* current_sequence;

  us.u_seq = (int*)malloc(sizeof(int)*ds.n_values);
  us.multiplicities = (int*)malloc(sizeof(int)*ds.n_values);

  for(i=0; i<ds.n_values;i++) {
    us.multiplicities[i] = 1;
  }
  
  u_seq_counter=0;
  for(i=0; i<ds.n_values;i++) {
    if (!get_value_in_binary_array_at_index(visited_sites,i)) {
      us.u_seq[u_seq_counter] = i;
      set_value_in_binary_array_at_index(visited_sites,i);
      for(j=0; j<ds.n_values;j++) {
	if(!get_value_in_binary_array_at_index(visited_sites,j)) {
	  if (!strcmp(ds.sequences[i],ds.sequences[j]) &&
	      (strlen(ds.sequences[i]) == strlen(ds.sequences[j]))) {
	    set_value_in_binary_array_at_index(visited_sites,j);
	    us.multiplicities[u_seq_counter]++;
	  }
	}
      }
      u_seq_counter++;
    }
  }

  us.u_seq = (int*)realloc(us.u_seq, sizeof(int)*u_seq_counter);
  us.multiplicities =
    (int*)realloc(us.multiplicities, sizeof(int)*u_seq_counter);
  us.n_seq = u_seq_counter;
  return(us);
}


void dataset_to_fasta(FILE* f, dataset ds) {

  int j, k;

  for(j=0;j < ds.n_values; j++) {
    fprintf(f,">sequence_%i\n",j);
    for ( k = 0;
	  k < (ds.sequence_lengths[j]-1);
	  k++) {
      if( k != 0 && k%50 == 0 ) {
	fprintf(f, "\n");
	fputc(ds.sequences[j][k],f);
      } else {
	fputc(ds.sequences[j][k],f);
      }
    }
    if ( (ds.sequence_lengths[j]-1) != 0 &&
	 (ds.sequence_lengths[j]-1) %50 == 0 ) {
      fprintf(f, "\n");
      fputc(ds.sequences[j][k],f);
    } else {
      fputc(ds.sequences[j][k],f);
    }
    fprintf(f, "\n");
  }
}

void reverse_complement_sequences(dataset* ds, char* binary_mask) {
  int i,j;
  char* buffer = (char*)malloc(sizeof(char)*ds->max_sequence_length);
  
  for (i = 0 ; i<ds->n_values; i++) {
    if (get_value_in_binary_array_at_index(binary_mask,i)) {
      memcpy(buffer, ds->sequences[i], ds->sequence_lengths[i]);
      for(j = 0; j<ds->sequence_lengths[i];j++) {
	switch(buffer[ds->sequence_lengths[i]-1-j]) {
	case 'A':
	  ds->sequences[i][j] = 'T';
	  break;
	case 'C':
	  ds->sequences[i][j] = 'G';
	  break;
	case 'G':
	  ds->sequences[i][j] = 'C';
	  break;
	case 'T':
	  ds->sequences[i][j] = 'A';
	  break;
	default:
	  ds->sequences[i][j] = buffer[ds->sequence_lengths[i]-1-j];
	  break;
	}
      }
    }
  }
  free(buffer);
}	

void reverse_sequences(dataset* ds, char* binary_mask) {
  int i,j;

  char* buffer = (char*)malloc(sizeof(char)*ds->max_sequence_length);
  
  for (i = 0 ; i<ds->n_values; i++) {
    if (get_value_in_binary_array_at_index(binary_mask,i)) {
      memcpy(buffer, ds->sequences[i], ds->sequence_lengths[i]);
      for (j = 0 ; j<ds->sequence_lengths[i]; j++) {
	ds->sequences[i][j] = buffer[ds->sequence_lengths[i]-1-j];
      }
    }
  }
  free(buffer);
}

void write_unique_sequences(FILE* outfile, dataset ds, unique_sequences us) {

  int i,k;
  for(i=0; i<us.n_seq; i++) {
    fprintf(outfile,">sequence_%i_mutiplicity_%i\n",
	    us.u_seq[i],
	    us.multiplicities[i]);
    for( k = 0;
	 k < (strlen(ds.sequences[us.u_seq[i]]))-1;
	 k++) {
      if( k != 0 && k%50 == 0 ) {
	fprintf(outfile, "\n");
	fputc(ds.sequences[us.u_seq[i]][k],outfile);
      } else {
	fputc(ds.sequences[us.u_seq[i]][k],outfile);
      }
    }
    if( (strlen(ds.sequences[us.u_seq[i]])-1) != 0 &&
	(strlen(ds.sequences[us.u_seq[i]])-1) %50 == 0 ) {
      fprintf(outfile, "\n");
      fputc(ds.sequences[us.u_seq[i]][k], outfile);
    } else {
      fputc(ds.sequences[us.u_seq[i]][k], outfile);
    }
    fprintf(outfile,"\n");
  }
}

void char_sequence_to_binary(char* c_seq, char* b_seq, size_t seq_len) {

  size_t i;
  int j;
  char bin;
  for(i = 0; i< seq_len; i++) {

    switch(c_seq[i]) {
    case 'A':
      bin = 0;
      break;
    case 'C':
      bin = 1;
      break;
    case 'G':
      bin = 2;
      break;
    case 'T':
      bin = 3;
      break;
    default:
      printf("unhandled character for binary conversion %c\n", c_seq[i]);
      _exit(1);
    }

    j=i%4;

    b_seq[i/4] |= (bin << 2*j);
  }
}

char get_char_from_binary_sequence_at_index(char* b_seq, size_t idx) {

  size_t b_idx = idx/4;
  size_t b_off = idx%4;

  char test_a = (1 << (b_off*2));
  char test_b = (1 << (b_off*2)+1);

  char res = 0;

  char ret_val;
  
  if(test_a == (b_seq[b_idx] & test_a)) res += 1;
  if(test_b == (b_seq[b_idx] & test_b)) res += 2;

  switch (res) {
    case 0:
      ret_val = 'A';
      break;
    case 1:
      ret_val = 'C';
      break;
    case 2:
      ret_val = 'G';
      break;
    case 3:
      ret_val = 'T';
      break;
    }
  return(ret_val);
}
  

void add_binary_sequences_to_dataset(dataset* ds) {

  size_t i;
  char** c_seq = ds->sequences;
  size_t* lengths = ds->sequence_lengths;
  char** b_seq = malloc(sizeof(ds->n_values)*sizeof(char*));
  
  for(i = 0 ; i< ds->n_values; i++ ) {
    b_seq[i] = (char*)malloc(sizeof(char)*(lengths[i]/4+1));
    char_sequence_to_binary(c_seq[i], b_seq[i], lengths[i]);
  }
}

dataset dataset_from_fasta(FILE* in) {

  char* line = NULL;
  size_t line_size = 0;

  int sequences = 0;

  dataset ds;

  char linebuffer[2000];
  size_t linebuffer_length;
  size_t sequence_length;
  int i;
  size_t length_before_addition_of_line;
  
  rewind(in);

  while ( -1 != getline(&line, &line_size, in) ) {
    if( line[0] == '>' ) sequences++;
  }

  rewind(in);

  ds.sequence_lengths = (size_t*)malloc(sizeof(size_t*)*sequences);

  ds.sequences = (char**)malloc(sizeof(char*)*sequences);
  ds.n_values = sequences;

  sequences = -1;

  ds.max_sequence_length = 0;
  
  while ( -1 != getline(&line, &line_size, in) ) {
    
    if ( line[0] == '>') {
      /* new sequence in file */
      sequences++;

      /* redress memory usage of previous sequence */
      if(sequences > 0) {
	ds.sequences[sequences-1] =
	  (char*)realloc(ds.sequences[sequences-1],
			 sizeof(char)*(ds.sequence_lengths[sequences-1]+1)); 
      }
			 						   
      ds.sequences[sequences] = (char*)malloc(sizeof(char)*1001);
      ds.sequences[sequences][0] = 0;
      ds.sequence_lengths[sequences] = 0;

    } else {
      /* continuation of current sequence */

      sscanf(line,"%s",linebuffer);

      linebuffer_length = strlen(linebuffer);
      length_before_addition_of_line = ds.sequence_lengths[sequences];
      ds.sequence_lengths[sequences] += linebuffer_length;
      if(ds.sequence_lengths[sequences] > 1000) {
	ds.sequences[sequences] =
	  (char*)realloc(ds.sequences[sequences],
			 sizeof(char)*(ds.sequence_lengths[sequences]+1));
      }
      for(i=0; i < linebuffer_length; i++) {
	linebuffer[i] = toupper(linebuffer[i]);
      }
      memcpy(ds.sequences[sequences]+length_before_addition_of_line,
	     linebuffer,
	     linebuffer_length);
      ds.sequences[sequences][length_before_addition_of_line
			      +linebuffer_length] = 0; 
    }
    
  }
  for(i=0;i<ds.n_values;i++) {
    if(ds.max_sequence_length < ds.sequence_lengths[i]) {
      ds.max_sequence_length = ds.sequence_lengths[i];
    }
  }
  free(line);
  return(ds);
}  
  
data_shape shape_from_kmer_file(int infile) {

  data_shape s;

  unsigned char current_character;

  off_t size;

  int i;

  s.n_features = 0;
  s.n_samples = 0;

  size = lseek(infile, 0, SEEK_END);
  lseek(infile, 0, SEEK_SET);

  char* f_buffer = (char*)malloc(sizeof(char)*size);

  if (size != read(infile, f_buffer, size)) {
    printf("Could not read file in order to read obtain data shape\n");
    _exit(1);
  }

  i = 0;
  while( f_buffer[i] != '\n') {
    if ( f_buffer[i] == '\t') s.n_features++;
    i++;
  }

  while( i < size ) {
    if ( f_buffer[i] == '\n' ) s.n_samples++;
    i++;
  }

  return(s);
}


dataset load_kmer_from_file_into_dataset(FILE* in_file, data_shape shape) {

  dataset ds;
  int i, j;

  char buffer[1024];

  ds.n_dimensions = shape.n_features;

  ds.n_values = shape.n_samples;

  ds.values = (float**)malloc(sizeof(float*)*ds.n_dimensions);

  for(i = 0 ; i < ds.n_dimensions; i++) {

#if defined(__AVX__) || defined(__SSE__)
    if(0 != posix_memalign((void**)&ds.values[i],
#if defined(__AVX__)
			   32,
#elif defined(__SSE__)
			   16,
#endif
			   sizeof(float)*ds.n_values)) {
      printf("Error allocating aligned memory for kmers \n");
      _exit(1);
    }
#else
    ds.values[i] =
      (float*)malloc(sizeof(float)*ds.n_values);
#endif
  }

  rewind(in_file);

  for(i=0;i<shape.n_samples;i++) {
    fscanf(in_file,"%s", buffer);
    for(j=0;j<shape.n_features;j++) {
      fscanf(in_file,"%f", ds.values[j]+i);
    }
  }
  return(ds);
}

void load_projections_from_file_into_dataset(FILE* projections,
					     size_t dimensions,
					     dataset* ds) {

  size_t i,j;

  size_t counter;

  ds->n_dimensions = dimensions;
  ds->values = (float**)malloc(sizeof(float*)*ds->n_dimensions);

  for(i = 0 ; i < ds->n_dimensions; i++) {

#if defined(__AVX__) || defined(__SSE__)
    if(0 != posix_memalign((void**)&ds->values[i],
#if defined(__AVX__)
			   32,
#elif defined(__SSE__)
			   16,
#endif
			   sizeof(float)*(ds->n_values+(8-(ds->n_values%8))))) {
      printf("Could not allocate memory to load projections \n");
      _exit(1);
    }
#else
    if (NULL == (ds->values[i] = (float*)malloc(sizeof(float)*ds->n_values))) {
      printf("Could not allocate memory to load projections \n");
      _exit(1);
    }    
#endif
  }

  for(i=0;i<ds->n_values;i++) {
    counter = 0;
    for(j=0;j<dimensions;j++) {
      counter += fscanf(projections, "%f", ds->values[j]+i);
    }
    if(counter != dimensions) {
      printf("Error reading projections file \n");
      _exit(1);
    }
  }
}

void free_consens(consens cs) {

  int i;
  for(i=0;i<4;i++) {
    free(cs.absolute_frequencies[i]);
    free(cs.relative_frequencies[i]);
  }
}

consens obtain_consens_from_dataset(dataset ds) {

  int i,j;
  consens cs;

  int bin,val;
  int maximum;

  cs.sequence = (char*)malloc(sizeof(char)*ds.max_sequence_length);
  
  for(i=0;i<4;i++) {
    cs.absolute_frequencies[i] =
      (size_t*)malloc(sizeof(size_t)*ds.max_sequence_length);
    cs.relative_frequencies[i] =
      (double*)malloc(sizeof(double)*ds.max_sequence_length);
    memset(cs.absolute_frequencies[i],
	   0, sizeof(size_t)*ds.max_sequence_length);
  }

  for(i=0;i<ds.n_values;i++) {
    for(j=0;j<ds.sequence_lengths[i];j++) {
      switch(ds.sequences[i][j]) {
      case 'A':
	 cs.absolute_frequencies[0][j]++;
	 break;
      case 'C':
	cs.absolute_frequencies[1][j]++;
	break;
      case 'G':
	cs.absolute_frequencies[2][j]++;
	break;
      case 'T':
	cs.absolute_frequencies[3][j]++;
	break;
      default:
	break;
      }
    }
  }

  for(i=0;i<4;i++) {
    for(j=0;j<ds.max_sequence_length;j++) {
      cs.relative_frequencies[i][j] =
	(double)cs.absolute_frequencies[i][j]/(double)ds.n_values;
    }
  }

  for(j=0;j<ds.max_sequence_length;j++) {
    maximum = 0;
    for(i=0;i<4;i++) {
      if ( cs.absolute_frequencies[i][j] > maximum ) {
	maximum = cs.absolute_frequencies[i][j];
	bin = i;
      }
    }

    switch (bin) {
    case 0:
      val = 'A';
      break;
    case 1:
      val = 'C';
      break;
    case 2:
      val = 'G';
      break;
    case 3:
      val = 'T';
      break;
    default:
      printf("unhandled integer for consensus convers %c\n", bin);
      _exit(1);
    }

    cs.sequence[j] = val;
  }
  cs.length = ds.max_sequence_length;
  return(cs);
}

void print_consensus_statistics(FILE* f, consens cs) {

  int i,j;

  double max, mean, sigma;

  fprintf(f,
	  "Base       A          C          G          T       A       C"
	  "       G       T\n");
  
  for(i=0;i<cs.length;i++) {

    max = 0;
    for(j=0;j<4;j++) {
      if (max < cs.relative_frequencies[j][i]) {
	max = cs.relative_frequencies[j][i];
      }
    }

    if (max < .85) {
      fprintf(f,ANSI_COLOR_RED);
    }
    
    fprintf(f,"%c %10lu %10lu %10lu %10lu %7.5lf %7.5lf %7.5lf %7.5lf\n",
	    cs.sequence[i],
	    cs.absolute_frequencies[0][i],
	    cs.absolute_frequencies[1][i],
	    cs.absolute_frequencies[2][i],
	    cs.absolute_frequencies[3][i],

	    cs.relative_frequencies[0][i],
	    cs.relative_frequencies[1][i],
	    cs.relative_frequencies[2][i],
	    cs.relative_frequencies[3][i]);

    if (max < .85) {
      fprintf(f,ANSI_COLOR_RESET);
    }
  }

  mean = 0;
  for(i=0;i<cs.length;i++) {
    max = 0;
    for(j=0;j<4;j++) {
      if (max < cs.relative_frequencies[j][i]) {
	max = cs.relative_frequencies[j][i];
      }
    }
    mean += max;
  }
  mean /= cs.length;

  sigma = 0;
  for(i=0;i<cs.length;i++) {
    max = 0;
    for(j=0;j<4;j++) {
      if (max < cs.relative_frequencies[j][i]) {
	max = cs.relative_frequencies[j][i];
      }
    }
    sigma += (max-mean)*(max-mean);
  }
  sigma /= cs.length;
  sigma = sqrt(sigma);

  printf("    MEAN MAX ACCURACY: %lf SIGMA: %lf\n", mean, sigma);
  
}

void free_values_from_dataset(dataset ds) {
  int i;
  for(i=0;i<ds.n_dimensions;i++) {
    free(ds.values[i]);
  }
  free(ds.values);
}

void free_sequences_from_dataset(dataset ds) {
  int i;
  for(i=0;i<ds.n_values;i++) {
    free(ds.sequences[i]);
  }
  free(ds.sequence_lengths);
  free(ds.sequences);
}

void free_dataset(dataset ds) {
  free_values_from_dataset(ds);
  free_sequences_from_dataset(ds);
}
