/*! \file dataset.h
 *  \brief The dataset header file containing interfaces for the
 *         main dataset structure used in our toolbox
 */

/*! \brief Structure for a n dimenstional dataset composed of
 *         non missing floating point values.
 *  This structure is the main dataset structure used for all
 *  kinds of computations in this toolset. In general it is created
 *  initially by converting a fasta file to this internal data structure
 *  augmenting it with supplimental values.
 *  \var sequences The sequencies of a dataset orgeniced by
 *                 sequences[sequence_index .. from 0 to n_values-1]
 *                          [base position index]
 *  \var sequence_lengths The lengths of the sequences in this dataset
 *  \var binary_sequences A place holder for the use of binary sequences
 *  \var n_values the number of sequences stored in this dataset
 *  \var values additinal values for each sequence. Multiple values
 *              can be stored using the n_dimensions switch for each sequence.
 *              values are organized by 
 *              values[dimension][value ... from 0 to n_value-1],
 *              for computational efficiency. 
 *  \var n_dimensions the dimensionality of the data stored in values
 *  \var max_sequence_length the length of the longest sequence in
 *                           sequences
 */
typedef struct {
  char** sequences; 
  size_t* sequence_lengths;
  void** binary_sequences;
  int n_values; 
  float** values; 
  int n_dimensions; 
  size_t max_sequence_length;
} dataset;

typedef struct {
  size_t n_features;
  size_t n_samples;
} data_shape;

/*! \brief Structure to hold a consensus caluclation.
 *  This structure contains the statistics that lead to, as well as the
 *  conensus sequence calculated from a dataset
 *  \var seqeunece A string holding the consensus sequence.
 *  \var absolute_frequencies holds the absolute frequencies for each 
 *                            of the 4 bases at an indexed position. 
 *  \var relative_frequencies holds the relative frequencies for each
 *                            of the 4 bases at an indexed position.
 *  \var length The length of the sequence and hence of all arrays in this
 *              strucutre
 */
typedef struct {
  char* sequence;
  size_t* absolute_frequencies[4];
  double* relative_frequencies[4];
  size_t length;
} consens;

/*! \brief a structure for unique sequences
 */
typedef struct {
  int* u_seq;            /*!< array with indexes to
                          *    unique sequences in a dataset */
  int* multiplicities;   /*!< array with the number of duplicates */
  int n_seq;             /*!< number of unique sequences in dataset */
} unique_sequences;

/*! \brief a function that generates a dataset structure from a 
 *         fasta file.
 *  \param in an opened readable file pointer pointing to a fasta file.
 *  \return a dataset that holds the sequences in the fasta file.
 */
dataset dataset_from_fasta(FILE* in);

/*! \brief a function to free the memory used by a dataset
 *  This function frees an entire dataset.
 *  \param ds the dataset to be freed
 */
void free_dataset(dataset ds);

/*! \brief a function that frees supplimental values to a dataset 
 *         from the dataset.
 *  This function only frees the data held by the supplimental
 *  values to the dataset.
 *  \param ds the dataset that the supplimental values should be freed from.
 */
void free_values_from_dataset(dataset ds);

/*! \brief a function that frees sequences from a dataset.
 *  This function only frees the data held by the sequences from a dataset.
 *  \param ds the dataset that the supplimental values should be freed from.
 */
void free_sequences_from_dataset(dataset ds);

/*! \brief a function that loads projections from a PCA calculatxion 
 *         into a dataset as supplimental data. 
 *  \param projections an opened readable file pointer to 
 *                     a file containing the projections
 *                     i.e. obtained from the kmer2pca tool.
 *  \param dimensions the dimensions that are available in the file
 *                    to be loaded. Warning! Defining less or more dimensions
 *                    then those held in the file results in corrupt results.
 *  \param ds a pointer to the dataset that the projections shall be
 *            loaded to.
 */
void load_projections_from_file_into_dataset(FILE* projections,
					     size_t dimensions,
					     dataset* ds);
/*! \brief a function to load kmer vectors generated with fasta2kmer 
 *         as supplimental values into the dataset.
 *  This function loads a kmer represenation generated with fasta2kmer
 *  into the dataset. In order to do this one has to know the shape of the
 *  kmer file. The shape can be determined with the shape_from_kmer_file()
 *  function. One can load sequences from an initial fasta
 *  \param in_file a readable, opened file pointer to file containing the
 *                 kmer representation
 *  \param shape the datashape of the kmer file. 
 *  \return a dataset with the values of the kmer file.
 */
dataset load_kmer_from_file_into_dataset(FILE* in_file, data_shape shape);

/*! \brief A function to obtain the data shape of a kmer representation.
 *  \param infile a file descriptor of an opened readable file 
 *                containing the a kmer representation
 *  \return the shape of the kmer file
 */
data_shape shape_from_kmer_file(int infile);

/*! \brief A function to calculate a consensus sequence and consensus
 *         statistics from the sequences in a dataset.
 *  \param ds the dataset to calculate the consensus sequence 
 *  \return the consesus sequence and the statistics of frequencies of bases
 *          stored in a consens structure. 
 */ 
consens obtain_consens_from_dataset(dataset ds);

void free_consens(consens cs);

/*! \brief A function to print the consensus statistics in a readable
 *         manner to a file.
 *  \param f a pointer to an opened writeable file.
 *  \param cs the consensus structure whose contents shall be printed out.
 */
void print_consensus_statistics(FILE* f, consens cs);

/*! \brief A function to sort unique sequencies by frequency.
 *  \param us the unique sequences structure to be sorted by frequency.
 */
void sort_unique_sequences(unique_sequences us);

/*! \brief A function to obtain an index and multiplicities of unique
 *         sequences in a dataset.
 *  \param ds the dataset to find unique sequences in.
 */
unique_sequences get_sequence_multiplicities(dataset ds);

/*! \brief Write unique sequences in fasta format to a file.
 *  \param outfile An opened writeable file pointer the file to print the
 *                 unique sequences of a dataset to.
 *  \param ds The dataset that the unique sequences were calculated from
 *  \param us The index to the unique sequences in the dataset obtained with
 *            get_sequence_multiplicities.
 */
void write_unique_sequences(FILE* outfile, dataset ds, unique_sequences us);

/*! \brief Writes a whole dataset to a fasta file.
 *  \param outfile a opened writeable file pointer to write to.
 *  \param ds the dataset to be written to the file 
 */
void dataset_to_fasta(FILE* outfile, dataset ds);

/*! \brief Reverses sequences in a dataset using a binary mask
 *  \param ds a dataset containing sequences that will be manipulated
 *  \param binary_mask a binary array containing 1 where sequences
 *                     shall be reversed
 */
void reverse_sequences(dataset *ds, char* binary_mask);

/*! \brief Reverse complements sequences in a dataset using a binary mask
 *  \param ds a dataset containing sequences that will be manipulated
 *  \param binary_mask a binary array containing 1 where sequences
 *                     shall be translated into their reverse complement
 */
void reverse_complement_sequences(dataset* ds, char* binary_mask);
