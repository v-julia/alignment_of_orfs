# alignment_of_orfs

This folder contains scripts for semi-authomatic generating of multiple sequence alignments of concatenated ORFs.


### *parser_gb.py* 

Converts file with nucleotide sequences in genbank format to fasta format with sequence names in the following format: GenbankAccessionNumber_annotation1_annotation2_..._annotationN, where annotation1, ..., annotationN are metadata from genbank qualifiers indicated in *-features* option. 


Sequences of length  beyond or above defined threshold will not be included in output file. 

Saves output file to the same directory as the input file.


#### Usage

```
parser_gb.py [-h] -input INPUT_FILE -min MIN_LENGTH -max MAX_LENGTH -f FEATURES
                    FEATURES


  -input <str>          Path to input file in GenBank format

                        
  -min   <int>          Minimal length of sequence. Sequences shorter than min
                        length will not be included in the final dataset
                        
  -max   <int>          Maximal length of sequence. Sequences longer than max
                        length will not be included in the final dataset
  -features   <str>     string with qualifiers to retrieve from GenBank annotation,\
                        e.g. 'country,host,collection_date'


```

E.g. Command
```
parser_gb.py -input seqs.gb -min 1000 -max 2000 -f 'host,country'

```
will produce *seqs.fasta* with sequences of length range 1000-2000 nt. The names of sequences in fasta file will be in the following format ">GenbankAccession_host name_country".

If *country_map.csv*  is in the local directory, sequence names will contain abbreviations rather than full countries names.

If *host_map.csv* is in the local directory, sequence names will contain short host names inducated in *host_map.csv*.

### genecds_field_distribution.py

Parses file with nucleotide sequences in genbank format, prints sorted list of values found in *product* and *gene* qualifiers in CDS field
#### Usage
```
genecds_field_distribution.py [-h] -input INPUT_FILE [-odir OUTPUT_DIR]
                                     -oname OUTPUT_NAME

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file in GenBank format
  -odir OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory to save the output file
  -oname OUTPUT_NAME, --output_name OUTPUT_NAME
                        Name of output file
```

### get_orfs_coord.py

Retrieves the coordinates of ORFs from file with nucleotide sequences in genbank format.

#### Usage
```
get_orfs_coord.py [-h] -input INPUT_FILE -orf_map ORF_MAP_FILE [-r]

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file
  -orf_map ORF_MAP_FILE, --orf_map_file ORF_MAP_FILE
                        Csv-file with short codes for ORFs
  -r, --remove_exceptions
                        Remove exceptions. 
```


### get_slice.py 

Creates the slice of multiple sequence alignment in fasta format using defined coordinates and saves it to fasta file
    
#### Usage
```
get_alignment_slice.py [-h] -input INPUT_FILE -s START -e END

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file
  -s START, --start START
                        Start of slice
  -e END, --end END     End of slice
```

### *remove_similar.py*

Calculates p-distances between nucleotide sequences pairs in loop. If p-distance < min_distance and p-distance > max_distance 
removes the sequence with higher serial number. 

Saves output file with reduced sequence set to the directory of input file.


#### Usage

```
remove_similar.py [-h] -input INPUT_FILE -min MIN_DISTANCE -max MAX_DISTANCE


  -input INPUT_FILE, --input_file INPUT_FILE
                        Path to input file in fasta format
  -min MIN_DISTANCE, --min_distance MIN_DISTANCE
                        Minimal pairwise distance between sequences. If
                        p-distance is lower than min_distance sequence with
                        higher serial number will be removed from the dataset
  -max MAX_DISTANCE, --max_distance MAX_DISTANCE
                        Maximal pairwise distance between sequences. If
                        p-distance is higher than max_distance sequence with
                        higher serial number will be removed from the dataset
```


### split_genome.py

Excises ORF sequences from sequences in *input_file* in genbank format according to coordinates in *coord_file*

#### Usage
```
split_genome.py [-h] -input INPUT_FILE -coord COORD_FILE

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file in fasta-format
  -coord COORD_FILE, --coord_file COORD_FILE
                        csv-file with coordinates of ORFs retrived from CDS field of genbank file
```

### trans_alignment.py

### gap_in_row.py

Removes sequences with many gaps in row from file in fasta format. Output file in fasta format is saved in the save folder as input file.

#### Usage
```
gap_in_row.py [-h] -input INPUT_FILE -gap_count GAP_COUNT

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file
  -gap_count GAP_COUNT, --gap_count GAP_COUNT
                        Amount of gaps in a row
```


### join_al.py

### create_meta.py




## Requirements

* python 3
* Biopython
* numpy
* matplotlib

