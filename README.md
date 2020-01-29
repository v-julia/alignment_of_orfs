# GenAlignment

This folder contains scripts for handling sets of biological sequences and generating alignments.

### *GenAlignment.py* 
Generates alignment of genomic region of interest in fasta-format with 
annotated sequence names (GenbankID_strain_isolate_country_collectionYear) from file with sequences in GenBank-format.

The input parameters should be defined in *config.txt* file located in the same directory as script:
* path to input file with fasta sequences in GenBank-format
* path to file with reference sequence
* positions of region of interest in reference sequence
* maximal and minimal length to retrieve from GenBank file
* method of removing some sequences from the dataset:
    * if p-distance between two sequences is beyond defined cut-off the one with higher serial order will be removed
    * set of sequences is divided into groups based on the first 5 characters in GenBank ID, then if the size of groups \
    exceeds m k% of randomly chosen sequences will be removed
* paths to mafft and blastn programs

File *country_map.csv* with countries' abbreviations located in the same directory as script can be used to replace countries names with abbreviations in sequences names.

#### Usage

```
python GenAlignment.py
```

*GenAlignment.py* uses several modules which can be used as independent scripts.

### *parser_gb.py* 
Converts file with sequences in GenBank-format to fasta-format with sequence names in the following format: GenbankAccessionNumber_strain_isolate_country_collectionYear. 

Sequences which lengths are beyond or above defined threshold wll not be included in ouput-file. 

Saves output-file in the directory of input-file.
 
If *country_map.csv* with countries' abbreviations is in the local directory, sequence names will contain abbreviations rather than full contries names.

#### Usage

```
parser_gb.py [-h] -input INPUT_FILE -min MIN_LENGTH -max MAX_LENGTH -f
                    FEATURES

  -input <str>          Path to input-file in GenBank-format
                        
  -min   <int>          Minimal length of sequence. Sequences shorter than min
                        length will not be included in the final dataset
                        
  -max   <int>          Maximal length of sequence. Sequences longer than max
                        length will not be included in the final dataset
  -features   <str>     Maximal length of sequence. Sequences longer than max
                        length will not be included in the final dataset


```

### *remove_similar.py*

Calculates p-distances between sequences pairs in loop, if p-distance < min_distance and p-distance > max_distance 
removes the sequence with higher serial number. 

Saves output-file with reduced sequence set in the directory of input-file.

#### Usage

```
remove_similar.py [-h] -input INPUT_FILE -min MIN_DISTANCE -max MAX_DISTANCE


  -input INPUT_FILE, --input_file INPUT_FILE
                        Path to input file in fasta-format
  -min MIN_DISTANCE, --min_distance MIN_DISTANCE
                        Minimal pairwise distance between sequences. If
                        p-distance is lower than min_distance sequence with
                        higher serial number will be removed from the dataset
  -max MAX_DISTANCE, --max_distance MAX_DISTANCE
                        Maximal pairwise distance between sequences. If
                        p-distance is higher than max_distance sequence with
                        higher serial number will be removed from the dataset
```

### *remove_random.py*

Divides all sequences into groups by the first 5 characters in GenBank accession number. 
Plots distribution of groups sizes. Randomly removes k% sequences in groups which size exceed m. k and m are defined by user.

Saves the resulting sequences in fasta-format in the directory of input file


#### Usage
```
remove_random.py [-h] -input INPUT_FILE

  -input INPUT_FILE, --input_file INPUT_FILE
                        Path to input file in fasta-format
```


### *cut_fasta.py*

Finds positions without gaps in reference sequence (the first in alignment), 
removes the columns with gaps in reference sequence from alignment. Reference sequences is supposed to be the first in input-file. 

Saves the resulting alignment in fasta-format in the directory of input file.



#### Usage
```
cut_fasta.py [-h] -input INPUT_FILE

  -input INPUT_FILE, --input_file INPUT_FILE
                        Path to input alignment in fasta-format
```

## Requirements

* python 3
* Biopython
* numpy
* matplotlib

