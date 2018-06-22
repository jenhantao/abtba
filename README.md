# TBA (a Transcription factor Binding Analysis)

## TBA Overview
TBA is a multi-functional machine learning tool for identifying transcription factors associated with genomic features. Specifically, TBA can be applied to:
* ChIP-seq targeting a transcription to identify collaborative binding partners for a given transcription factor. 
* DNAse-seq and ATAC-seq to identify transcription factors associated with open chromatin
* GRO-seq and other assays measuring enhancer activity to identify transcription factors associated with enhancer activty
* Predict the effect of genetic variation in any of the above contexts.

## TBA Algorithm
TBA takes the genomic sequence of sites of interest as input and selects a set of GC-matched background loci. For each locus of interest and background locus, TBA calculates the best match to hundreds of DNA binding motifs, and quantifies the quality of the match as the motif score (aka log likelihood ratio score). To allow for degenerate motifs, all motif matches scoring over zero are considered. The motif scores are then used to train the TBA model to distinguish loci of interest from background loci. TBA scores the probability of observing binding at a sequence by computing a weighted sum over all the motif scores for that sequence. The weight for each motif is learned by iteratively modifying the weights until the model’s ability to differentiate binding sites from background loci no longer improves. The final motif weight measures whether the presence of a motif is correlated with TF binding. 

## Motif Library
TBA uses a programatically curated library of motifs to reduce the effects of multiple collinearity, which can be problematic for machine learning models. You can view and download the motifs at: homer.ucsd.edu/jtao/merged_motifs/allList.html](http://homer.ucsd.edu/jtao/merged_motifs/allList.html "Motif Library")

## Installing TBA
Content coming soon - please refer to our BioRxiv manuscript for now.

## Usage
TBA is accessible as a series of command line (aka termina) scripts. The easiest way to train a TBA model is to use the train_model_default.sh command. train_model_default.sh will run all TBA commands necessary to train a TBA with default parameters. You can invoke the command like this:
```bash
train_model_default.sh mouse_pu1_peaks.bed mm10 path_to_output 
```
The script will create a script at path_to_output/run.sh and execute it. run.sh will have correctly formatted TBA commands for each step. Output files will be created at path_to_output

## Interpreting Results
Content coming soon - please refer to our BioRxiv manuscript for now.

## Visualizing Results
Content coming soon - please refer to our BioRxiv manuscript for now.

## Sample Data
Content coming soon - please refer to our BioRxiv manuscript for now.

## Example Analysis
Content coming soon - please refer to our BioRxiv manuscript for now.

## TBA Parameters
If you're ever unsure how to use a TBA command simply run the command without any parameters and help text should be displayed. All TBA commands and their associated parameters are listed here. Optional parameters are indicated with the default parameters

**annotate_results_with_genes.py <result_path> <output_path>**

Given a TBA coefficients file or a TBA significance file, maps the motif names to gene names
```
arguments:
* result_path - path to a TBA coefficients or significance file
* output_path - file path where output should be writtente_results_with_genes.py <result_path>
```
**calc_feature_significance.py <feature_path> <label_path> <output_path> -num_iterations 5 -test_fraction 0.2 -num_procs 4**

Performs an in silico mutagenesis test to assign significance to each motif
```
arguments:
* feature_path - path to a standardized feature set created by create_features.py
* label_path - path to a fasta_file containing negative sequences to score
* output_path - directory where output file should be written
optional arguments:
* -num_iterations - number of iterations to train classifier
* -test_fraction - fraction of data to use for testing classifier
* -num_procs - number of processors to use
```
**calculate_all_motif_scores.py <fasta_path> <output_path> <motif_file_1> <motif_file_2>... -num_procs 4 -psuedocount 0.01**

given a fasta file and a list and a list of motif files, calculates all matches to each motif for each sequence
```
arguments:
* fasta_path - path to a fasta_file containing sequences to score
* output_path - directory where output file should be written
* motif_files - list of motif files

optional arguments:
* -num_procs - number of processor cores to use
* -pseudocount - pseudocount for calculating motif scores
```
**usage: create_features.py <positive_sequences_path> <negative_sequences_path> <output_path> <motif_file_1> <motif_file_2>... -num_procs 4 -pseudocount 0.01**

Given a set of negative sequences and positive sequences (in FASTA format) as
well as a list of motifs, calculates motif scores and sequence labels suitable
for training classifier
```
arguments:
* positive_sequences_path - path to a fasta_file containing positive sequences to score
* negative_sequences_path - path to a fasta_file containing negative sequences to score
* output_path - directory where output file should be written
* motif_files - list of motif files

optional arguments:
* -num_procs - number of processor cores to use
* -pseudocount - pseudocount for calculating motif scores
```                     
**extract_sequences.py <bed_file> <genome> <output_file_path>**

Extracts the genomic sequence at positions indicated by input bed file
```
arguments:
* bed_file - path to a bed file
* genome - build of genome to use (eg. hg38, mm10)
* output_file_path - directory where output file should be written
```
**generate_background_coordinates.py <inputPath> <genome> <outputPath> -sizeRatio 1 -numBins 10 -nThreshold 0.1 -filterChromosomes chrM chrY**

Constructs random GC matched background regions
```
arguments:
* inputPath - path to a bed file containing a chr, start, end, and strand column
* genome - genome from which to construct background regions
* outputPath - directory where output files should be written

optional arguments:
* -sizeRatio - size of the background region with respect to the target region
* -numBins - number of bins to use for GC normalization
* -nThreshold - maximum fraction of background sequences that can be N
* -filterChromosomes - chromosomes to ignore
```
**train_classifier.py  <feature_path> <label_path> <output_path> -num_iterations 5-test_fraction 0.2**
                          
Given standardized motif features and matching labels trains a classifier and returns performance metrics and model coefficients
```
arguments:
* feature_path - path to a standardized feature set created by create_features.py
* label_path - path to a fasta_file containing negative sequences to score
* output_path - directory where output file should be written

optional arguments:
* -num_iterations - number of iterations to train classifier
* -test_fraction - fraction of data to use for testing classifier
```
**train_model_default.sh <bed_file> <genome> <output directory>**

Wrapper script for training a TBA model to distinguishes sequences in input bed file from genomic background with default parameters
```
arguments:
* bed_file - input bed file
* genome - genome build
* output_directory - directory where file outputs should be written

optional arguments:
* -t - generate scripts but do not execute them
```
**train_comparison_model.sh <bed_file_1> <bed_file_2> <genome> <output directory>**

Wrapper script for training a TBA model to distinguishes sequences in bed_file_1 from bed_file_2 from genomic background with default parameters
```
arguments:
* bed_file_1 - input bed file
* bed_file_2 - input bed file
* genome - genome build
* output_directory - directory where file outputs should be written
  
optional arguments:
* -t - generate scripts but do not execute them
```
## Authors
TBA was created by Jenhan Tao with feedback from Gregory Fonseca and Christopher Benner. If you have any questions, please send an email to jenhantao@gmail.com. We would be glad to help you apply TBA to your research problem

## Citing TBA
If you use TBA for research, please cite the following manuscript:
```
Fonseca, G. J., Tao, J. Diverse motif ensembles specify non-redundant DNA binding activities of AP-1 family members in macrophages. Preprint at: https://www.biorxiv.org/content/early/2018/06/13/345835 (2018). 
```
## License
```
MIT License

Copyright (c) 2017 Jenhan Tao

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
