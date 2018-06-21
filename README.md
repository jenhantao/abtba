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
TBA uses a programatically curated library of motifs to reduce the effects of multiple collinearity, which can be problematic for machine learning models. You can view and download the motif library here: http://homer.ucsd.edu/jtao/merged_motifs/allList.html

## Installing TBA

## Usage

## Interpreting Results

## Visualizing Results

## Sample Data

## Custom Parameters

## Contributors
TBA was created by Jenhan Tao with feedback from Gregory Fonseca and Christopher Benner.

## Citing TBA
If you use TBA for research

## License
Software is subject to a non-exclusive, revocable, non-transferable, and limited right to use the Code for the exclusive purpose of undertaking academic, governmental, or not-for-profit research. Use of the Code or any part thereof for commercial or clinical purposes is strictly prohibited in the absence of a Commercial License Agreement from Deep Genomics. 4
