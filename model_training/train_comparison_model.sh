#! /bin/bash
################################################################################

### SUMMARY OF FUNCTIONS ###

###


### OPTIONS AND ARGUMENTS ###
# -t generate qsub scripts but do not execute them
### 

### set default options ###

testing=false

# check number of arguments
if [ $# -lt 3 ] 
then
    echo "Usage: "
    echo "train_comparison_model.sh <bed_file_1> <bed_file_2> <genome> \
<output directory> [optional arguments]"
    echo "Options:
-t    generate scripts but do not execute them
"
    exit 1
fi

### parse the input ###

OPTIND=4
while getopts "t" option ; do # set $o to the next passed option
    case "$option" in  
    t)  
        testing=true
    ;;  
    esac
done

bed_file_1=$1
bed_file_2=$2
genome=$3
output_dir=$4
###

script_directory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
motif_directory=${script_directory/model_training/default_motifs}

# if directories don't exist, create them
if [ ! -d $output_dir/ ]; 
    then mkdir -p $output_dir; 
else
    # delete the existing files
    rm $output_dir/*
    # create a script file
    touch $output_dir/run.sh
fi

# generate peak sequence fasta file name
seq_file_1=${bed_file_1##*/}
seq_file_1=${seq_file_1/.bed/.fasta};
seq_file_1=$output_dir/${seq_file_1}

seq_file_2=${bed_file_2##*/}
seq_file_2=${seq_file_2/.bed/.fasta};
seq_file_2=$output_dir/${seq_file_2}

## execute command to extract sequences
echo "$script_directory/extract_sequences.py $bed_file_1 $genome $seq_file_1">> $output_dir/run.sh
echo "$script_directory/extract_sequences.py $bed_file_2 $genome $seq_file_2">> $output_dir/run.sh
#

echo "$script_directory/create_features.py -num_procs 12 $seq_file_1 $seq_file_2  $output_dir $motif_directory/*">> $output_dir/run.sh

# calculate motif scores for peaks and background
combined_features=$output_dir/combined_features.tsv
labels=$output_dir/labels.txt
echo "$script_directory/train_classifier.py $combined_features $labels $output_dir/">> $output_dir/run.sh

# perform insilico mutagenesis
echo "$script_directory/calc_feature_significance.py -num_procs 12 -num_iterations 5 $combined_features $labels $output_dir/">> $output_dir/run.sh

if ! $testing;
then 
bash $output_dir/run.sh
fi

### LICENSE STATEMENT ###
# Copyright (c) 2018, Jenhan Tao 
# All rights reserved. 
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met: 
#
# * Redistributions of source code must retain the above copyright notice, 
#   this list of conditions and the following disclaimer. 
# * Redistributions in binary form must reproduce the above copyright 
#   notice, this list of conditions and the following disclaimer in the 
#   documentation and/or other materials provided with the distribution. 
# * Neither the name of UC San Diego nor the names of its contributors may 
#   be used to endorse or promote products derived from this software 
#   without specific prior written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE. 
###
###############################################################################
