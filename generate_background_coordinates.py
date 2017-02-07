#!/usr/bin/env python
"""
Given a set of genomic coordinates in BED format:
chr start end
...

calculates a GC matched set of genomic coordinates.
Only first 3 columns of input file will be used- all other columns are ignored.
"""

### imports ###
import os
import sys
import numpy as np
import argparse


def read_target_positions(file_path):
    """
    reads a bed file and returns a list of tuples containing genomic coordinates
    """
    with open(file_path) as f:
        data = f.readlines()
    positions = []
    for line in data:
        tokens = line.strip().split()
        chrom = tokens[0]
        start = int(tokens[1])
        end = int(tokens[2])
        name = tokens[3]
        positions.append((chrom, start, end, name))
    return positions 

def get_random_background(target_positions, 
                        size_ratio = 1.0, 
                        tolerance = 0.01, 
                        N_threshold = 0.5,
                        genome = 'mm10'
                        ):
    """
    target_sequences: 2D numpy array, list of genomic coordinates for target 
                      sequences [[chr,start,end],...]
    size_ratio: float, ratio of background sequences to target sequences
    tolerance: float, max difference in GC content between target and background
    """
    
    ###load genome into memory
    
    # index target positions
    # {chr:[]}, value is chromosome length boolean array
    # largest chromosome has 200 million bps 

    genome_path = os.path.dirname(__file__) + '/mm10/'

    if genome == 'mm10':
        _chromosomes = ['chr1' , 'chr2' , 'chr3' , 'chr4' , 'chr5' , 
                        'chr6' , 'chr7' , 'chr8' , 'chr9' , 'chr10', 
                        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
                        'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
    _chrom_size_dict = {}
    _chrom_seq_dict = {}

    print('reading genome')
    for chrom in _chromosomes:
        with open(genome_path + chrom + '.fa') as f:
            data = f.readlines()
        seq = ''.join(x.upper().strip() for x in data[1:])
        size = len(seq)
        _chrom_size_dict[chrom] = size
        _chrom_seq_dict[chrom] = seq
    _numChromosomes = len(_chromosomes)
    print('done reading genome')
    
    target_chr_position_dict = {x:np.zeros(200000000) for x in _chromosomes} 
    ### initialize target_chr_position_dict using target positions
    ### retreive target sequences
    target_sequences = []
    for pos in target_positions:
        chrom = pos[0]        
        start = int(pos[1])
        end = int(pos[2])
        # use 0 indexing of position, versus 1 indexing used in fasta
        target_chr_position_dict[chrom][start-1:end] = 1         
        seq = _chrom_seq_dict[chrom][start:(end)]
        target_sequences.append(seq)
    ### calculate GC content and average length of the target sequences
    target_gc_count = 0
    target_length_count = 0
    for s in target_sequences:
        target_gc_count += s.count('G')
        target_gc_count += s.count('C')
        target_length_count += len(s)
    # GC content of target sequences
    target_gc_content = target_gc_count/(target_length_count+0.0000001) 
    # average length of target sequences
    mean_target_length = target_length_count/len(target_sequences)     
    mean_target_length = int(np.floor(mean_target_length))
    ### select random genomic loci such that they do no overlap target sequences
    numSelected = 0
    # candidate pool of background seqs is size_ratio X larger
    numToSelect = len(target_positions) * size_ratio 
    candidate_positions = []
    numNallowed = int(N_threshold * mean_target_length) # number of allowable Ns
    counter = 0
    while numSelected < numToSelect:
        if counter % 100000 == 0:
            print(counter, numSelected)
        # select random chromsome
        chromIndex = np.random.randint(_numChromosomes)
        randChrom = _chromosomes[chromIndex]
        randChromSize = _chrom_size_dict[randChrom]
        # must find non overlapping segment on this chromosome before moving on
        selectedSequence = False
        while not selectedSequence:
            counter += 1
            randStart = np.random.randint(randChromSize)
            randEnd = randStart + mean_target_length
            overlap_sum = np.sum(target_chr_position_dict[randChrom][randStart:(randEnd + 1)])
            
            if not overlap_sum > 0:
                randSeq = _chrom_seq_dict[randChrom][randStart:(randEnd+1)]
                numN = randSeq.count('N')
                if numN <= numNallowed:
                    rand_gc_count = randSeq.count('G')+ randSeq.count('C')
                    rand_gc = rand_gc_count/mean_target_length
                    if abs(target_gc_content - rand_gc) <= tolerance:
                        selectedSequence = True
                        numSelected+=1
                        candidate_positions.append([randChrom, randStart, randEnd, randSeq])
    # calcuate GC content of background samples
    background_gc_count = 0
    background_length = 0
    for cp in candidate_positions:
        s = cp[3]
        background_gc_count += s.count('G')
        background_gc_count += s.count('C')
        background_length += len(s)
    background_gc_content = background_gc_count/(background_length+0.0000001)
    print('target GC:', target_gc_content, 
          'background GC:', background_gc_content, 
          'target length:', mean_target_length,
          'numTargetPositions',len(target_positions),
          'backgroundPositions', len(candidate_positions))
    return candidate_positions

def write_background_positions(background_positions, output_dir):
    """
    converts background positions into a bed file and a 
    """

    bed_file = open(output_dir + '/background.bed', 'w')
    fasta_file = open(output_dir + '/background.fasta', 'w')
    counter = 0
    for pos in background_positions:
        chrom = pos[0]
        start = str(pos[1])
        end = str(pos[2])
        seq = str(pos[3])
        randID = 'bg_' + str(np.random.randint(100000)) + '_' + str(counter)
        counter += 1
        bed_file.write('\t'.join([chrom, start, end, randID, '\n']))
        fasta_file.write('>' + randID + '\n')
        fasta_file.write(seq + '\n')
    bed_file.close()
    fasta_file.close()

if __name__ == '__main__':
    input_path = sys.argv[1]
    output_path = sys.argv[2]

    parser = argparse.ArgumentParser(description='Constructs random GC matched'+
                                     'background regions')
    parser.add_argument("samples",
        help="space separated listed of Homer peak files", nargs='+')
    parser.add_argument("outputPath",
        help="directory where output files should be written",
        default="~/", type=str)
    parser.add_argument("-threshold",
        help="idr threshold to use",
        default = "0.05", type=float)
    parser.add_argument("-scoreColumn",
        help="column to use for ranking peaks",
        default = ['findPeaks','Score'],
        type=str,
        nargs='+')
    parser.add_argument("-print", 
        help="just print commands", 
        default = False, action = "store_true")
    # parse arguments
    args = parser.parse_args()

    samples = args.samples
    outPath = args.outputPath
    threshold = args.threshold
    scoreColumn = ' '.join(args.scoreColumn)
    justPrint = False

    target_positions = read_target_positions(input_path)
    
    background_positions = get_random_background(target_positions, 
                                                size_ratio = 1.0, 
                                                tolerance = 0.01, 
                                                N_threshold = 0.5,
                                                genome = 'mm10'
                                                )
    write_background_positions(background_positions, output_path) 
