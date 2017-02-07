#!/usr/bin/env python
"""
Given a set of genomic coordinates in BED format:
chr start end
...

calculates a GC matched set of genomic coordinates.
Only the first 3 columns of input file will be used - all other columns are ignored.
"""

### imports ###

def read_target_positions()
    """
    """
    return None

def getRandomBackground(target_positions, 
                        size_ratio = 1.0, 
                        tolerance = 0.01, 
                        N_threshold = 0.5 ):
    """
    target_sequences: 2D numpy array, list of genomic coordinates for target sequences [[chr,start,end],...]
    size_ratio: float, number of background sequences relative to target sequences
    tolerance: float, max difference in GC content between True and background labelled samples
    """
    
    ###load mm10 genome into memory
    
    # index target positions
    # {chr:[]}, value is chromosome length boolean array
    # largest chromosome has 200 million bps 
    _chromosomes = ['chr1' , 'chr2' , 'chr3' , 'chr4' , 'chr5' , 
                    'chr6' , 'chr7' , 'chr8' , 'chr9' , 'chr10', 
                    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
                    'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
    _chrom_size_dict = {}
    _chrom_seq_dict = {}
    for chrom in _chromosomes:
        with open('./mm10_genome/' + chrom + '.fa') as f:
            data = f.readlines()
        seq = ''.join(x.upper().strip() for x in data[1:])
        size = len(seq)
        _chrom_size_dict[chrom] = size
        _chrom_seq_dict[chrom] = seq
    _numChromosomes = len(_chromosomes)
    
    target_chr_position_dict = {x:np.zeros(200000000) for x in _chromosomes} 
    ### initialize target_chr_position_dict using target positions
    ### retreive target sequences
    target_sequences = []
    for pos in target_positions:
        chrom = pos[0]        
        start = pos[1]
        end = pos[2]
        target_chr_position_dict[chrom][start-1:end] = 1 # use 0 indexing of position, versus 1 indexing used in fasta
        seq = _chrom_seq_dict[chrom][start:(end)]
        target_sequences.append(seq)
    ### calculate GC content and average length of the target sequences
    target_gc_count = 0
    target_length_count = 0
    for s in target_sequences:
        target_gc_count += s.count('G')
        target_gc_count += s.count('C')
        target_length_count += len(s)
    target_gc_content = target_gc_count/(target_length_count+0.0000001) # GC content of target sequences
    mean_target_length = target_length_count/len(target_sequences) # average length of target sequences
    mean_target_length = int(np.floor(mean_target_length))
    ### select random genomic loci such that they do no overlap target sequences
    numSelected = 0
    numToSelect = len(target_positions) * size_ratio # candidate pool of background seqs is size_ratio X larger
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
