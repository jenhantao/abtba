### imports ###
import numpy as np
from scipy.stats import pearsonr
import Bio
from Bio import motifs

# scoring matrix
#   A  C  G  T
# A 1 -1 -1 -1
# C-1  1 -1 -1
# G-1 -1  1 -1
# T-1 -1 -1  1
scoringMatrix = np.array([[1.0,-1.0,-1.0,-1.0], [-1.0,1.0,-1.0,-1.0], [-1.0,-1.0,1.0,-1.0], [-1.0,-1.0,-1.0,1.0]]) # favor matches and mismatches equally
gapPenalty = -5.0

def create_logo(motif, output_path, fmt='svg'): 
    motif.weblogo(output_path,  
        format = fmt, 
        show_errorbars = False, 
        show_xaxis= False, 
        show_yaxis= False, 
        show_ends = False, 
        show_fineprint = False, 
        unit_name='probability',
        color_scheme = 'color_classic',
        ) 

# reads all motif files in a directory 
# inputs: path to a directory containing homer motif files
# outputs: an array of tuples representing each motif
def readMotifFile(motifPath):
    with open(motifPath) as f:
        m = motifs.read(f, 'jaspar')
    name = m.name
    matrix = [m.pwm['A'], m.pwm['C'], m.pwm['G'], m.pwm['T']]
    matrix = np.array(matrix).T
        
    return (name,np.array(matrix), m.matrix_id)
# inputs: two motif data objects
# outputs: an alignment data object
def local_align_motifs(motif1, motif2):
    length1 = motif1[1].shape[0]
    length2 = motif2[1].shape[0]

    # initialize alignment score matrix
    scoreMatrix = np.zeros((length1+1,length2+1))

    # populate score matrix
    for i in range(1, length1+1):
        for j in range(1, length2+1):
            matchScore = scoreMatrix[i-1][j-1] + score_column_pearson(motif1[1][i-1], motif2[1][j-1])
            deleteScore = scoreMatrix[i-1][j] + gapPenalty  # penalize deletions
            insertScore = scoreMatrix[i][j-1] + gapPenalty  # penalize insertions
            #scoreMatrix[i][j] = np.max([0,matchScore, deleteScore, insertScore])
            scoreMatrix[i][j] = np.max([matchScore, deleteScore, insertScore])
    # perform traceback
    alignMatrix1 = []
    alignMatrix2 = []
    i = length1
    j = length2
    #print motif1[0], motif1[0]
    #print np.array_str(scoreMatrix, max_line_width=200)
    while i > 0 or j > 0:
    #       print i,j
        if i > 0 and j > 0 and  scoreMatrix[i][j] == scoreMatrix[i-1][j-1] + score_column_pearson(motif1[1][i-1], motif2[1][j-1]):
            alignMatrix1.append(motif1[1][i-1])
            alignMatrix2.append(motif2[1][j-1])
            i -= 1
            j -= 1
        elif i > 0 and scoreMatrix[i][j] == scoreMatrix[i-1][j] + gapPenalty:
            alignMatrix1.append(motif1[1][i-1])
            alignMatrix2.append([0.25,0.25,0.25,0.25])
            #alignMatrix2.append([0.0,0.0,0.0,0.0])
            i -= 1
        elif j > 0 and scoreMatrix[i][j] == scoreMatrix[i][j-1] + gapPenalty:
            alignMatrix1.append([0.25,0.25,0.25,0.25])
            #alignMatrix1.append([0.0,0.0,0.0,0.0])
            alignMatrix2.append(motif2[1][j-1])
            j -= 1
        else:   
            if i > 0 and j > 0:
                alignMatrix1.append(motif1[1][i-1])
                alignMatrix2.append(motif2[1][j-1])
                i -= 1
                j -= 1
            elif i > 0:
                alignMatrix1.append(motif1[1][i-1])
                #alignMatrix2.append([0.25,0.25,0.25,0.25])
                alignMatrix2.append([0.0,0.0,0.0,0.0])
                i -= 1
            elif j > 0:
                #alignMatrix1.append([0.25,0.25,0.25,0.25])
                alignMatrix1.append([0.0,0.0,0.0,0.0])
                alignMatrix2.append(motif2[1][j-1])
                j -= 1
    alignMatrix1 = np.array(alignMatrix1[::-1])
    alignMatrix2 = np.array(alignMatrix2[::-1])
    
    return (alignMatrix1, alignMatrix2), scoreMatrix[-1][-1]

# inputs: motif pwm
# outputs: write position weight matrix file
def writePWMMatrix(matrix, name, outputPath, header_line = None):
    outFile = open(outputPath, "w")
    if header_line == None:
        outFile.write(">"+name+"\t"+name+"\n")
    else:
        outFile.write(header_line)
    for freq in matrix:
        normFreq = freq/np.sum(freq)
        outFile.write("\t".join([str(x) for x in normFreq])+"\n")
    outFile.close()

# given the frequencies at two positions, computes the match score 
def scoreMatch(freq1, freq2): 
    # A, C, G, T 
    score = 0.0 
    for ind1 in range(4): 
        for ind2 in range(4): 
            frac = freq1[ind1] * freq2[ind2]  
            score += frac * scoringMatrix[ind1][ind2]
    return score 

# given the frequencies at two positions, computes the pearson correlation
def score_column_pearson(freq1, freq2):
    pearson, p = pearsonr(freq1, freq2)
    return pearson

# reverse complements a pwm
def revCompMotif(motif):
    name = motif[0]
    scores = []
    for freqs in motif[1]:
        scores.append(freqs[::-1])
    scores = np.array(scores[::-1])
    toReturn = (name, scores)
    return toReturn

# calculates the distance between two motifs
# inputs: two motifs
# outputs: a distance value
def calcDistance(motif1, motif2):
    dist = 0.0
    # assumes that the motifs are the same size - do an alignment to add blanks for motifs of different size
    for i in range(motif1.shape[0]):
        for j in range(len(motif1[0])):
            dist += ((motif1[i][j] - motif2[i][j]) ** 2.0)
    dist = dist ** 0.5
    dist = dist / 4.0
    return dist

# calculates the Pearson correlation between two motifs
# inputs: two motifs
# outputs: a Pearson correlation value
def calcCorrelation(motif1, motif2):
    a = 0.0 
    b = 0.0 
    c = 0.0 
    for i in range(motif1.shape[0]):
        for j in range(len(motif1[0])):
            dif1 = motif1[i][j] - np.mean(motif1[i])
            dif2 = motif2[i][j] - np.mean(motif2[i])
            a += dif1 * dif2
            b += dif1 ** 2.0
            c += dif2 ** 2.0
    toReturn = a/((b*c) ** 0.5)
    return toReturn

# cleans ambigous positions from the ends of a pwm
def cleanMatrix(matrix):
    startIndex = -1
    endIndex = matrix.shape[0] + 1
    for i in range(matrix.shape[0]):
        freqs = matrix[i]
        ambiguous = True
        for freq in freqs:
            if freq > 0.4:
                ambiguous = False
        if ambiguous and (i == startIndex + 1 or i==startIndex):
            startIndex = i
        else:  
            break
    for i in range(matrix.shape[0])[::-1]:
        freqs = matrix[i]
        ambiguous = True
        for freq in freqs:
            if freq > 0.4:
                ambiguous = False
        if ambiguous and (i == endIndex - 1 or i==endIndex):
            endIndex = i
        else:  
            break
    endIndex = np.min([matrix.shape[0] + 1, endIndex]) - 1
    startIndex = np.max([-1, startIndex]) + 1
    return matrix[startIndex:endIndex]


# given two motif objects, aligns the motifs using a needleman-wunsch derivative
# using the pearson correlation to score columns
# inputs: two motif data objects
# outputs: an alignment data object
def global_align_motifs(motif1, motif2):

    length1 = motif1[1].shape[0]    
    length2 = motif2[1].shape[0]    

    # initialize alignment score matrix
    scoreMatrix = np.zeros((length1+1,length2+1))
    for i in range(length1+1):
        scoreMatrix[i][0] = i * -1.0
    for j in range(length2+1):
        scoreMatrix[0][j] = j * -1.0
    
    # populate score matrix
    for i in range(1, length1+1):
        for j in range(1, length2+1):
            matchScore = scoreMatrix[i-1][j-1] + score_column_pearson(motif1[1][i-1], motif2[1][j-1])
            deleteScore = scoreMatrix[i-1][j] + gapPenalty  # penalize deletions
            insertScore = scoreMatrix[i][j-1] + gapPenalty  # penalize insertions
            scoreMatrix[i][j] = np.max([matchScore, deleteScore, insertScore])
    # perform traceback
    alignMatrix1 = []
    alignMatrix2 = []
    i = length1 
    j = length2
    while i > 0 or j > 0:
        if i > 0 and j > 0 and  scoreMatrix[i][j] == scoreMatrix[i-1][j-1] + score_column_pearson(motif1[1][i-1], motif2[1][j-1]):
            alignMatrix1.append(motif1[1][i-1])
            alignMatrix2.append(motif2[1][j-1])
            i -= 1
            j -= 1
        elif i > 0 and scoreMatrix[i][j] == scoreMatrix[i-1][j] + gapPenalty:
            alignMatrix1.append(motif1[1][i-1])
            alignMatrix2.append([0.25,0.25,0.25,0.25])
            i -= 1
        elif j > 0 and scoreMatrix[i][j] == scoreMatrix[i][j-1] + gapPenalty:
            alignMatrix1.append([0.25,0.25,0.25,0.25])
            alignMatrix2.append(motif2[1][j-1])
            j -= 1
        else:
            if i > 0:
                alignMatrix1.append(motif1[1][i-1])
                #alignMatrix2.append([0.25,0.25,0.25,0.25])
                alignMatrix2.append([0.0,0.0,0.0,0.0])
                i -= 1
            elif j > 0:
                alignMatrix2.append(motif2[1][j-1])
                #alignMatrix1.append([0.25,0.25,0.25,0.25])
                alignMatrix1.append([0.0,0.0,0.0,0.0])
                j -= 1
    alignMatrix1 = np.array(alignMatrix1[::-1])
    alignMatrix2 = np.array(alignMatrix2[::-1])
    return (alignMatrix1, alignMatrix2), scoreMatrix[-1][-1]
