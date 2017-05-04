'''
given a path to a score matrix and a directory containing motifs, produces plots
 showing the effect of merging motifs, creates separate folders for each group 
of motifs
'''

### NOTES ###
## HOMER MUST BE INSTALLED http://homer.salk.edu/homer/ ###

### imports ###
import sys
import numpy as np
import os
from os import listdir
from os.path import isfile, join
from subprocess import call
import argparse
from motif_utilities import *
import shutil



def mergeMotifs(motifArray):
    '''
    given an array of tuples representing motifs, produces a consensus motif
    inputs: an array of motifs
    outputs: motif representing the consensus motif (name, matrix)
    '''
    names = []    
    for i in range(len(motifArray)):
        names.append(motifArray[i][0])
    name = "_".join(sorted(list(set(names))))+"_merged"
    if len(motifArray) < 2:
        return None
    alignedMotifs = []
    alignAgainst = None

    #find the longest motif, move it to the front
    maxLength = -1
    maxLengthMotif = None
    for motif in motifArray:
        if motif[1].shape[0] > maxLength:
            maxLength = motif[1].shape[0]
            maxLengthMotif = motif
    motifArray.remove(maxLengthMotif)
    motifArray.insert(0, maxLengthMotif)

    orientedMotifArray = [ motifArray[0] ] # array of motifs with all motifs in the best orientation
    directions = [orientedMotifArray[0][0]+"_fwd"]

    # determine the orientation of each motif relative to the longest motif
    for motif in motifArray[1:]:
        # calc scores for the forward direction
        alignment_fwd, alignScore_fwd = localAlignMotifs(orientedMotifArray[0], 
                                                         motif)
        r_fwd = calcCorrelation(alignment_fwd[0], alignment_fwd[1])
        # calc scores for one motif reversed
        rcMotif = revCompMotif(motif)
        alignment_rev, alignScore_rev = localAlignMotifs(motifArray[0], rcMotif)
        r_rev = calcCorrelation(alignment_rev[0], alignment_rev[1])
        if r_rev > r_fwd:
            orientedMotifArray.append(rcMotif)
            directions.append(motif[0]+"_rev")
        else:
            orientedMotifArray.append(motif)
            directions.append(motif[0]+"_fwd")

    # iteratively align motifs against the longest motif
    alignAgainst = orientedMotifArray[0][1]
    if len(orientedMotifArray) > 1:
        for i in range(1, len(orientedMotifArray)):
            cleanedMotif = (orientedMotifArray[i][0], cleanMatrix(orientedMotifArray[i][1]))
            alignment, score= localAlignMotifs(("name", alignAgainst), cleanedMotif)
            align1 = alignment[0]
            align2 = alignment[1]
            alignAgainst = (align1 + align2)/2.0
            alignAgainst = cleanMatrix(alignAgainst)
    alignAgainst = cleanMatrix(orientedMotifArray[0][1])
    consensus = (name, alignAgainst)
    orientedMotifArray.insert(0,consensus)
    return orientedMotifArray

def thresholdClusterMotifs(scoreArray, allMotifs, motifNames, outputPath):
    '''
    given a score matrix for an array of motifs, merges motifs and writes a new 
    files for the new set of motifs
    inputs: score matrix, array of motifs, threshold, outputPath
    '''
    mergeDict = {} # key: motif index, value: set of motifs that should be merged together
    # copy heatmap.js file
    heatmap_script_path = os.path.dirname(__file__) + '/heatMap.js'
    shutil.copy(heatmap_script_path, outputPath +'/html_files/heatMap.js')
    # create list page 
    motifListFile = open(outputPath+"/motifList.txt", "w")
    listFile = open(outputPath+"/allList.html", "w")
    listFile.write("<html><head><style>table, th, tr, td {border: 1px solid black;} .nameCol{word-wrap: break-word;max-width: 250px;} table {border-collapse:collapse;}</style></head><body>\n")
    listFile.write('<table><thead><tr><th>Motif Number</th><th>Motif Name</th><th>Full Motif Name</th><th>Logo</th><th>PWM</th></tr></thead><tbody>\n')
    
    # based on table, compute which motifs to merge
    for i in range(scoreArray.shape[0] - 1):
        for j in range(i + 1, scoreArray.shape[0]):
            if scoreArray[i][j] > threshold:
                mergeSet = None
                if i in mergeDict:
                    mergeSet = mergeDict[i]
                elif j in mergeDict:
                    mergeSet = mergeDict[j]    
                else:
                    mergeSet = set()
                mergeSet.add(i)
                mergeSet.add(j)
                
                mergeDict[i] = mergeSet
                mergeDict[j] = mergeSet
    toMergeSets = set(frozenset(i) for i in mergeDict.values())
    unmergedMotifIndices = set(range(scoreArray.shape[0]))
    
    seenNames = set()
    motif_count = 0
    for ms in toMergeSets:
        motif_count+=1
        unmergedMotifIndices -= ms

        mergeNames = []
        toMerge = []    
            
        # get names of motifs being merged
        for ind in ms:
            toMerge.append(allMotifs[ind])
            mergeNames.append(motifNames[ind])

        # create table from merged indices
        mergeNames.sort()
        consensusName = "_".join(sorted(list(set(mergeNames))))+ "_merged"
        if not consensusName in seenNames:
            # don't add repeats
            seenNames.add(consensusName)
            # merge consensus motif
            motifs = mergeMotifs(toMerge) # list of all motifs associated with merge
            consensusMotif = motifs[0]

            # write position weight matrix
            writePWMMatrix(consensusMotif[1], 
                           consensusMotif[0], 
                           outputPath+"/html_files/"+consensusMotif[0]+".motif")
            writePWMMatrix(consensusMotif[1], 
                           consensusMotif[0], 
                           outputPath+"/clustered_motifs/"+consensusMotif[0]+".motif")
            # call homer to create logos
            call(["motif2Logo.pl" ,outputPath+"/html_files/"+consensusMotif[0]+".motif"])

            # create merged motif page
            mergedMotifFile = open(outputPath+"/html_files/"+consensusName+".html", "w")
            mergedMotifFile.write("<html><head><style> td {border: 1px solid black;} .rotate{-webkit-transform:rotate(-90deg); writing-mode: tb-rl;filter: flipv fliph;white-space:nowrap;display:block} table {border-collapse:collapse;}</style><script src='http://code.jquery.com/jquery-2.1.1.min.js'></script><script src='heatMap.js'></script></head><body>\n")
            mergedMotifFile.write("<h1>"+consensusName+"</h2>\n")
            # show logo
            mergedMotifFile.write("<h2>Logo</h2>\n")
            mergedMotifFile.write("<img src = '" + motifs[0][0]+".motif.png'>\n")
            # show pwm
            mergedMotifFile.write("<h2>Position Weight Matrix</h2>\n")
            mergedMotifFile.write("<table><thead><tr><th>Position</th><th>A</th><th>C</th><th>G</th><th>T</th></tr></thead>\n<tbody>\n")
            for i in range(motifs[0][1].shape[0]):
                mergedMotifFile.write("<tr><td>"+str(i+1)+"</td>")
                for j in range(motifs[0][1].shape[1]):
                    mergedMotifFile.write("<td>"+str(np.round(motifs[0][1][i][j], 3))+"</td>")
                mergedMotifFile.write("</tr>\n")
            mergedMotifFile.write("</tbody></table>\n")
                    
            # show download link for pwm
            mergedMotifFile.write("<a href='"+consensusName+".motif'>Download Position Weight Matrix</a>")
            
            # list merged motifs in table
            mergedMotifFile.write("<h2>Contributing Motifs</h2>\n")    
            mergedMotifFile.write("<table><thead><tr><th>Motif Name</th><th>Full Motif Name</th><th>Logo</th><th>PWM</th></tr></thead><tbody>\n")
            for motif in motifs[1:]:
                # write position weight matrix
                writePWMMatrix(motif[1], motif[0], outputPath+"/html_files/"+motif[0]+".motif")
                # call homer to create logos
                call(["motif2Logo.pl" ,outputPath+"/html_files/"+motif[0]+".motif"])

                mergedMotifFile.write("<tr><td><a href='"+motif[0]+".html'>" +motif[0]+"</a></td><td>"+motif[0]+"</td><td><img src = '" + motif[0]+".motif.png'></td><td><a href='"+motif[0]+".motif' target='_blank'>Download</a></tr>\n")
            mergedMotifFile.write("</tbody></table>\n")
            # find related motifs
            mergedMotifFile.write("<h2>Related Motifs</h2>\n")
            relatedScores = scoreArray[list(ms)[0]]
            rankings = sorted(range(len(relatedScores)), key=lambda x: relatedScores[x])
            start = np.min([len(ms), len(relatedScores)])
            relatedIndices = rankings[start: np.min([start+10, len(relatedScores)])]
            mergedMotifFile.write("<br><br><br><br>\n")
            mergedMotifFile.write("<table class='heat-map'><thead><tr><th></th>")
            for ri in relatedIndices:
                mergedMotifFile.write("<th><span class='rotate'><a href='"+allMotifs[ri][0]+".html'>"+allMotifs[ri][0]+"</a></span></th>")
            mergedMotifFile.write("</tr></thead>\n<tbody>")
            for i in range(len(relatedIndices)):
                mergedMotifFile.write("<tr class='stats-row'><th class='stats-title'><a href='"+allMotifs[relatedIndices[i]][0]+".html'>"+allMotifs[relatedIndices[i]][0]+"</a></th>")
                for j in range(len(relatedIndices)):
                    mergedMotifFile.write("<td>"+str(int(np.round(scoreArray[relatedIndices[i]][relatedIndices[j]]*100,3)))+"</td>")
                mergedMotifFile.write("</tr>\n")
            mergedMotifFile.write("</tbody></table>")

            mergedMotifFile.write("<script>$(function(){$('table th').height('10px');$('.rotate').width('10px');});</script>")
            mergedMotifFile.write("</body></html>")
            mergedMotifFile.close()
            # add merged motif to list page
            listFile.write("<tr><td>"+str(motif_count) + "</td><td class='nameCol'><a href='html_files/"+consensusName+".html'>" +motifs[0][0]+"</a></td><td class='nameCol'>"+motifs[0][0]+"</td><td><img src = 'html_files/" + motifs[0][0]+".motif.png'></td><td><a href='html_files/"+motifs[0][0]+".motif' target='_blank'>Download</a></td></tr>\n")
            motifListFile.write(motifs[0][0]+".motif\n")

    # add unmerged motifs to list file

    # sort unmerged indices
    unmergedMotifIndexTuples = tuple(zip([allMotifs[x][0] for x in unmergedMotifIndices], unmergedMotifIndices))
    unmergedMotifIndexTuples = sorted(unmergedMotifIndexTuples)
    unmergedMotifIndices_sorted = [x[1] for x in unmergedMotifIndexTuples]
    
    for ind in unmergedMotifIndices_sorted:
        motif_count+=1
        listFile.write("<tr><td>" + str(motif_count) + "</td><td class='nameCol'><a href='html_files/"+allMotifs[ind][0]+".html'>" +allMotifs[ind][0]+"</a></td><td class='nameCol'>"+allMotifs[ind][0]+"</td><td><img src = 'html_files/" + allMotifs[ind][0]+".motif.png'></td><td><a href='html_files/"+allMotifs[ind][0]+".motif' target='_blank'>Download</a></td></tr>\n")
        motifListFile.write(allMotifs[ind][0]+".motif\n")
        writePWMMatrix(allMotifs[ind][1], 
                       allMotifs[ind][0], 
                       outputPath+"/clustered_motifs/"+allMotifs[ind][0]+".motif")
        

    # write files for all individual motifs
    for ind in range(len(allMotifs)):
        # write pwm matrix file
        writePWMMatrix(allMotifs[ind][1], allMotifs[ind][0], outputPath+"/html_files/"+allMotifs[ind][0]+".motif")
        # call homer to create logos
        call(["motif2Logo.pl" ,outputPath+"/html_files/"+allMotifs[ind][0]+".motif"])
        # write html file
        indMotifFile = open(outputPath+"/html_files/"+allMotifs[ind][0]+".html", "w")
        indMotifFile.write("<html><head><style> td {border: 1px solid black;} .rotate{-webkit-transform:rotate(-90deg); writing-mode: tb-rl;filter: flipv fliph;white-space:nowrap;display:block} table {border-collapse:collapse;}</style><script src='http://code.jquery.com/jquery-2.1.1.min.js'></script><script src='html_files/heatMap.js'></script></head><body>\n")
        indMotifFile.write("<h1>"+allMotifs[ind][0]+"</h1>\n")
        # show logo
        indMotifFile.write("<h2>Logo</h2>\n")
        indMotifFile.write("<img src = '" + allMotifs[ind][0]+".motif.png'>\n")
        # show pwm
        indMotifFile.write("<h2>Position Weight Matrix</h2>\n")
        indMotifFile.write("<table><thead><tr><th>Position</th><th>A</th><th>C</th><th>G</th><th>T</th></tr></thead>\n<tbody>\n")
        for i in range(allMotifs[ind][1].shape[0]):
            indMotifFile.write("<tr><td>"+str(i+1)+"</td>")
            for j in range(allMotifs[ind][1].shape[1]):
                indMotifFile.write("<td>"+str(np.round(allMotifs[ind][1][i][j], 3))+"</td>")
            indMotifFile.write("</tr>\n")
        indMotifFile.write("</tbody></table>\n")
        indMotifFile.write("<a href='"+allMotifs[ind][0]+".motif'>Download Position Weight Matrix</a>")

        # find related motifs
        indMotifFile.write("<h2>Related Motifs</h2>\n")
        relatedScores = scoreArray[ind]
        rankings = sorted(range(len(relatedScores)), key=lambda x: relatedScores[x], reverse=True)
        start = 0
        relatedIndices = rankings[start: np.min([start+10, len(relatedScores)])]
        indMotifFile.write("<br><br><br><br>\n")
        indMotifFile.write("<table class='heat-map'><thead><tr><th></th>")
        for ri in relatedIndices:
            indMotifFile.write("<th><span class='rotate'><a href='"+allMotifs[ri][0]+".html'>"+allMotifs[ri][0]+"</a></span></th>")
        indMotifFile.write("</tr></thead>\n<tbody>")
        for i in range(len(relatedIndices)):
            indMotifFile.write("<tr class='stats-row'><th class='stats-title'><a href='"+allMotifs[relatedIndices[i]][0]+".html'>"+allMotifs[relatedIndices[i]][0]+"</a></th>")
            for j in range(len(relatedIndices)):
                indMotifFile.write("<td>"+str(int(np.round(scoreArray[relatedIndices[i]][relatedIndices[j]]*100,3)))+"</td>")
            indMotifFile.write("</tr>\n")
        indMotifFile.write("</tbody></table>")
        # add script to fix heights
        indMotifFile.write("<script>$(function(){$('table th').height('10px');$('.rotate').width('10px');});</script>")
        indMotifFile.close()

    # create index file
    scoreFile = open(outputPath+"/allScores.html", "w")
    scoreTsvFile = open(outputPath+"/allScores.tsv", "w")
    scoreFile.write("<html><head><style> td {border: 1px solid black;} .rotate{-webkit-transform:rotate(-90deg); writing-mode: tb-rl;filter: flipv fliph;white-space:nowrap;display:block} table {border-collapse:collapse;}</style><script src='http://code.jquery.com/jquery-2.1.1.min.js'></script><script src='heatMap.js'></script></head><body>\n")
    # add in some blank spaces
    scoreFile.write("<br><br><br><br>\n")
    # write score array as a matrix
    scoreFile.write("<table class='heat-map'><thead><tr><td></td>")
    for name in motifNames:
        scoreFile.write("<th class='stats-title'><span class='rotate'><a href='html_files/"+name+".html'>" + name + "</a></span></th>")
        scoreTsvFile.write("\t"+name)
    scoreTsvFile.write("\n")
    scoreFile.write("</tr></thead>")
    scoreFile.write("<tbody>\n")
    for i in range(scoreArray.shape[0]):
        scoreFile.write("<tr class='stats-row'><th class='stats-title'><a href='html_files/"+motifNames[i]+".html'>" + motifNames[i] + "</a></th>")
        scoreTsvFile.write(motifNames[i])
        for j in range(scoreArray.shape[0]):
            scoreFile.write("<td>" + str(int(np.round(scoreArray[i][j]*100, 3))) + "</td>")
            scoreTsvFile.write("\t"+str(scoreArray[i][j]))
        scoreFile.write("</tr>\n")
        scoreTsvFile.write("\n")
    scoreFile.write("</tbody></table>\n")
    # add script to fix heights
    scoreFile.write("<script>$(function(){$('table th').height('10px');$('.rotate').width('10px');});</script>")
    scoreFile.write("</body></html>")
    listFile.write("</tbody></table>\n")
    listFile.write("</body></html>")
    motifListFile.close()
    listFile.close()
    scoreFile.close()
    scoreTsvFile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='using scores calculated by \
                                     the scoreMotifs.py script, clusters merges\
                                      similar motifs and creates an html \
                                     representation')

    parser.add_argument("scorePath",
        help="path to a npz file containing motif similarity scores",
        type = str)
    parser.add_argument("outputPath",
        help="path to directory where output will be written",
        type=str)
    parser.add_argument("threshold",
        help="threshold for clustering motifs",
        default = 0.8,
        type=float)
    parser.add_argument("motifFiles",
        help="list of moti files to cluster",
        type=str,
        nargs="+")

    # parse arguments
    args = parser.parse_args()

    scorePath = args.scorePath
    outputPath = args.outputPath
    threshold = args.threshold
    motifFiles = args.motifFiles

    if not os.path.isdir(outputPath):
        os.mkdir(outputPath)

    if not os.path.isdir(outputPath + '/html_files/'):
        os.mkdir(outputPath + '/html_files')
    if not os.path.isdir(outputPath + '/clustered_motifs/'):
        os.mkdir(outputPath + '/clustered_motifs')

    # read in motifs
    # find all motifs in input directory

    allMotifs = []
    motifNames = []
    for mf in sorted(motifFiles):
        motif = readMotifFile(mf) # (name, PWM)
        allMotifs.append(motif)
        motif_name = motif[0]
        motifNames.append(motif_name)


    # read in scores
    scoreArray = np.load(scorePath)['arr_0']

    # if output directory doesn't exist, create it
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

    thresholdClusterMotifs(scoreArray, allMotifs, motifNames, outputPath)
