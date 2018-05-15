#!/usr/bin/env python
'''
given a path to a score matrix and a directory containing motifs, produces plots
 showing the effect of merging motifs, creates separate folders for each group 
of motifs
'''

### NOTES ###

### imports ###
import sys
import numpy as np
import os
from os import listdir
from os.path import isfile, join
import argparse
from motif_utilities import *
import shutil
import inspect
import Bio
import pandas as pd
import scipy.cluster
from itertools import combinations

def mergeMotifs(motifArray):
    '''
    given an array of tuples representing motifs, produces a consensus motif
    inputs: an array of motifs
    outputs: motif representing the consensus motif (name, matrix)
    '''
    if len(motifArray) < 2:
        return None

    # determine the orientation of each motif relative to the longest motif
    mergedMotifMatrix = motifArray[0][1]
    for motif in motifArray[1:]:
        # align motifs in both orientations
        alignment_fwd, alignScore_fwd = local_align_motifs(('merged_name', mergedMotifMatrix, 'merged_id'), 
                                                            motif)
        rev_comp_motif = revCompMotif(motif)
        alignment_rev, alignScore_rev = local_align_motifs(('merged_name', mergedMotifMatrix, 'merged_id'), 
                                                            rev_comp_motif)

        mergedMotif_alignment_fwd = alignment_fwd[0]
        mergedMotif_alignment_rev = alignment_rev[0]
        compared_motif_alignment_fwd = alignment_fwd[1]
        compared_motif_alignment_rev = alignment_rev[1]

        # calc scores for aligned motifs in both orientations
        pearson_fwd = calcCorrelation(mergedMotif_alignment_fwd, 
                                      compared_motif_alignment_fwd)
        pearson_rev = calcCorrelation(mergedMotif_alignment_rev, 
                                      compared_motif_alignment_rev)
        if pearson_rev > pearson_fwd:
            mergedMotifMatrix = mergedMotif_alignment_rev + compared_motif_alignment_rev
        else:
            mergedMotifMatrix = mergedMotif_alignment_fwd + compared_motif_alignment_fwd
        mergedMotifMatrix = mergedMotifMatrix/2.0
        mergedMotifMatrix = cleanMatrix(mergedMotifMatrix)

    names = []    
    aligned_motif_array = motifArray
    for i in range(len(motifArray)):
        mn = motifArray[i][0]
        names.append(mn)
    name = "_".join(sorted(list(set(names)))[:10])+"_merged"

    consensus = (name, mergedMotifMatrix, name)
    aligned_motif_array.insert(0,consensus)

    return aligned_motif_array

def thresholdClusterMotifs(scoreArray, 
    threshold, 
    allMotifs, 
    motifNames, 
    outputPath,
    create_html=False):
    '''
    given a score matrix for an array of motifs, merges motifs and writes a new 
    files for the new set of motifs
    inputs: score matrix, array of motifs, threshold, outputPath
    '''

    metadata_path= os.path.dirname(
        os.path.abspath(
        inspect.getfile(
        inspect.currentframe()))).replace('motif_tools', 'motif_metadata.tsv')

    metadata_frame = pd.read_csv(metadata_path, sep='\t')
    nameRoot_count_dict = {}
    motifName_family_dict = dict(zip(metadata_frame['Name'].values, 
        metadata_frame['Family'].values))
    motifName_gene_dict = dict(zip(metadata_frame['Name'].values, 
        metadata_frame['Gene'].values))

    alphabet = Bio.Seq.IUPAC.Alphabet.IUPAC.IUPACUnambiguousDNA()

    mergeDict = {} # key: motif index, value: set of motifs that should be merged together
    # create list page 
    motifListFile = open(outputPath+"/motif_list.txt", "w")
    mergedMetadataFile = open(outputPath + '/motif_metadata.tsv', 'w')
    if create_html:
        # copy heatmap.js file
        heatmap_script_path = os.path.dirname(__file__) + '/heatMap.js'
        shutil.copy(heatmap_script_path, outputPath +'/html_files/heatMap.js')

        listFileLines = []
        listFile = open(outputPath+"/allList.html", "w")
        listFile.write("<html><head><style>table, th, tr, td {border: 1px solid black;} .nameCol{word-wrap: break-word;max-width: 250px;} table {border-collapse:collapse;}</style></head><body>\n")
        listFile.write('<table><thead><tr><th>Motif Number</th><th>Motif Name</th><th>Full Motif Name</th><th>Logo</th><th>PWM</th></tr></thead><tbody>\n')
    
    # based on table, compute which motifs to merge
    scoreArray = np.clip(scoreArray, -1, 1)
    dissimilarity = 1-np.abs(scoreArray)
    coords = list(combinations(range(len(motifNames)),2))
    dissimilarity_as_pdist = [dissimilarity[x[0]][x[1]] for x in coords]

    Z=scipy.cluster.hierarchy.linkage(dissimilarity_as_pdist, 
        method = 'centroid')

    tree_cut = scipy.cluster.hierarchy.cut_tree(Z, height=0.1).flatten()
    cluster_set_dict = {}
    for i in range(len(tree_cut)):
        cluster = tree_cut[i]
        if cluster in cluster_set_dict:
            cluster_set_dict[cluster].add(i)
        else:
            cluster_set_dict[cluster] = set([i])
    toMergeSets = [x for x in cluster_set_dict.values() if len(x) > 1]

    unmergedMotifIndices = set(range(scoreArray.shape[0])) # initialize unmerged indices
    seenNames = set()
    motif_count = 0
    for ms in toMergeSets:
        motif_count+=1
        unmergedMotifIndices -= ms

        mergeNames = []
        geneNames = []
        toMerge = []    
            
        # get names of motifs being merged
        for ind in ms:
            toMerge.append(allMotifs[ind])
            mergeNames.append(motifNames[ind])

            if allMotifs[ind][0] in motifName_gene_dict:
                gene = motifName_gene_dict[allMotifs[ind][0]]
            else:
                gene = 'Unknown'
            geneNames.append(gene)
            
        # merged motif doesn't have JASPAR id - use space to track merging instead
        consensus_id_string = '|'.join(mergeNames)
        gene_string = '|'.join(geneNames)

        # create table from merged indices
        mergeNames.sort()

        consensusFamily = 'Unknown'
        for tm in toMerge:
            tm_name = tm[0]
            if tm_name in motifName_family_dict:
                consensusFamily = motifName_family_dict[tm_name]
                if not consensusFamily == 'Unknown':
                    break
        consensusNameRoot = consensusFamily

        
        if consensusNameRoot in nameRoot_count_dict:
            nameRoot_count_dict[consensusNameRoot] += 1
        else:
            nameRoot_count_dict[consensusNameRoot] = 1

        # get the TF family of the first motif
        consensusName = consensusNameRoot+ '_' + str(nameRoot_count_dict[consensusFamily]) + '_merged'
        consensusName = consensusName.replace('/','').replace(' ', '_').replace('(','_').replace(')', '_').replace('__', '_')
        
        if not consensusName in seenNames:
            # don't add repeats
            seenNames.add(consensusName)
            # merge consensus motif
            motifs = mergeMotifs(toMerge) # list of all motifs associated with merge
            consensusMotif = (consensusName, motifs[0][1])

            # write position weight matrix
            counts_dict = {x[0]:x[1] for x in zip(list('ACGT'), 
                consensusMotif[1].T)}
            bio_motif = Bio.motifs.jaspar.Motif(alphabet = alphabet, 
                counts = counts_dict, 
                name = consensusMotif[0],
                matrix_id = consensus_id_string)

            pwm_file = open(outputPath+"/clustered_motifs/"+consensusName+".motif", "w")
            pwm_file.write(Bio.motifs.write([bio_motif], format='jaspar'))
            pwm_file.close()

            # call weblogo to create logos
            if create_html:
                create_logo(bio_motif, 
                    outputPath+'/html_files/'+consensusName+'.motif.svg')

                pwm_file = open(outputPath+"/html_files/"+consensusName+".motif", "w")
                pwm_file.write(Bio.motifs.write([bio_motif], format='jaspar'))
                pwm_file.close()

                # create merged motif page
                mergedMotifFile = open(outputPath+"/html_files/"+consensusName+".html", "w")
                mergedMotifFile.write("<html><head><style> td {border: 1px solid black;} .rotate{-webkit-transform:rotate(-90deg); writing-mode: tb-rl;filter: flipv fliph;white-space:nowrap;display:block} table {border-collapse:collapse;}</style><script src='http://code.jquery.com/jquery-2.1.1.min.js'></script><script src='heatMap.js'></script></head><body>\n")
                mergedMotifFile.write("<h1>"+consensusName+"</h2>\n")
                # show logo
                mergedMotifFile.write("<h2>Logo</h2>\n")
                mergedMotifFile.write("<img width='500px' src = '" + consensusName +".motif.svg'>\n")
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
                    counts_dict = {x[0]:x[1] for x in zip(list('ACGT'),
                        motif[1].T)}
                    bio_motif = Bio.motifs.jaspar.Motif(alphabet = alphabet, 
                        counts = counts_dict, 
                        name = motif[0],
                        matrix_id = motif[2]
                    ) 

                    pwm_file = open(outputPath+"/html_files/"+bio_motif.name+".motif", "w")
                    pwm_file.write(Bio.motifs.write([bio_motif], format='jaspar'))
                    pwm_file.close()

                    # call weblogo to create logos
                    create_logo(bio_motif,
                        outputPath+'/html_files/'+bio_motif.name+'.motif.svg')


                    mergedMotifFile.write("<tr><td><a href='"+motif[0]+".html'>" +motif[0]+"</a></td><td>"+motif[0]+"</td><td><img src = '" + motif[0]+".motif.svg'></td><td><a href='"+motif[0]+".motif' target='_blank'>Download</a></tr>\n")
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
                listFileLines.append((consensusName, "<tr><td>MOTIF_COUNT</td><td class='nameCol'><a href='html_files/"+consensusName+".html'>" +consensusName+"</a></td><td class='nameCol'>"+consensusName+"</td><td><img src = 'html_files/" + consensusName +".motif.svg'></td><td><a href='html_files/"+consensusName+".motif' target='_blank'>Download</a></td></tr>\n"))

            motifListFile.write(consensusName + '\t' + consensus_id_string + '\n')
            mergedMetadataFile.write('\t'.join([consensusName, gene_string, consensusFamily])+'\n')

    # add unmerged motifs to list file

    # sort unmerged indices
    unmergedMotifIndexTuples = tuple(zip([allMotifs[x][0] for x in unmergedMotifIndices], unmergedMotifIndices))
    unmergedMotifIndexTuples = sorted(unmergedMotifIndexTuples)
    unmergedMotifIndices_sorted = [x[1] for x in unmergedMotifIndexTuples]
    
    for ind in unmergedMotifIndices_sorted:
        motif_count+=1
        motif_name = allMotifs[ind][0]
        geneName = 'Unknown'
        family = 'Unknown'
        if motif_name in motifName_gene_dict:
            geneName = motifName_gene_dict[motif_name]
        if motif_name in motifName_family_dict:
            family = motifName_family_dict[motif_name]

        motifListFile.write(allMotifs[ind][0] + '\t' + allMotifs[ind][0] +'\n' )
        mergedMetadataFile.write('\t'.join([motif_name, geneName, family])+ '\n')
        if create_html:
            listFileLines.append((allMotifs[ind][0], "<tr><td>MOTIF_COUNT</td><td class='nameCol'><a href='html_files/"+allMotifs[ind][0]+".html'>" +allMotifs[ind][0]+"</a></td><td class='nameCol'>"+allMotifs[ind][0]+"</td><td><img src = 'html_files/" + allMotifs[ind][0]+".motif.svg'></td><td><a href='html_files/"+allMotifs[ind][0]+".motif' target='_blank'>Download</a></td></tr>\n"))
        counts_dict = {x[0]:x[1] for x in zip(list('ACGT'),
            cleanMatrix(allMotifs[ind][1]).T)}
        bio_motif = Bio.motifs.jaspar.Motif(alphabet = alphabet, 
            counts = counts_dict, 
            name = allMotifs[ind][0],
            matrix_id = allMotifs[ind][2]
            )
        pwm_file = open(outputPath+"/clustered_motifs/"+allMotifs[ind][0]+".motif", 'w')
        pwm_file.write(Bio.motifs.write([bio_motif], format='jaspar'))
        pwm_file.close()
        

    # write files for all individual motifs
    for ind in range(len(allMotifs)):

        # write pwm matrix file
        counts_dict = {x[0]:x[1] for x in zip(list('ACGT'),
            allMotifs[ind][1].T)}
        bio_motif = Bio.motifs.jaspar.Motif(alphabet = alphabet, 
            counts = counts_dict, 
            name = allMotifs[ind][0],
            matrix_id = allMotifs[ind][2]
            )
        
        if create_html:
            pwm_file = open(outputPath+"/html_files/"+bio_motif.name+".motif", "w")
            pwm_file.write(Bio.motifs.write([bio_motif], format='jaspar'))
            pwm_file.close()

            # call weblogo to create logos
            create_logo(bio_motif, 
                outputPath+'/html_files/'+bio_motif.name+'.motif.svg')

            # write html file
            indMotifFile = open(outputPath+"/html_files/"+allMotifs[ind][0]+".html", "w")
            indMotifFile.write("<html><head><style> td {border: 1px solid black;} .rotate{-webkit-transform:rotate(-90deg); writing-mode: tb-rl;filter: flipv fliph;white-space:nowrap;display:block} table {border-collapse:collapse;}</style><script src='http://code.jquery.com/jquery-2.1.1.min.js'></script><script src='html_files/heatMap.js'></script></head><body>\n")
            indMotifFile.write("<h1>"+allMotifs[ind][0]+"</h1>\n")
            # show logo
            indMotifFile.write("<h2>Logo</h2>\n")
            indMotifFile.write("<img src = '" + allMotifs[ind][0]+".motif.svg'>\n")
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
    if create_html:
        scoreFile = open(outputPath+"/allScores.html", "w")
        scoreFile.write("<html><head><style> td {border: 1px solid black;} .rotate{-webkit-transform:rotate(-90deg); writing-mode: tb-rl;filter: flipv fliph;white-space:nowrap;display:block} table {border-collapse:collapse;}</style><script src='http://code.jquery.com/jquery-2.1.1.min.js'></script><script src='heatMap.js'></script></head><body>\n")
        # add in some blank spaces
        scoreFile.write("<br><br><br><br>\n")
        # write score array as a matrix
        scoreFile.write("<table class='heat-map'><thead><tr><td></td>")
        for name in motifNames:
            scoreFile.write("<th class='stats-title'><span class='rotate'><a href='html_files/"+name+".html'>" + name + "</a></span></th>")
        scoreFile.write("</tr></thead>")
        scoreFile.write("<tbody>\n")
        for i in range(scoreArray.shape[0]):
            scoreFile.write("<tr class='stats-row'><th class='stats-title'><a href='html_files/"+motifNames[i]+".html'>" + motifNames[i] + "</a></th>")
            for j in range(scoreArray.shape[0]):
                scoreFile.write("<td>" + str(int(np.round(scoreArray[i][j]*100, 3))) + "</td>")
            scoreFile.write("</tr>\n")
        scoreFile.write("</tbody></table>\n")
        # add script to fix heights
        scoreFile.write("<script>$(function(){$('table th').height('10px');$('.rotate').width('10px');});</script>")
        scoreFile.write("</body></html>")

        # sort list file lines
        listFileLines = sorted(listFileLines, key = lambda x:x[0])
        for i in range(len(listFileLines)):
            count = str(i+1)
            line = listFileLines[i][1].replace('MOTIF_COUNT', count)
            listFile.write(line)
        
        listFile.write("</tbody></table>\n")
        listFile.write("</body></html>")
        listFile.close()
        scoreFile.close()
    motifListFile.close()
    mergedMetadataFile.close()

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
        default = 0.9,
        type=float)
    parser.add_argument("motifFiles",
        help="list of moti files to cluster",
        type=str,
        nargs="+")
    parser.add_argument('-familyBasedName',action='store_true', default=True)
    parser.add_argument('-createHTML',action='store_true', default=False)

    # parse arguments
    args = parser.parse_args()

    scorePath = args.scorePath
    outputPath = args.outputPath
    threshold = args.threshold
    motifFiles = args.motifFiles
    create_html = args.createHTML

    if not os.path.isdir(outputPath):
        os.mkdir(outputPath)
    if create_html:
        if not os.path.isdir(outputPath + '/html_files/'):
            os.mkdir(outputPath + '/html_files')
        else:
            for f in os.listdir(outputPath + '/html_files'):
                os.remove(outputPath + '/html_files/' + f)
    if not os.path.isdir(outputPath + '/clustered_motifs/'):
        os.mkdir(outputPath + '/clustered_motifs')
    else:
        for f in os.listdir(outputPath + '/clustered_motifs'):
            os.remove(outputPath + '/clustered_motifs/' + f)

    # read in motifs
    # find all motifs in input directory

    allMotifs = []
    motifNames = []
    for mf in sorted(motifFiles):
        (motif_name, motif, motif_id) = readMotifFile(mf)
        allMotifs.append((motif_name, motif, motif_id))
        motifNames.append(motif_name)

    # read in scores
    score_frame = pd.read_csv(scorePath, sep='\t',index_col=0)
    scoreArray = score_frame.values

    # if output directory doesn't exist, create it
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

    thresholdClusterMotifs(scoreArray, 
        threshold, 
        allMotifs, 
        motifNames, 
        outputPath, 
        create_html = create_html)
