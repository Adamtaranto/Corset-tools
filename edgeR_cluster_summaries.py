#!/usr/bin/env python
#python 2.7.5 requires biopython
#edgeR_cluster_summaries.py
#Version 1. Adam Taranto, January 2016
#Contact, Adam Taranto, adam.taranto@anu.edu.au

'''
Report misc annotation records for members of significant clusters from EdgeR DE topTags reports.

Example usage:
./edgeR_cluster_summaries.py \
--infile edgeR_topTags_SetA_1vsSetB_1.csv \
--outDir temp \
--clust clusters.txt \
--seqFiles SetA_transcripts.fa SetB_transcripts.fa \
--labels SetA SetB \
--grepme AB_annotdescriptions.txt AB_annotgaff.txt AB_kegg.txt AB_topblast.txt \
--topx 50 \
--maxfdr 0.05 \
--minfc 1 \
--writeClusters
'''

import os
import re
import csv
import sys
import copy
import gzip
import logging
import argparse
import pandas as pd
from os.path import basename

def iterFasta(path):
    """Iterates over the sequences of a fasta file"""
    logging.info("Loading fasta files from %s" % path)
    name = None
    seq = []
    if path.endswith('.gz') or path.endswith('.gz"'):
        handle = gzip.open(path)
    else:
        handle = open(path)
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name: #then first seq has been read and it is time to yield
                yield (name, ''.join(seq))
            # Update the name to current seq and blank the sequence
            name = line.strip('>').split()[0]
            seq = []
        else:
            # Extend the current sequence with line
            seq.append(line)
    if name:
        # Yield the final sequence at end of file.
        yield (name, ''.join(seq))
    handle.close()

def getFasta(fastaPath):
    """Write fasta to dictionary, key by scaffold name."""
    seqDict = dict()
    for name, seq in iterFasta(fastaPath):
        seqDict[name] = seq
    return seqDict

def getMultiFasta(args):
    multifasta = dict()
    for x in args.seqFiles:
        name = args.labels[args.seqFiles.index(x)]
        multifasta[name] = getFasta(x)
    return multifasta

def outPathCheck(args):
    """Check if the temporary folder exists, else creates it"""
    logging.info('Checking for output directory.')
    outFolder = os.path.abspath(args.outDir)
    if not os.path.isdir(outFolder):
        logging.info('Creating output directory.')
        os.makedirs(outFolder)
    else:
        logging.info('Output directory exists.')

def getclusters(mapFile):
    '''Read clusters to dictionary'''
    #Read transcript-to-cluster mapping file into dict object.
    #Sample data row:
    #TranscriptID   ClusterID
    #nnt3Ldvymb Cluster-0.0
    logging.info('Mapping transcripts to parent clusters.')
    mapFile = open(mapFile, 'rt')
    readMap = csv.reader(mapFile,delimiter='\t')
    clustMem={}
    #Write records for seqs in name file to new fasta
    for row in readMap:
        transID=row[0]
        clustID=row[1]
        if clustID not in clustMem:
            clustMem[clustID] = {}
            clustMem[clustID]['Seqs'] = list()
        clustMem[clustID]['Seqs'].append(transID)
    mapFile.close()
    return clustMem

def makeDEtable(args):
    '''Convert EdgeR results to pandas table and filter'''
    df = pd.read_csv(args.infile)
    df.columns = ['ClusterID','logFC','logCPM','LR','PValue','FDR']
    df = df[(df.logFC.abs() >= args.minfc) & (df.FDR <= args.maxfdr)]
    df = df.sort_values(by='FDR', ascending=True).reset_index().drop('index', axis=1)
    if args.topx:
        df = df[:args.topx]
    return df

def getSeqLabels(seqDict):
    '''Map transcript names to the set label '''
    logging.info('Mapping transcript names to parent transcriptome.')
    setLabels = dict()
    for x in seqDict.keys():
        allLabels = setLabels.copy()
        setLabels = dict.fromkeys(seqDict[x].keys(),x)
        allLabels.update(setLabels)
        setLabels = allLabels.copy()
    return allLabels

def getClustID(resultTable):
    clustList = list()
    for x in resultTable['ClusterID']:
        clustList.append(str(x))
    return clustList

def getMemSum(clusterMap,seqLabels,topTags,labels):
    '''Generate transcriptome membership summary by clustID.'''
    logging.info('Calculating transcriptome membership counts for transcript clusters.')
    setCounts = dict()
    blankSets = dict()
    for label in labels:
        blankSets[label] = 0
    clustList = getClustID(topTags)
    for clust in clustList:
        setCounts[clust] = copy.deepcopy(blankSets)
        for x in clusterMap[clust]['Seqs']: #lookup transcripts belonging to clusters
            fromSet = seqLabels[x]
            setCounts[clust][fromSet] += 1
    return setCounts

def addMemSum(resultTable,memSum,labels):
    ''' Add membership count columns to result table '''
    logging.info('Append membership counts to DE dataframe.')
    counts = pd.DataFrame.from_dict(memSum,orient='index',dtype=int).reset_index()
    counts['Totals'] = counts[labels].sum(axis=1)
    merged_results = pd.merge(resultTable, counts, how='left',left_on='ClusterID', right_on='index').drop('index', axis=1)
    return merged_results

def writeAnnotLines(clusterMap, annotFile, outFile, queryClustList):
    handle = open(outFile, "w")
    for queryClust in queryClustList:
        handle.write('#' + str(queryClust) + '\n')
        with open(annotFile) as f:
            strings = clusterMap[queryClust]['Seqs']
            for line in f:
                if any(s in line for s in strings):
                    handle.write(line)
    handle.close()
    pass

def formatFasta(string, chunk_size=60):
    num_chunks = len(string)/chunk_size
    if (len(string) % chunk_size != 0):
        num_chunks += 1
    output = []
    for i in range(0, num_chunks):
        output.append(string[chunk_size*i:chunk_size*(i+1)])
    return output
    
def writeClustFasta(args,allTranscripts,clusterMap,queryClustList):
    logging.info('Writing member transcripts to fasta by cluster.')
    for clust in queryClustList:
        outFile = os.path.join(args.outDir, clust + '.fa')
        if not os.path.isfile(outFile):
            handle = open(outFile, "w")
            seqNames = clusterMap[clust]['Seqs']
            for x in allTranscripts:
                for y in seqNames:
                    try: 
                        if allTranscripts[x][y]:
                            header = '>' + '\t'.join((y,x,clust)) + '\n'
                            handle.write(header)
                            seqLines = formatFasta(allTranscripts[x][y], chunk_size=60)
                            for line in seqLines:
                                handle.write(line + "\n")
                    except:
                        pass
        else:
            logging.info('Previous cluster fasta found for %s - skipp to next.' % os.path.isfile(outFile))
            continue
        handle.close()
    pass

def mainArgs():
    """Process command-line arguments"""
    
    parser = argparse.ArgumentParser(
        description='Filter DE cluster report from EdgeR and report lines mentioning \
        members of significant clusters from misc annotation files provided .'
        )
    parser.add_argument("--clust",
                        type=str,
                        help="Corset cluster file 'transcript_name Cluster_ID'"
                        )
    parser.add_argument("-x","--topx", 
                        default=None, 
                        type=int,
                        help="If set returns top 'x' DE clusters that pass cut-offs, sorted by FDR."
                        )
    parser.add_argument("--maxfdr",
                        type=float,
                        default=0.05,
                        help="Only records below this FDR threshold will be considered."
                        )
    parser.add_argument("--minfc",
                        type=float,
                        default=1, 
                        help="Only records above this log2 FC threshold will be considered. \
                        Note: Default of 1 corresponds to 2-fold change as log2(2) = 1")
    parser.add_argument("-i","--infile",
                        type=str,
                        required=True,
                        help="EdgeR topTags result table. Format: 'ClusterID,logFC,logCPM,LR,PValue,FDR'"
                        )
    parser.add_argument("-o","--outDir",
                        type=str, 
                        default="cluster_reports", 
                        help="Output directory."
                        )
    parser.add_argument("--seqFiles",
                        type=str,
                        nargs='+',
                        required=True,
                        default=None,
                        help="List of fasta files containing transcript sequences. Must be uniquely named."
                        )
    parser.add_argument("--labels",
                        nargs='+',
                        type=str,
                        required=True,
                        help="List of labels for provided transcript fastas."
                        )
    parser.add_argument("--grepme",
                        nargs='+',
                        default=None,
                        required=True,
                        help="Space delimited list of misc annotation files. Will attempt to collect and print any line \
                        mentioning a given cluster that passes threshold."
                        )
    parser.add_argument("--writeClusters",
                        default=False,
                        action='store_true',
                        help="If set will write cluster members to fasta for each cluster."
                        )

    args = parser.parse_args()
    return args

def main():
    logging.basicConfig(level=logging.INFO, format=("%(asctime)s - %(funcName)s - %(message)s"))
    args = mainArgs()
    outPathCheck(args)

    if len(args.labels) != len(args.seqFiles):
            logging.info('Label count different to fasta count. Exiting.')
            sys.exit(0)

    ''' Import cluster file '''
    #clusterMap[clusterID] = ['list of transcripts']
    clusterMap = getclusters(args.clust)

    ''' Import fasta files '''
    # {Set_label : {seq_name : seq} }
    allTranscripts = getMultiFasta(args)

    '''Map transcript names to transcriptome set label '''
    #ta_037465786 : Grandin
    seqLabels = getSeqLabels(allTranscripts)
    
    ''' Import and filter results csv '''
    topTags = makeDEtable(args)

    '''Get Transcriptome Membership summary for each cluster '''
    # { clusterID :{set_labelA : 1, set_labelB : 2}, ...}
    memSummary = getMemSum(clusterMap,seqLabels,topTags,args.labels)

    ''' Append membership count columns to result table '''
    # columns now = ['ClusterID','logFC','logCPM','LR','PValue','FDR'] + [args.labels] + 'Totals'
    tagsWithCounts = addMemSum(topTags,memSummary,args.labels)

    '''Write cluster membership report'''
    # Create inputfile base name
    infileBase = os.path.basename(args.infile)
    infileBase = os.path.splitext(infileBase)[0]
    # Optional write significant cluster report with membership counts
    if args.writeClusters:
        reportOut = os.path.join(args.outDir, infileBase + "_significant_cluster_report" + '.tab')
        if os.path.isfile(reportOut):
            logging.info('Overwriting existing report file: %s' % str(reportOut))
        tagsWithCounts.to_csv(reportOut, sep='\t', encoding='utf-8')

    ''' Report transcript annotations '''
    for annotFile in args.grepme:
        pathbase = os.path.basename(annotFile)
        pathbase = os.path.splitext(pathbase)[0]
        reportOut = os.path.join(args.outDir, pathbase + '_' + infileBase + "_significant_cluster_annotation_records" + '.txt')
        writeAnnotLines(clusterMap,annotFile,reportOut,queryClustList=tagsWithCounts['ClusterID'])

    '''Write transcripts to cluster fasta files'''
    writeClustFasta(args,allTranscripts,clusterMap,queryClustList=tagsWithCounts['ClusterID'])

if __name__== '__main__':
    main()