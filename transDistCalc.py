#!/usr/bin/env python
#python 2.7.5 requires biopython
#transDist.py
#Version 1. Adam Taranto, April 2015
#Contact, Adam Taranto, adam.taranto@anu.edu.au

#################################################################################################
# ___________  ___   _   _  _____       ______ _____ _____ _____      _____   ___   _     _____ #
#|_   _| ___ \/ _ \ | \ | |/  ___|      |  _  \_   _/  ___|_   _|    /  __ \ / _ \ | |   /  __ \#
#  | | | |_/ / /_\ \|  \| |\ `--. ______| | | | | | \ `--.  | |______| /  \// /_\ \| |   | /  \/#
#  | | |    /|  _  || . ` | `--. \______| | | | | |  `--. \ | |______| |    |  _  || |   | |    #
#  | | | |\ \| | | || |\  |/\__/ /      | |/ / _| |_/\__/ / | |      | \__/\| | | || |___| \__/\#
#  \_/ \_| \_\_| |_/\_| \_/\____/       |___/  \___/\____/  \_/       \____/\_| |_/\_____/\____/#
#################################################################################################
# 
# Finds the alignment score between matched transcript pairs, and estimates bowtie2 read mapping 
# success between sets.
#
# To be used for estimating the average pairwise distance between transcript sets and appropriate 
# bowtie2 --score-min functions for mapping reads across transcriptomes for downstream Corset 
# clustering.
#################################################################################################

from __future__ import print_function
import os
import csv
import sys
import numpy as np
import pandas as pd
#import mosaic
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from os.path import basename
#import matplotlib.pyplot as plt
import scipy.stats as stats 
import pylab as pl

def main(fastaA, fastaB, gapOpen=-6, gapExtend=-3, mismatch=-5, readLength=100, scoreMinIntercept=-0.6, scoreMinSlope=-0.6, percentile=95, outFig="readPenaltyDist.pdf", outFile="alignmentStats.txt", verbose=False, recipFile=False, minLen=0, eVal=0.001, blastAvB=None, blastBvA=None, pairNames=None):

	#Read pairs in as list of tuples
	if pairNames:
		pairs = readPairs(pairNames)
	elif blastAvB and blastBvA:
		pairs = readBlast(minLen, eVal,blastAvB,blastBvA,recipFile)
	else:
		sys.exit('Must provide list of transcript pairs OR blast output files to perform reciprocal blast analysis')

	#Read in transcript seqs as dictionary keyed by name
	seqMaster = readSeqs(fastaA, fastaB)

	alignmentInfo = alignTrimStats(pairs, seqMaster, gapOpen, gapExtend, mismatch, readLength, verbose)

	#Estimate average read alignment score
	minScore = scoreMinIntercept + (scoreMinSlope * readLength)

	#Summary stats for all aligned pairs
	summaryDict = summaryStats(alignmentInfo)

	#Write output file for cluster counts from each transcriptome
	summaryFile = open(outFile,'w')

	header = "\t".join(['transNameA', 'transNameB', 'align_len', 'gaps', 'mismatches', 'align_score', 'mean_read_penalty', 'cross_map' ])
	summaryFile.write(header + "\n")
	
	totalAligns = 0
	crossMapPass = 0
	for n in alignmentInfo.iterkeys():
		totalAligns += 1
		if minScore <= alignmentInfo[n]['readScore']:
			mappable = 'SUCCESS'
			crossMapPass +=1
		if minScore > alignmentInfo[n]['readScore']:
			mappable = 'FAIL'

		alignmentScore = alignmentInfo[n]['score']
		
		summaryString = "\t".join([alignmentInfo[n]['SeqA'], alignmentInfo[n]['SeqB'], str(alignmentInfo[n]['len']), str(alignmentInfo[n]['gaps']), str(alignmentInfo[n]['mismatch']), str(alignmentScore), str(alignmentInfo[n]['readScore']), mappable])
		#print('writing summary string:')
		#print(summaryString + "\n")
		summaryFile.write(summaryString + "\n")

	summaryFile.close()

	#Print(summary stats to screen)
	print('Summary stats:')
	print("Total pairwise comparisons made: " + str(summaryDict['recordCount']))
	print("Median alignment length: " + str(summaryDict['len']))
	print("Median gaps: " + str(summaryDict['gaps']))
	print("Median mismatch: " + str(summaryDict['mismatch']))
	print("Median alignment score: " + str(summaryDict['score']))
	print("Median predicted read penalty: " + str(summaryDict['readScore']))

	goldenPenalty = writeGraph(percentile, summaryDict['readScoreList'], outFig)

	print('Proposed penalty threshold of ' + str(minScore) + ' predicted to allow crossmapping for ' + str(crossMapPass) + ' of ' + str(totalAligns) + ' pairs.' )
	print(str(percentile) + '%' ' of crossmappings are predicted to be successful with a --score-min threshold <= ' + str(goldenPenalty))

def writeGraph(pcen, readScoreList, outFig):
	a = np.array(sorted(readScoreList))
	recip = 100 - pcen
	print(a)
	print(recip)
	p = np.percentile(a, recip)
	#Fit normal distribution
	fit = stats.norm.pdf(a, np.mean(a), np.std(a))
	#Plot gaussian curve with blue dashed line
	pl.plot(a,fit,'b-')
	#Not sure about this normed setting.. normed=True
	pl.hist(a, bins=20, color='c') 
	#Vertical line at nominated percentile
	pl.axvline(p, color='b', linestyle='dashed', linewidth=2)
	#Render graphic
	pl.savefig(outFig, bbox_inches='tight')

	return p

def readPairs(pairNames):
	with open(pairNames) as src:
		data = [rec.rstrip().split('\t') for rec in src]
		cleaned_data = [map(str.strip, line) for line in data]
		pairs = [tuple(line[0:2]) for line in cleaned_data if line[0][0] != "#"]

	return pairs

def readBlast(minLen, eVal, blastTabA, blastTabB, recipFile):
	#Reads in two blast tab output files and returns a list of pairs which are reciprocal best blast hits.
	from collections import Counter
	pairs = list()

	#Find best A-B pairs
	with open(blastTabA) as src:
		data = [rec.rstrip().split('\t') for rec in src]
		cleaned_data = [map(str.strip, line) for line in data]
		lastQuery = None
		for line in cleaned_data:
			#Ignore lines begining with '#'
			if line[0][0] == "#":
				continue
			#Ignore if eVal is above threshold
			if float(line[10]) >= float(eVal):
				continue
			#Ignore in hit length is less than threshold
			if line[3] <= int(minLen):
				continue
			#If first record for this query, store the query:target pair
			if line[0] != lastQuery:	
				pairs.append(tuple(line[0:2]))
			#Update previous query ID
			lastQuery = line[0]

	#Append best B-A pairs as A-B
	with open(blastTabB) as src:
		data = [rec.rstrip().split('\t') for rec in src]
		cleaned_data = [map(str.strip, line) for line in data]
		lastQuery = None
		for line in cleaned_data:
			if line[0][0] == "#":
				continue
			#Ignore if eVal is above threshold
			if float(line[10]) >= float(eVal):
				continue
			#Ignore in hit length is less than threshold
			if line[3] <= int(minLen):
				continue
			if line[0] != lastQuery:	
				flipPair = (line[1],line[0])
				pairs.append(flipPair)
			lastQuery = line[0]

	#Retrieve list items that occur more than once
	recipPairs = [k for (k,v) in Counter(pairs).iteritems() if v > 1]

	if recipFile:
		baseA = os.path.splitext(basename(blastTabA))[0]
		baseB = os.path.splitext(basename(blastTabB))[0]
		outName = baseA + "_" + baseB + "_reciprocal_pairs.tab"
		outFile = open(outName,'w')
		header = "\t".join(['#SetA','SetB'])
		outFile.write(header + "\n")
		
		for a,b in recipPairs:
			abPair = "\t".join([a,b])
			outFile.write(abPair+ "\n")

		outFile.close()
	return recipPairs

def readSeqs(fastaA, fastaB):
	#Create empty dictionary
	seqMaster={}

	#Populate dictionary with master set of fasta records
	for seq_record in SeqIO.parse(fastaA, "fasta"):
		seqMaster[seq_record.id]=str(seq_record.seq)

	for seq_record in SeqIO.parse(fastaB, "fasta"):
		seqMaster[seq_record.id]=str(seq_record.seq)

	return seqMaster

def median(valList):
    return np.median(np.array(valList))

def format_alignment(align1, align2, scoreBitscore): 
	#Remake the Pairwise2 alignment formatter such that it works.
	score = scoreBitscore[0][0]
	bitscore = scoreBitscore[0][1]
	#Initialize new list
	s = [] 
	#Add aligned seq A
	s.append("%s\n" % align1) 

	#Add linker line
	for a,b in zip(align1,align2):
		if a != b and a != "-" and b != "-":
			s.append("%s" % ".")
		elif a == b and a != "-" and b != "-":
			s.append("%s" % "|")
		elif a == "-" or b == "-":
			s.append("%s" % " ")

	#Add aligned seq B
	s.append("\n%s\n" % align2)

	#Report scores
	s.append("Score= %g\n" % score)
	s.append("BitScore= %g\n" % bitscore)
	
	return ''.join(s)

def getBitScore(alignedA, alignedB):
	score = 0
	scoremat = pd.read_csv(os.path.join(os.path.dirname(__file__),'EDNAFULL.txt'),header=0, index_col=0, sep="\s+")
	for xx, yy in zip(alignedA, alignedB):
		xx = '*' if xx == '-' else xx
		yy = '*' if yy == '-' else yy
		score += scoremat.ix[xx, yy]
	
	bitscore = 1.0 * score / len(alignedA)
	scoreBitscore = [(score, bitscore)]

	return scoreBitscore

def alignTrimStats(pairs, seqMaster, gapOpen, gapExtend, mismatch, readLength, verbose):
	#Initialize dict and counter
	n = 0
	alignDict = {}

	#Check for len diff > 50%
	for a,b in pairs:
		#Retrieve sequences by name
		aSeq = seqMaster[a]
		bSeq = seqMaster[b]

		#Test for min comparable length
		if len(aSeq) < 0.5 * len(bSeq):
			continue
		if len(bSeq) < 0.5 * len(aSeq):
			continue

		alignDict[n] = {}

		alignDict[n]['SeqA'] = list()
		alignDict[n]['SeqB'] = list()
		alignDict[n]['aligned'] = list()
		alignDict[n]['len'] = int()
		alignDict[n]['gaps'] = int()
		alignDict[n]['mismatch'] = int()
		#alignDict[n]['score'] = float()
		alignDict[n]['readScore'] = float()

		#Set sequence names
		alignDict[n]['SeqA'] = a
		alignDict[n]['SeqB'] = b

		#Generate primary alignment
		#Returns tuple of Sequence 1 aligned, sequence 2 aligned
		#firstAlign = [mosaic.pairwisealign(aSeq, bSeq, gapopen=gapOpen, gapextend=gapExtend, AA='false')]

		#Pairwise2. Global align. Custom penalties. Same penalties for both seqs.
		#Identical characters are given 2 points, 1 point is deducted for each non-identical character.
		#5 points are deducted when opening a gap, and 3 points are deducted when extending it.
		allGlobal = pairwise2.align.globalms(aSeq, bSeq, 2, -1, -5, -3)

		#Package resulting global alignment as a tuple
		firstAlign = [(allGlobal[0][0],allGlobal[0][1])]

		#Trims gap positions at start and end of alignment, return alignment tuple
		trimmedAlign = trimAlign(firstAlign)

		#Returns tuple with score and bitscore
		scoreBitscore = getBitScore(trimmedAlign[0][0],trimmedAlign[0][1])

		if verbose:
			print('\nAlignment: ' + a + '-' + b)
			print(format_alignment(trimmedAlign[0][0],trimmedAlign[0][1], scoreBitscore))

		#Final alignment
		alignDict[n]['aligned'] = trimmedAlign

		# Calculate the number of mismatches in alignment
		alignDict[n]['mismatch'] = countMismatch(trimmedAlign[0][0],trimmedAlign[0][1])

		# Calculate the number of gaps
		alignDict[n]['gaps'] = countGaps(trimmedAlign[0][0],trimmedAlign[0][1])

		#Store trimmed alignment length
		alignDict[n]['len'] = len(trimmedAlign[0][0])

		#Calculate average penalty score for read length k
		alignDict[n]['readScore'] = readScore(alignDict[n]['gaps'], alignDict[n]['mismatch'], len(trimmedAlign[0][0]), gapExtend, gapOpen, mismatch, readLength)

		#Whole alignment score
		alignDict[n]['score'] = scoreBitscore[0][0]
		
		#Increment key counter
		n += 1

	return alignDict

def trimAlign(firstAlign):
	#Scrub trailing and leading gaps
	#Scrub double gap positions
	nonGapPosition = 0
	newA = list()
	newB = list()
	#For each position pair
	for seqA,seqB in firstAlign:
		for (a,b) in zip(seqA, seqB):
			#If one character is not a gap
			if a != "-" or b != "-":
				#Log if this is the first zero gap position
				if a != "-" and b != "-":
					nonGapPosition =+ 1
				#If this is any position after and including the first zero-gap
				#write to new position list
				if nonGapPosition > 0:
				 	newA.append(a)
				 	newB.append(b)

	#Convert list of positions back to string
	trimmedA = ''.join(newA)
	trimmedB = ''.join(newB)

	#Flip and repeat from other end
	revAlign = [(trimmedA[::-1],trimmedB[::-1])]

	nonGapPosition = 0
	revlistA = list()
	revlistB = list()

	#For each position pair
	for revA,revB in revAlign:
		for (a,b) in zip(revA, revB):
			#If one character is not a gap
			if a != "-" or b != "-":
				#Log if this is the first zero gap position
				if a != "-" and b != "-":
					nonGapPosition =+ 1
				#If this is any position after and including the first zero-gap
				#write to new position list
				if nonGapPosition > 0:
				 	revlistA.append(a)
				 	revlistB.append(b)

	finalA = ''.join(revlistA)
	finalB = ''.join(revlistB)

	#Make tuple
	finalAlign = [(finalA[::-1],finalB[::-1])]

	return finalAlign

def readScore(gapCount, mismatchCount, alignLen, gapExtend, gapOpen, mismatch, readLength):
	meanReadGap = (float(gapCount)/float(alignLen)) * readLength
	meanReadMismatch = (float(mismatchCount)/float(alignLen)) * readLength
	scoreMismatch = mismatch * meanReadMismatch
	scoreGap = (gapOpen + gapExtend) * meanReadGap
	score = scoreMismatch + scoreGap
	
	#print('mismatches found: ' + str(mismatchCount))
	#print('alignment length: ' + str(alignLen))
	#print('mismatches/len: ' + str(float(mismatchCount)/float(alignLen)))
	#print('Gaps found: ' + str(gapCount))
	#print('Testing for read of length: ' + str(readLength))
	#print('\nCalculating mean read penalty:')
	#print('mean gaps per read: ' + str(meanReadGap))
	#print('mean mismatches per read: ' + str(meanReadMismatch))
	#print('Mismatch penalty * mean mismatch per read: ' + str(scoreMismatch))
	#print('Gap penalty * mean gaps per read: ' + str(scoreGap))
	#print('Gap penalty + mismatch penalty per read: ' + str(score))

	return score
	
def countGaps(seqA, seqB):
	gaps = sum(a == "-" or b == "-" for (a,b) in zip(seqA, seqB))
	return gaps

def countMismatch(seqA, seqB):
	mismatches = sum(a != b and a != "-" and b != "-" for (a,b) in zip(seqA, seqB))
	return mismatches

def summaryStats(alignDict):
	summary = {}
	summary['len'] = list()
	summary['gaps'] = list()
	summary['mismatch'] = list()
	summary['score'] = list()
	summary['readScore'] = list()

	count = 0

	for n in alignDict:
		summary['len'].append(alignDict[n]['len'])
		summary['gaps'].append(alignDict[n]['gaps'])
		summary['mismatch'].append(alignDict[n]['mismatch'])
		summary['score'].append(alignDict[n]['score'])
		summary['readScore'].append(alignDict[n]['readScore'])
		count += 1

	summaryStats = {}
	summaryStats['len'] = float()
	summaryStats['gaps'] = float()
	summaryStats['mismatch'] = float()
	summaryStats['score'] = float()
	summaryStats['readScore'] = float()
	summaryStats['recordCount'] = int()
	summaryStats['readScoreList'] = list()

	summaryStats['len'] = median(summary['len'])
	summaryStats['gaps'] = median(summary['gaps'])
	summaryStats['mismatch'] = median(summary['mismatch'])
	summaryStats['score'] = median(summary['score'])
	summaryStats['readScore'] = median(summary['readScore'])
	summaryStats['recordCount'] = count
	summaryStats['readScoreList'] = summary['readScore']

	#Should probably summarise total pass / fail
	return summaryStats

if __name__== '__main__':
	###Argument handling / alphabet soup.
	arg_parser = argparse.ArgumentParser(description='Finds the alignment score between matched transcript pairs, and estimates bowtie2 read mapping success between sets.')
	arg_parser.add_argument("-r", "--readLength", default=100, type=int, help="Length of reads intended for mapping to transcriptomes (after quality trimming).")
	arg_parser.add_argument("-i", "--scoreMinIntercept", default=-0.6, type=float, help="Intercept value for Bowtie2 --score-min function.")
	arg_parser.add_argument("-s", "--scoreMinSlope", default=-0.6, type=float, help="Slope value for Bowtie2 --score-min function.")
	arg_parser.add_argument("-n", "--pairNames", help="Two column tab delimited file containing names for sequences from fastaA and fastaB to be aligned.")
	arg_parser.add_argument("-x", "--blastAvB", help="Blast tab result file for fastaA query against fastaB subject")
	arg_parser.add_argument("-y", "--blastBvA", help="Blast tab result file for fastaB query against fastaA subject")
	arg_parser.add_argument("-a", "--fastaA", required=True, help="Multifasta set A.")
	arg_parser.add_argument("-b", "--fastaB", required=True, help="Multifasta set B.")
	arg_parser.add_argument("-g", "--gapOpen", default=-5, type=int, help="Penalty for opening a gap.")
	arg_parser.add_argument("-e", "--gapExtend", default=-3, type=int, help="Penalty for extending a gap.")
	arg_parser.add_argument("-m", "--mismatch", default=-6, type=float, help="Pentaly for a mismatched position.")
	arg_parser.add_argument("-o", "--outFile", default="alignmentStats.txt", help="Write output table to file.")
	arg_parser.add_argument("-w", "--recipFile", action='store_true', default=False, help="Write reciprocal blast pairs to file. True if flag is set.")
	arg_parser.add_argument("-l","--minLen", default=0, help="Minimum length of hit to consider valid.")
	arg_parser.add_argument("-E","--eVal", default=0.001, type=float, help="Minimum eval to consider valid pair.")
	arg_parser.add_argument("-v", "--verbose", action='store_true', default=False, help="Write formatted alignments to screen.")
	arg_parser.add_argument("-p", "--percentile", default=95, type=int, help="Report the estimated read penalty at the lower bound of this percentile.")
	arg_parser.add_argument("-f", "--outFig", default="readPenaltyDist.pdf", help="Write read penalty distribution histogram to this file.")

	#Parse arguments
	args = arg_parser.parse_args()

	###Variable Definitions
	readLength = args.readLength
	scoreMinIntercept = args.scoreMinIntercept
	scoreMinSlope = args.scoreMinSlope
	pairNames = args.pairNames
	blastAvB = args.blastAvB
	blastBvA = args.blastBvA
	fastaA = args.fastaA
	fastaB = args.fastaB
	gapOpen = args.gapOpen
	gapExtend = args.gapExtend
	mismatch = args.mismatch
	outFile = args.outFile
	recipFile = args.recipFile
	minLen = args.minLen
	eVal = args.eVal
	verbose = args.verbose
	percentile = args.percentile
	outFig = args.outFig

	main(fastaA, fastaB, gapOpen, gapExtend, mismatch, readLength, scoreMinIntercept, scoreMinSlope, percentile, outFig, outFile, verbose, recipFile, minLen, eVal, blastAvB, blastBvA, pairNames)
