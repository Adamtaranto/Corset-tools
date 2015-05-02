#!/usr/bin/env python
#python 2.7.5 requires biopython
#crossClustCount.py
#Version 1. Adam Taranto, April 2015
#Contact, Adam Taranto, adam.taranto@anu.edu.au

#Take two transcriptomes and a Corset cluster map.
#Determine number of member transcripts in each cluster that belong to each of the input transcriptomes.

import os
import csv
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import basename

def main(transFastaX=None, transFastaY=None, transNameX=None, transNameY=None, clustMap=None, outFile='CountClusterMembers.txt'):

	#If main is imported to python term as function, check that required inputs are provided
	if transFastaX is None:
		sys.exit('Missing transcriptome fasta X')

	if transFastaY is None:
		sys.exit('Missing transcriptome fasta Y')

	if clustMap is None:
		sys.exit('No cluster mapping file provided')

	#Set transcriptome set names same as input fastas if none supplied 
	if transNameX is None:
		transNameX = basename(os.path.splitext(transFastaX)[0])
		
	if transNameY is None:
		transNameY = basename(os.path.splitext(transFastaY)[0])

	#Returns dict keyed by transcript name, mapping transcriptome set names
	transSets = tranDict(transFastaX, transFastaY, transNameX, transNameY)

	#Returns dict keyed by cluster name, with first level dict 'Seq' containing a list of transcripts belonging to cluster.
	#clustMem[clustID]['Seqs']
	clustMem = mapDict(clustMap)

	#Appends sub-dictionaries to clust mem to store member counts for each transcriptome input file
	clustMemCounts = memStats(clustMem, transSets, transNameX, transNameY)

	#Returns dict with summary of clusters containing zero members from one or more of the reference transcriptomes
	summaryDict = summaryStats(clustMemCounts)

	#Write output file for cluster counts from each transcriptome
	summaryFile = open(outFile,'w')

	header = "\t".join(['clusterID', transNameX, transNameY, 'Total_Members'])
	summaryFile.write(header + "\n")
	
	for cluster in clustMemCounts.iterkeys():
		#print 'this is the cluster: '
		#print clustMemCounts[cluster]
		#print 'x count =' + str(clustMemCounts[cluster]['countX'])
		#print 'y count =' + str(clustMemCounts[cluster]['countY'])

		summaryString = "\t".join([cluster, str(clustMemCounts[cluster]['countX']), str(clustMemCounts[cluster]['countY']), str(clustMemCounts[cluster]['countT'])])
		#print 'writing summary string:' 
		#print summaryString + "\n"
		summaryFile.write(summaryString + "\n")

	summaryFile.close()

	#Print summary stats to screen
	print "\n" + 'Summary stats:'
	print 'Clusters with 0 members from ' + transNameX + ': ' + str(summaryDict['zeroSetX'])
	print 'Clusters with 0 members from ' + transNameY + ': ' + str(summaryDict['zeroSetY'])
	print 'Clusters with 0 members from either ' + transNameX + ' or ' + transNameY + ': ' + str(summaryDict['zeroBoth'])


def tranDict(transFastaX, transFastaY, transNameX, transNameY):
	#Create dictionary keyed by transcript name, mapping transcripts to their host transcriptome
	#Create empty dictionary
	transSets={}

	#Read transcript names in as keys, add set name as value.
	for seq_record in SeqIO.parse(transFastaX, "fasta"):
		#Check if transcript name/key is unique
		if seq_record.id not in transSets:
			transSets[seq_record.id]=str(transNameX)
		#If transcript name already exists exist with error
		else:
			sys.exit('Duplicate transcript name: ' + seq_record.id + ' in ' + transNameX)

	for seq_record in SeqIO.parse(transFastaY, "fasta"):
		if seq_record.id not in transSets:
			transSets[seq_record.id]=str(transNameY)
		else:
			sys.exit('Duplicate transcript name: ' + seq_record.id + ' in ' + transNameY)

	#print "This in the transcript to set dict:"
	#print transSets
	return transSets

def mapDict(clustMap):
	#Read transcript-to-cluster mapping file into dict object.
	#Sample data row:
	#TranscriptID	ClusterID
	#nnt3Ldvymb	Cluster-0.0

	mapFile = open(clustMap, 'rt')
	readMap = csv.reader(mapFile,delimiter='\t')

	clustMem={}

	#Write records for seqs in name file to new fasta
	for row in readMap:
		transID=row[0]
		clustID=row[1]
		#print 'reading transid'
		#print row[0]
		#print 'reading clustid'
		#print row[1]
		
		if clustID not in clustMem:
			clustMem[clustID] = {}
			clustMem[clustID]['Seqs'] = list()

		clustMem[clustID]['Seqs'].append(transID)

	mapFile.close()

	#print 'clustmem dict, no counts'
	#print clustMem
	return clustMem


def memStats(clustMem, transSets, transNameX, transNameY):
	
	for clustID in clustMem:
		#For each cluster check each transcripts parent set.
		clustMem[clustID]['countT'] = int()
		clustMem[clustID]['countX'] = int()
		clustMem[clustID]['countY'] = int()

		for transID in clustMem[clustID]['Seqs']:
			
			clustMem[clustID]['countT'] += 1
			
			if transID not in transSets:
				print 'No set name found for transcript: ' + transID

			else:
				if transSets[transID] == transNameX:
					clustMem[clustID]['countX'] += 1
				#
				if transSets[transID] == transNameY:
					clustMem[clustID]['countY'] += 1

	clustMemCounts = clustMem

	return clustMemCounts

def summaryStats(clustMemCounts):
	summary = {}
	summary['zeroSetX'] = int(0)
	summary['zeroSetY'] = int(0)
	summary['zeroBoth'] = int(0)

	for cluster in clustMemCounts:
		if clustMemCounts[cluster]['countX'] == 0:
			summary['zeroSetX'] += 1
		if clustMemCounts[cluster]['countY'] == 0:
			summary['zeroSetY'] += 1
		if clustMemCounts[cluster]['countX'] == 0 and clustMemCounts[cluster]['countY'] == 0:
			summary['zeroBoth'] += 1

	return summary

if __name__== '__main__':
	###Argument handling.
	arg_parser = argparse.ArgumentParser(description='For Corset clusters generated from two transcriptomes this program will report the number of transcripts from each transciptome that are members of each cluster.')
	arg_parser.add_argument("-X","--transFastaX", required=True, help="Fasta file for transcriptome X")
	arg_parser.add_argument("-Y","--transFastaY", required=True, help="Fasta file for transcriptome Y")
	arg_parser.add_argument("-x","--transNameX", help="Unique lable for transcriptome X")
	arg_parser.add_argument("-y","--transNameY", help="Unique lable for transcriptome Y")
	arg_parser.add_argument("-o","--outFile", default='CountClusterMembers.txt', help="Location for summary file to be written to.")
	arg_parser.add_argument("-c","--clustMap", required=True, help="Corset transcript-to-cluster mapping file.")
	args = arg_parser.parse_args()

	###Variable Definitions
	transFastaX=args.transFastaX
	transFastaY=args.transFastaY
	transNameX=args.transNameX
	transNameY=args.transNameY
	outFile=args.outFile
	clustMap=args.clustMap

	main(transFastaX, transFastaY, transNameX, transNameY, clustMap, outFile)