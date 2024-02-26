#!/usr/bin/env python

import cairo
import bioinfo
import re
import argparse
#def get_args():
	#parser = argparse.ArgumentParser(description="")
	#parser.add_argument("-f", "--filename", help="Path to input FASTA file", required=True)
	#parser.add_argument("-m", "--motiffile", help="Path to motifs text file", type=str)
	#return parser.parse_args()

#args=get_args()	

#f = args.filename  # input file name
#o = args.motiffile # output file name

# Assess the fasta file sequence: is gDNA? or transcript? if gDNA is motifs containing T, if not is motifs containing U

# Y is for C and U nucleotide

# bioinfo.oneline_fasta(f): tweak this functon so that it does not output a different file



#class Motif:
#   def __repr__(self) -> str:

class DNAsequence(object):
	def __init__(self,fa_sequence:str):
		self.sequence = fa_sequence
	def get_intron(self):
		introns = re.split('[A-Z]+',self.sequence)
		[introns.remove(x) for x in introns if x == '']
		return introns
	def get_exon(self):
		exons = re.split('[a-z]+',self.sequence)
		[exons.remove(x) for x in exons if x == '']
		return exons
	
test = DNAsequence('aaaaacgctaacgcctgctAAAAAAAGCATGAAAACAtagaaaaaaaa')



