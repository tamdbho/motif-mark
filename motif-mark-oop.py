#!/usr/bin/env python

# Author: Tam Ho
# Collaborator: Sydney Hamilton, Anna Grace Welch

import cairo
import re
import argparse
def get_args():
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-f", "--filename", help="Path to input FASTA file", required=True)
	parser.add_argument("-m", "--motiffile", help="Path to motifs text file", type=str)
	return parser.parse_args()

args=get_args()	

f = args.filename  # input file name
m = args.motiffile # motif file name

def oneline_fasta(f):
    '''Takes a FASTA file with wrapped sequence line and output same FASTA file but with sequences as one line'''
    with open (f,"r") as f_in, open (f'{f}_oneline',"w") as f_out:
        aaseq = ""
        line = f_in.readline()
        f_out.write(line)
        for line in f_in:
            if line[0] == ">":
                f_out.write(aaseq+"\n")
                f_out.write(line)
                aaseq = ""
            else: 
                aaseq += line.strip()
        f_out.write(aaseq+"\n")

oneline_fasta(f)

class DNAsequence: # class DNAsequence is initiated as we iterate through FASTA file. Argument: sequence line
    def __init__(self, sequence:str):
        '''This function takes in the DNA sequence in FASTA file and initiate class DNA sequence with the sequence string and sequence length'''
        self.sequence = sequence
        self.length = len(sequence)
    def get_exon_start(self):
        '''This function uses self.sequence string to find the start position of exon (defined by capitalized nucleotide)'''
        exon_start = re.search('[A-Z]+',self.sequence)
        exon_start_pos = exon_start.start()
        return exon_start_pos
    def get_exon_end(self):
        '''This function uses self.sequence string to find the end position of exon (defined by capitalized nucleotide)'''
        exon_end = re.search('[A-Z]+',self.sequence)
        exon_end_pos = exon_end.end()
        return exon_end_pos
    def get_exon_length(self):
        '''This function uses self.sequence string to find the end position of exon (defined by capitalized nucleotide)'''
        exon = re.search('[A-Z]+',self.sequence)
        exon_strt_pos = exon.start()
        exon_end_pos = exon.end()
        exon_length = exon_end_pos - exon_strt_pos
        return exon_length

class Motif:    # class Motif is mainly created to store the translate_motif() function. The idea is so that we can call on this class and function in class Transcript. 
    def __init__(self,motif:str):
        '''This function takes in a single motif string as argument and initiate class Motif with the motif sequence'''
        self.motif = motif
    def translate_motif (self):
        '''This function translate the motif string into its regular expression format. Capitalize the motif string for simplicity since motif would bind to DNA sequence in both intron and exon regions.
        (Later on we will also have to capitalize all nuceotides in the DNA sequence so that we can find these motifs in sequence regardless of intron/exon region)'''
        NA_notation = { 
        'A':'[A]',
        'C':'[C]',
        'G':'[G]',
        'U':'[T]',
        'T':'[TU]',
        'W':'[AT]',
        'S':'[CG]',
        'M':'[AC]',
        'K':'[GT]',
        'R':'[AG]',
        'Y':'[CT]',
        'B':'[CGT]',
        'D':'[AGT]',
        'H':'[ACT]',
        'V':'[ACG]',
        'N':'[ACGT]'} 
        motif = self.motif.upper()  # another way to do this is {'A':[Aa]} so we can just leave both motif and sequence as is and not capitalizing them. 
        motif_regex = ""
        for nt in motif:
            motif_regex += NA_notation[nt]
        return motif_regex 

class Transcript: # class Transcript takes in DNA sequence in FASTA file and motif sequence in motif file and will be used to find motif in DNA sequence
    def __init__(self, sequence:str, all_motifs:list):
        '''This function takes in DNA sequence string and a list of all motifs (need to loop through motif file first to get this list)'''
        self.dnaseq = sequence
        self.motifs = all_motifs
    def find_motif(self): #list of motif
        '''This function takes in list of all motifs and output the position of each motif in the sequence line. Return a dictionary of motif as key and all the positions (start,end) in a list as value'''
        matches = {}
        for motif in self.motifs:
            motif_class = Motif(motif) # for every motif in the list of all motifs, initiate class Motif so we can "translate" the motif into each regex format
            motif_pattern = motif_class.translate_motif() 
            matched = [m.span() for m in re.finditer(motif_pattern, self.dnaseq.upper())] # returns list of (start,end) positions of the current motifs in DNA sequence
            if not matches.get(motif):
                matches[motif] = matched # return dictionary of motif name as key and list of (start,end) positions as value
            else:
                matches[motif].append(matched)
        return matches

# Create image surface - scaling for sequence length
with open(f'{f}_oneline',"r") as input:
    lengths = []
    seq_width = 100
    for line in input:
        line = line.strip()
        if line[0] == ">" :
            header = line.split()
            gene_name = header[0][1:]
        else:
            sequence = DNAsequence(line)
            lengths.append(len(line))
            seq_width+=100

width, height = max(lengths) + 50, seq_width 
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width , height)
context = cairo.Context(surface)
context.set_source_rgb(1, 1, 1)
context.paint()

class Drawing: # This class is used to store all pycairo function needed to draw sequence, exons, and motifs
    def __init__ (self,x,y):
        '''Initiate class Drawing with x-coordiate of all sequences (should be the same for all) and y-coordinate (should be different for each sequence)'''
        self.x = x
        self.y = y
    def draw_gene(self, name, length):
        '''Draw the the entire gene with one line (doesn't matter exon or intron): input x- and y- coordinate and length of each sequence'''
        context.set_line_width(2)
        context.set_source_rgb(0,0,0,)
        context.move_to(self.x,self.y) 
        context.line_to(self.x+length, self.y)
        context.stroke()
        context.set_source_rgba(0,0,0,1)
        context.move_to(self.x + 1, self.y+20)
        context.show_text(name)
        context.stroke()
    def draw_exon(self,start,end):
        '''Draw the exons using start, end, and length of exon.'''
        context.set_line_width(10)
        context.set_source_rgb(0,0,0)
        context.move_to(self.x + start, self.y + 4) # x is the beginning of the sequence (this is hard coded since we don't want x to start at the very beginning of the page so the start of exon = x + start position of exon)
        context.line_to(self.x + end, self.y+ 4) # y + 4 since I want the exon to be annotated below the sequence line
        context.stroke() 
    def draw_motif(self,start,end, color):
        '''Draw the motifs using start, end and color that is distinct to each motif '''
        context.set_line_width(10)
        context.set_source_rgb(*color)
        context.move_to(self.x + start, self.y - 6)
        context.line_to(self.x + end, self.y - 6) # y - 6 because I want the motifs to be annotated above the sequence line
        context.stroke()

# Iterate through motif text file and store all motifs in a list
with open(m, "r") as m_input:
    all_motifs = []
    for l in m_input:
        l = l.strip()
        all_motifs.append(l)

# Create a dictionary with each motif as key and a corresponding color
color_list = [(255, 0, 0), (0, 255, 0), (0, 0, 255), (0, 102, 102), (153, 0, 153)]
motif_colors = color_list[0:len(all_motifs)]
color_dict = {}
for i in range(0,len(all_motifs)):
    if not color_dict.get(all_motifs[i]):
        color_dict[all_motifs[i]] = motif_colors[i]

# MAIN code: iterate through FASTA file and draw gene name, sequences, exons, and motifs
with open(f'{f}_oneline',"r") as input:
    motif_searches=[]
    lengths = []
    i = 20
    for line in input:
        line = line.strip()
        if line[0] == ">" :
            header = line.split()
            gene_name = header[0][1:] # Only want the motif name as label
        else:
            sequence = DNAsequence(line)
            lengths.append(len(line))
            result = Drawing(5, i) 
            result.draw_gene(gene_name,len(line)) # draw the sequence
            result.draw_exon(sequence.get_exon_start(), sequence.get_exon_end()) #draw the exons
            motifs = Transcript(sequence.sequence,all_motifs) # Initate class Transcript
            found_motifs = motifs.find_motif() # found_motifs is a dictionary containing motif as keys and list of (start,end) positions as values
            for motif,positions in found_motifs.items():
                for position in positions:
                    start = position[0]
                    end = position [1]
                    result.draw_motif(start,end, color_dict[motif])
        i += 40

# Define legend labels and colors
color_dict["Exon"] = (0,0,0)
# Draw legend
legend_x = 50
legend_y = height - (len(color_dict)*21) 
legend_spacing = 20
legend_font_size = 12
for label, color in color_dict.items():
    # Draw colored square
    context.set_source_rgb(*color)
    context.rectangle(legend_x, legend_y, 10, 10)
    context.fill()
    # Draw legend label
    context.set_source_rgb(0, 0, 0)  # Set text color to black
    context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(legend_font_size)
    context.move_to(legend_x + 15, legend_y + 8)
    context.show_text(label)
    # Move to the next legend item
    legend_y += legend_spacing
# Save the surface to a PNG file
output_name = f.split(".")[0]
surface.write_to_png(f'{output_name}.png')
