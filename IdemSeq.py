#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Takes a nucleic acid sequence and converts it into another, encoding the same
# amino acids sequence.
#
# For now, sequences have to start with ATG and finish with TAA/TAG/TGA.
#
# We assert that the same amino acids sequence is encoded, but there's no real
# error verifications. No warranty of actually working.
#
# Example use:
#
# >>> from IdemSeq import idemseq
# >>> seq = 'ATGATCGATCGATCGATCTAG'
# >>> idemseqIseq)
# # Output: 'ATGATCGATAGGAGTATATGA'
# >>> idemseqIseq)
# # Output: 'ATGATAGATCGCTCCATATAA'

# aa  --> Amino acids
# na  --> Nucleic acids
# seq --> Sequence
# cod --> Codons

#from pylab import *
from random import choice

# Used to map the codons to the amino acids then encode
aas_map = dict(
    F=["TTT","TTC"],
    L=["TTA", "TTG", "CTT","CTC","CTA", "CTG"],
    I=["ATT", "ATC", "ATA"],
    S=["TCT", "TCC", "TCA", "TCG","AGT","AGC"],
    P=["CCT", "CCC", "CCA","CCG"],
    Y=["TAT", "TAC"],
    C=["TGT", "TGC"],
    H=["CAT", "CAC"],
    W=["TGG"],
    Q=["CAA", "CAG"],
    R=["CGT", "CGC", "CGA","CGG", "AGA","AGG"],
    T=["ACT", "ACC", "ACA","ACG"],
    N=["AAT", "AAC"],
    K=["AAA", "AAG"],
    V=["GTT", "GTC", "GTA","GTG"],
    A=["GCT", "GCC", "GCA","GCG"],
    M=["ATG"],
    G=["GGT", "GGC", "GGA","GGG"],
    D=["GAT", "GAC"],
    E=["GAA", "GAG"]
)

# Reversing the preceding dict to map amino acids to the codons encoding them
cod_map = {}
for k, v in aas_map.iteritems():
    for w in v:
        cod_map[w] = k

# Codon sequence to amino acids sequence
def codseq_to_aaseq(c):
    a = [choice(cod_map[i]) for i in c]
    return a

# Amino acid sequence to codon sequence
def aaseq_to_codseq(a):
    c = [choice(aas_map[i]) for i in a]
    return c

# Trim both ends of nucleic acids sequence in string format
def trim_ends(s):
    i = 3 if s[:3]=="ATG" else None
    j = -3 if s[-3:] in ["TAA", "TAG", "TGA"] else None
    return s[i:j]

# Put some ends back to our sequence of nucleic acids in string format
def add_ends(s):
    return "ATG"+s+choice(["TAA", "TAG", "TGA"])

# Sequence of nucleic acids in string format -->  list of codons
def string_to_seq(s):
    n=3
    tmp = trim_ends(s)
    return [tmp[i:i+3] for i in range(0,len(tmp),n)]

# List of codons --> sequence of nucleic acids in string format
def seq_to_string(s):
    o = ''
    for i in s:
        o+=i
    return add_ends(o)

# From a NA sequence in string format, generates a new one (idem AA sequence)
def idemseq(seqstring):
    codseq = string_to_seq(seqstring)
    aaseq = codseq_to_aaseq(codseq)
    tmp = aaseq_to_codseq(aaseq)
    # Ensures it's the same amino acids sequence
    assert(codseq_to_aaseq(tmp)==codseq_to_aaseq(codseq))
    return seq_to_string(tmp)
    
