#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *01.04.2015
"""

import sys, os
from optparse import OptionParser
import pysam
import string
from collections import defaultdict

table = string.maketrans('TGCA','ACGT') # COMPLEMENT DNA

parser = OptionParser()
parser.add_option("-f","--fasta", dest="reference", help="Fasta index reference genome (default /net/shendure/vol1/home/mkircher/sequencedb/genome/grch37_1000g_phase2/whole_genome.fa)",default="/net/shendure/vol1/home/mkircher/sequencedb/genome/grch37_1000g_phase2/whole_genome.fa")
parser.add_option("-m","--max", dest="maxreads", help="Maximum number of reads to consider (default 0 = Off)",default=0,type="int")
parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
parser.add_option("-o","--outfile", dest="outfile", help="Write output to file (def STDOUT)",default="")
parser.add_option("-w","--window", dest="window", help="Position extension in both directions (def 150)",type="int",default=150)
(options, args) = parser.parse_args()

if not os.path.exists(options.reference) or not os.path.exists(options.reference+".fai"):
  sys.stderr.write("Fasta indexed reference genome is not available.\n")
  sys.exit()

reference = pysam.Fastafile(options.reference)

data = {'A':defaultdict(int),'C':defaultdict(int),'G':defaultdict(int),'T':defaultdict(int)}
dinucleotides = {}
for base1 in ['A','C','G','T']:
  for base2 in ['A','C','G','T']:
    dinucleotides[base1+base2]=defaultdict(int)

for line in sys.stdin:
  fields = line.split()
  chrom = fields[0].replace("chr","")
  pos = int(fields[1])

  seq = ""
  try:
    seq = reference.fetch(chrom,pos-options.window,pos+options.window+2)
  except:
    continue
      
  for ind,base in enumerate(seq[:-1]):
    try: data[base][ind]+=1
    except: pass

    try: dinucleotides[base+seq[ind+1]][ind]+=1
    except: pass


  if (options.maxreads > 0) and (readCount > options.maxreads): 
    sys.stderr.write("Maximum number of reads reached.\n")
    break 

if options.outfile != "":
  outfile = open(options.outfile,"w")
else:
  outfile = sys.stdout

outfile.write("Base\t"+"\t".join(map(str,range(-1*options.window,0)+range(0,options.window+1)))+"\n")
helper = map(lambda (x,y):y,sorted(data['A'].iteritems()))
for base in ['C','G','T']:
  for ind,val in data[base].iteritems():
    helper[ind]+=val

for base,subDict in sorted(data.iteritems()):
  outfile.write(base+"\t"+"\t".join(map(lambda (ind,x): "NA" if helper[ind] == 0 else "%.4f"%(x/float(helper[ind])),sorted(subDict.iteritems())))+"\n")

helper = map(lambda (x,y):y,sorted(dinucleotides['AA'].iteritems()))
for base in dinucleotides.keys():
  if base == "AA": continue
  for ind,val in dinucleotides[base].iteritems():
    helper[ind]+=val

for base,subDict in sorted(dinucleotides.iteritems()):
  outfile.write(base+"\t"+"\t".join(map(lambda (ind,x): "NA" if helper[ind] == 0 else "%.4f"%(x/float(helper[ind])),sorted(subDict.iteritems())))+"\n")

if (options.outfile != ""):
  outfile.close()