#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *21.07.2016
"""

import sys, os
from optparse import OptionParser
import gzip
import pysam
import math
from collections import defaultdict
from AnalysisLib import get_from_tabix,eval_1000G_frequencies
from bx.intervals.intersection import Intersecter, Interval

genomeBuild = "GRCh37"
#commonVars = set(['rs6048','rs6049','rs1800291','rs1800292','rs1800297','rs1050705','rs1396947','rs440051'])


def prefix(alleles):
  if len(alleles) > 1:
    check_shared = alleles[0]
    while len(check_shared) > 0:
      if reduce(lambda x,y: x and y.startswith(check_shared),alleles,True):
        return check_shared
      else: 
        check_shared = check_shared[:-1]
    return ""
  else:
    return ""


def sharedPrefix(s1,s2):
  minLength = min(len(s1),len(s2))
  shared = 0
  for ind in range(minLength-1):
    if s1[ind] == s2[ind]: shared+=1
    else: break
  if minLength == 1:
    return max(0,shared-1)
  else:
    return shared


def sharedSuffix(s1,s2):
  minLength = min(len(s1),len(s2))-1
  shared = 0
  for ind in range(minLength*-1,0)[::-1]:
    if s1[ind] == s2[ind]: shared+=1
    else: break
  return shared 


def splitFields(x):
  helper = x.partition("=")
  return helper[0],helper[2]


def eval_sex_check(filename):
  res = {}
  infile = open(filename)
  for line in infile:
    fields = line.split()
    if len(fields) == 3:
      sample = fields[0]
      if sample.endswith(".bam"): sample = fields[0][:-4]
      if sample.endswith(".M"): sample = ".".join(sample.split(".")[:-1])
      total = int(fields[2])
      sry = int(fields[1])
      state = "?"
      info = "SRY/Total: %d/%d = %.2f%%"%(sry,total,0 if total == 0 else sry/float(total)*100)
      if (total > 0) and (sry > total*0.001):
        state = 'M'
      elif (sry == 0) and (total > 1000):
        state = 'F'
      elif (total > 1000) and (sry < total*0.0001):
        state = 'F?'
      res[sample] = state,info
  infile.close()
  return res


def median(vals):
  sorted_values = list(vals)
  sorted_values.sort()
  if len(sorted_values) == 0: return None
  elif len(sorted_values) % 2 == 0:
    return (sorted_values[len(sorted_values)//2-1]+sorted_values[len(sorted_values)//2])*0.5
  else:
    return sorted_values[(len(sorted_values)-1)//2]


def percentile(vals,percentile):
  sorted_values = list(vals)
  sorted_values.sort()
  if len(sorted_values) == 0 or (0 > percentile) or (percentile > 1.0): return None
  else:
    return sorted_values[min(int(round((len(sorted_values)-1)*percentile)),len(sorted_values)-1)]


parser = OptionParser("%prog [options]")
parser.add_option("--vcf", dest="vcf", help="Filename of input multi-sample VCF file with all sites (def 'realign_all_samples.all_sites.vcf.gz')",default="realign_all_samples.all_sites.vcf.gz")
parser.add_option("--vep", dest="vep", help="VEP results (def 'realign_all_samples.vep.tsv.gz')",default="realign_all_samples.vep.tsv.gz")
parser.add_option("-i","--inversions", dest="inversions", help="Analysis results for inversion MIPs (def 'inversion_mips/inversion_summary_counts.txt')",default="inversion_mips/inversion_summary_counts.txt")
parser.add_option("-s","--sample_sex", dest="sample_sex", help="Analysis results for sex check (def 'samples_sex_check.txt')",default="samples_sex_check.txt")
parser.add_option("-t","--target", dest="target", help="BED file of target regions (def 'target_coords.bed')",default="target_coords.bed")
parser.add_option("-f", "--factor", dest="factor", help="Allowed deviation for MIP performance (def 10)",default=10,type="int")
parser.add_option("-m","--mipstats", dest="mipstats", help="File with MIP performance counts (def 'realign_all_samples.MIPstats.tsv')",default="realign_all_samples.MIPstats.tsv")
parser.add_option("-c","--indelCheck", dest="indelCheck", help="Only report indels with count evidence (def 'realign_all_samples.indel_check.txt')",default="realign_all_samples.indel_check.txt")
parser.add_option("-d", "--design", dest="design", help="MIP design file (default hemomips_design.txt)",default="hemomips_design.txt")
parser.add_option("--TG", dest="TG", help="1000 Genomes variant tabix file" )
parser.add_option("-b", "--benign",dest="benign", help="List of benign variants" )
#parser.add_option("--freq", dest="freq", help="Maximum 1000 Genomes allele frequency (def 0.05)",type="float",default=0.05)
(options, args) = parser.parse_args()

#benignVars = set(['X_138633280_A/G','X_138623355_A/G','X_154158285_G/C','X_154158201_T/G','X_154064200_C/T','X_154064580_T/C','X_138644917_G/A','X_154194886_C/T','X_154159851_G/A','X_154159104_C/T','X_154158444_G/A','X_154157565_C/T','X_154157330_T/C','X_154132301_G/T','X_154088758_G/A','X_154065843_G/A','X_154065794_G/A','X_154065446_C/T','X_154065069_T/G','X_138642995_T/C','X_138643939_A/G','X_138644836_G/A','X_138645058_GT/-','X_138645060_-/GT','X_138645149_T/C','X_138645157_G/C','X_154088838_T/C','X_154221432_G/A'])

benignVars = set()

if os.path.exists(options.benign):
  infile = open(options.benign)
  for line in infile:
    benignVars.add(line.rstrip()) 
  infile.close()

print benignVars


sex_check = eval_sex_check(options.sample_sex)

inversion_names = ["inv22_ID+IU","inv22_ED+2U","inv22_ED+3U","inv22_ID+2U","inv22_ID+3U","inv22_ED+IU","inv1_1IU+1ID","inv1_1IU+1ED"]
INT22_inversion_types = [
   ("INT22-1#1"  , [False,True,False,False,True,False,  None,None] ),
   ("INT22-1#2"  , [False,True,False,False,False,True,  None,None] ),
   ("INT22-1#3"  , [False,True,False,False,True,True,   None,None] ),
   ("INT22-1#4"  , [False,False,False,False,True,True,  None,None] ),
   ("INT22-2#5"  , [False,False,True,True,False,False,  None,None] ),
   ("INT22-2#6"  , [False,False,True,True,False,True,   None,None] ),
   ("INT22-2#7"  , [False,False,True,False,False,True,  None,None] ),
   ("INT22-2#8"  , [False,False,False,True,False,True,  None,None] ),
("INT22-unknown" , [False,False,False,False,False,True, None,None] )  ]

noINT22_inversion_types = [
  ("benign_dup" , [True,True,True,False,False,True,    None,None] ),
  ("noINT22#1" , [True,True,True,False,False,False,    None,None] ), 
  ("noINT22#2" , [True,True,False,False,False,False,   None,None] ), 
  ("noINT22#3" , [True,False,False,False,False,False,  None,None] ), 
  ("noINT22#4" , [True,False,True,False,False,False,   None,None] ), 
  ("noINT22#5" , [False,True,True,False,False,False,   None,None] )  ]

INT22failed_inversion_types = [
  ("INT22-FAILED#1" , [False,True,False,False,False,False,  None,None] ),
  ("INT22-FAILED#2" , [False,False,True,False,False,False,  None,None] ),
  ("INT22-FAILED#3" , [False,False,False,True,False,False,  None,None] ),
  ("INT22-FAILED#4" , [False,False,False,False,True,False,  None,None] ),
  ("INT22-FAILED#5" , [False,False,False,False,False,False, None,None] )  ]

INT1_inversion_types = [
  ("INT1"      , [None,None,None,None,None,None,      False,True] ),
  ("noINT1"    , [None,None,None,None,None,None,      True,False] ),
("INT1-FAILED" , [None,None,None,None,None,None,      False,False] ),
("Conflict: INT1" , [None,None,None,None,None,None,   True,True] )   ]

if not os.path.exists(options.vep) or not os.path.exists(options.vcf):
  sys.stderr.write("Error: VEP and/or VCF input files not available!\n")
  sys.exit()

TGTabix = None
if not os.path.exists(options.TG+".tbi"):
  sys.stderr.write("1000 Genomes tabix: Require valid path to compressed tabix file and tabix index file.\n")
  sys.exit()
else:
  sys.stderr.write('1000 Genomes variants tabix file (%s)...\n'%(options.TG))
  TGTabix = pysam.Tabixfile(options.TG,'r'),None,None,None,"1000 Genomes variants"

bedanno = None
coverage_stats_by_region = {}
sorted_regions = []
name2region = {}
if options.target != "" and os.path.exists(options.target):
  bedanno = {}
  infile = open(options.target)
  for line in infile:
    fields = line.rstrip().split('\t')
    if len(fields) > 3:
      chrom = fields[0]
      start = int(fields[1])
      end = int(fields[2])
      name = fields[3]
      if chrom not in bedanno: bedanno[chrom] = Intersecter()
      bedanno[chrom].add_interval( Interval(start+1, end+1, value = name) )
      sorted_regions.append(name)
      name2region[name] = (chrom,start,end)
      if name not in coverage_stats_by_region:
        coverage_stats_by_region[name] = defaultdict(int)
  infile.close()
  
inversion_obs = {}
if os.path.exists(options.inversions):
  infile = open(options.inversions)
  for line in infile:
    fields = line.split()
    if len(fields) >= 1:
      individual = fields[0]
      if individual.endswith(".M"): individual = ".".join(individual.split(".")[:-1])
      inversion_obs[individual] = defaultdict(int)
      for counts in fields[1:]:
        count,name = counts.split(':')
        count = int(count)
        inversion_obs[individual][name]=count
  infile.close()


MIPcoords = {}
if os.path.exists(options.design):
  # DESIGN FILE ARM RANGES ARE 1-BASED, CLOSED INTERVALS
  fchrom,ligStart,ligEnd,extStart,extEnd,fstrand,fmipname = 2,7,8,3,4,17,-1
  infile = open(options.design)
  for line in infile:
    fields = line.rstrip().split("\t")
    if len(fields) > 8:
      if line.startswith(">") or line.startswith('#'):
        for ind,elem in enumerate(fields):
          if elem == "chr" or elem == "chrom": fchrom = ind
          elif elem == "ext_probe_start": extStart = ind
          elif elem == "ext_probe_stop": extEnd = ind
          elif elem == "lig_probe_start": ligStart = ind
          elif elem == "lig_probe_stop": ligEnd = ind
          elif elem == "probe_strand": fstrand = ind
          elif elem == "mip_name": fmipname = ind
      else:
        chrom,lstart,lend,estart,eend,strand,mipname = fields[fchrom],int(fields[ligStart])-1,int(fields[ligEnd]),int(fields[extStart])-1,int(fields[extEnd]),fields[fstrand],fields[fmipname]
        if strand == "+": MIPcoords[mipname] = chrom,lend,estart+1
        else: MIPcoords[mipname] = chrom,lstart,eend
  infile.close()
else:
  sys.stderr.write("MIP design file (%s) not available.\n"%(options.design))
  sys.exit()

allsamples = []
sample2ind = {}
TotalSample = []
MIPcounts = {}
failedMIPs = defaultdict(list)
failedMIPs_summary = defaultdict(list)

if os.path.exists(options.mipstats):
  ######################
  # Reading MIP counts #
  ######################

  infile = open(options.mipstats)
  for line in infile:
    if line.startswith("#"):
      allsamples = map(lambda x: x if not x.endswith(".M") else ".".join(x.split(".")[:-1]),line.rstrip().split("\t")[1:])
      TotalSample = len(allsamples)*[0]
      for ind,sample in enumerate(allsamples):
        sample2ind[sample] = ind
      #print allsamples
      #print TotalSample
    else:
      fields = line.rstrip().split("\t")
      if len(fields) == len(allsamples)+1:
        MIPcounts[fields[0]] = map(int,fields[1:])
        for ind,count in enumerate(MIPcounts[fields[0]]):
          TotalSample[ind]+=count
      #else:
        #print len(fields), len(allsamples)+1
  infile.close()

  ######################################
  # By-plate MIP performance analysis  #
  ######################################

  plates = set()
  for i in allsamples:
    plates.add("_".join(i.split("_")[:2]))
  plates = list(plates)
  plates.sort()

  for plate in plates:
    for mip,counts in sorted(MIPcounts.iteritems()):
      if mip not in MIPcoords: continue
      if mip.startswith("Y"): continue 

      vals = []
      psamples = []
      for ind,count in enumerate(counts):
        if allsamples[ind].startswith(plate):
          vals.append(count/float(TotalSample[ind]))
          psamples.append(allsamples[ind])

      ## Infer variance across samples
      med,lowsig,uppersig = median(vals),percentile(vals,0.341),percentile(vals,0.682)
      if med == 0: continue
      #sig = ((med-lowsig) + (uppersig-med))*0.5
      siglow,sighigh = med-lowsig, uppersig-med
      factor = options.factor # qnorm(.975) = 1.959964, qnorm(.995) = 2.575829
      #intlow,inthigh = max(0.0,med-sig*factor),min(med+sig*factor,1.0)
      intlow,inthigh = max(0.0,med-siglow*factor),min(med+sighigh*factor,1.0)
      outliers = 0
      for ind,individual in enumerate(psamples):
        if (intlow > vals[ind]) or (vals[ind] > inthigh): 
          outliers+=1
          if vals[ind] > inthigh: 
            failedMIPs[individual].append("+:"+mip)
            failedMIPs_summary["+:"+mip].append(individual)
          elif vals[ind] < intlow: 
            failedMIPs[individual].append("-:"+mip)
            failedMIPs_summary["-:"+mip].append(individual)
else:
  sys.stderr.write("Error: MIP stats file (%s) not available.\n"%(options.design))

GT_stats = {}
GT_stats_gene = defaultdict(int)
genes = set()
coverage_stats = {}
coverage_holes = {}
het_stats = {}
var_stats = {}
variants = {}

VEP = {}

if options.vep.endswith('.gz'):
  infile = gzip.open(options.vep)
else:
  infile = open(options.vep)
VEPheader = None
for line in infile:
  if line.startswith('##'): continue
  elif line.startswith('#Chrom'):
    VEPheader = line[1:].rstrip().split("\t")
    VEPheader_html = map(lambda x: x.replace("_","<br>"),VEPheader)
  else:
    if VEPheader == None: 
      sys.stderr.write("Error: VEP file misses header.\n")
      sys.exit()
    else:
      fields = line.rstrip().replace("%3D","=").split("\t")
      vepline = dict(zip(VEPheader,fields))
      if 'Uploaded_variation' in vepline and genomeBuild+"_"+vepline['Uploaded_variation'] not in VEP:
        chrom,pos,alleles = vepline['Uploaded_variation'].split("_")
        pos = int(pos)
        alleles =alleles.split("/")
        TGTabix,variantLines = get_from_tabix(TGTabix,chrom,pos)
        #print variantLines,chrom,pos,alleles
        F1000g,ASN_AF,AMR_AF,AFR_AF,EUR_AF,isLowCov = eval_1000G_frequencies(variantLines,alleles[0],alleles[-1])
        #print F1000g,vepline['Uploaded_variation']
        #if F1000g <= options.freq:
        if F1000g > 0: 
          fields[-1]+=';1000G_AF=%.5f'%F1000g
          fields[-1]=fields[-1].lstrip(";")
        isBenign = False
        if vepline['Uploaded_variation'] in benignVars: 
          isBenign = True
        VEP[genomeBuild+"_"+vepline['Uploaded_variation']] = fields,isBenign
      else:
        sys.stderr.write("Error: Unexpected duplication of annotation line or misformed line in VEP.\n")
        sys.exit()
infile.close()

fchrom = 0
fpos = 1
fname = 2
fref = 3
falt = 4
fqual = 6
finfo = 7
fformat = 8
findividual = 9

sites = set()
sites_gene = defaultdict(set)

InDelCheck_counts = {}
if os.path.exists(options.indelCheck):
  infile = open(options.indelCheck)
  for line in infile:
    fields = line.split()
    fdict = dict(map(lambda x: (":".join(x.split(":")[1:]),int(x.split(":")[0])),map(lambda x: x if not x.endswith(".M") else ".".join(x.split(".")[:-1]),fields[1:])))
    InDelCheck_counts["%s_%s"%(genomeBuild,fields[0])] = fdict
  infile.close()

if options.vcf.endswith('.gz'):
  infile = gzip.open(options.vcf)
else:
  infile = open(options.vcf)
header = None
individuals = []
is_gatk4 = False
varcalls_gatk4 = {}
for line in infile:
  if line.startswith('##'): continue
  elif line.startswith('#CHROM'):
    header = line[1:].rstrip().split("\t")
    individuals = map(lambda x: x if not x.endswith(".M") else ".".join(x.split(".")[:-1]), header[findividual:])
    for individual in individuals:
      GT_stats[individual] = 0
      coverage_stats[individual] = 0
      coverage_holes[individual] = []
      het_stats[individual] = 0
      var_stats[individual] = 0
      variants[individual] = []
  else:
    if header == None: 
      sys.stderr.write("Error: VCF file misses header.\n")
      sys.exit()
    else:
      fields = line.rstrip().split("\t")
      chrom,pos = fields[0],int(fields[1])
      if bedanno != None: 
        if chrom not in bedanno: continue
        region_name = None
        for cinterval in bedanno[chrom].find(pos,pos+1):
          region_name = cinterval.value
        if region_name == None: continue
        gene_name = region_name.split("/")[0]
        genes.add(gene_name)
        
      #VCFline = dict(zip(header,fields))
      if (fields[falt].endswith('<NON_REF>')):
        is_gatk4 = True
        alleles = [fields[fref]]+fields[falt].split(',')[:-1]
      elif (fields[falt] != '.'):
        alleles = [fields[fref]]+fields[falt].split(',')
      else:
        alleles = [fields[fref]]
      allele_dict = dict(map(lambda (x,y):(str(x),y), enumerate(alleles)))
      
      if is_gatk4 and len(varcalls_gatk4) == 0:
        finalVariantsFilename = options.vcf.replace(".all_sites","")
        if os.path.exists(finalVariantsFilename) and options.vcf.endswith('.gz'):
          finalVars = gzip.open(finalVariantsFilename)
        elif os.path.exists(finalVariantsFilename):
          finalVars = open(finalVariantsFilename)
        for vline in finalVars:
          if vline.startswith('#'): continue
          vfields = vline.rstrip().split("\t")
          varcalls_gatk4[(vfields[fchrom],vfields[fpos],vfields[fref],vfields[falt].replace(",<NON_REF>",""))]=vfields
        finalVars.close()
        if (len(varcalls_gatk4) > 0):
          sys.stderr.write("Read %d variants from filtered variant output file (%s, assuming GATK 4).\n"%(len(varcalls_gatk4),finalVariantsFilename))
      
      count_GT = 0
      sample_GTs = []
      sample_coverages = []
      formatfields = fields[fformat].split(':')
      for ind,individual in enumerate(individuals):
        values = dict(zip(formatfields,fields[ind+findividual].split(':')))
        if (fields[4] == '<NON_REF>') or (fields[4].endswith('<NON_REF>') and (fields[fchrom],fields[fpos],fields[fref],fields[falt].replace(",<NON_REF>","")) not in varcalls_gatk4): # HOMOZYGOTE REFERENCE CALLS GATK 4
            if ('DP' in values) and (int(values['DP']) >= 3):
              sample_GTs.append(1)
              count_GT += 1
              sample_coverages.append(int(values['DP']))
            else:
              sample_GTs.append(0)
              sample_coverages.append(0)
        elif ((fields[fchrom],fields[fpos],fields[fref],fields[falt].replace(",<NON_REF>","")) in varcalls_gatk4) or fields[ind+findividual] != "./.": # COMPATIBILITY WITH GATK 3
          callQual = fields[fqual]
          if (fields[fchrom],fields[fpos],fields[fref],fields[falt].replace(",<NON_REF>","")) in varcalls_gatk4:
            vfields = varcalls_gatk4[(fields[fchrom],fields[fpos],fields[fref],fields[falt].replace(",<NON_REF>",""))]
            if vfields[falt].endswith("<NON_REF>"):
              alleles = [vfields[fref]]+vfields[falt].split(',')[:-1]
            else:
              alleles = [vfields[fref]]+vfields[falt].split(',')
            allele_dict = dict(map(lambda (x,y):(str(x),y), enumerate(alleles)))
            callQual = vfields[fqual]
            values = dict(zip(vfields[fformat].split(":"),vfields[ind+findividual].split(':')))
          if (('AD' in values) and ('GT' in values) and (sum(map(lambda x: 0 if not x.isdigit() else int(x),values['AD'].split(","))) >= 3)) or (('AD' not in values) and ('GT' in values) and ('DP' in values) and values['DP'].isdigit() and (int(values['DP']) >= 3)):
            obsalleles = []
            sample_GTs.append(1)
            for gt in values['GT'].split('/'):
              if gt in allele_dict: 
                obsalleles.append(allele_dict[gt])
            if len(obsalleles) == 2:
              count_GT += 1
            if ('GQ' in values) and (values['GQ'] != ".") and (float(values['GQ']) >= 30) and ('DP' in values) and (int(values['DP']) >= 8) and (len(set(obsalleles)) != 1): 
              het_stats[individual]+=1

            for alt in set(obsalleles):
              if alt != fields[fref]: 
                ppos,pref,palt = fields[fpos],fields[fref],alt
                if not(len(palt) == len(pref) and len(pref) == 1):
                  trimValue = sharedPrefix(pref,palt)
                  if trimValue != 0:
                    pref = pref[trimValue:]
                    palt = palt[trimValue:]
                    ppos = str(int(ppos)+trimValue)
                  trimValue2 = sharedSuffix(pref,palt)
                  if trimValue2 != 0:
                    pref = pref[:-trimValue2]
                    palt = palt[:-trimValue2]
                    
                pfix = len(prefix([pref,palt]))
                varstring = "%s_%s_%s_%s/%s"%(genomeBuild,fields[fchrom],ppos,pref,palt)
                varObs = None if varstring not in InDelCheck_counts else (0 if individual not in InDelCheck_counts[varstring] else InDelCheck_counts[varstring][individual])
                
                if varstring not in VEP:
                  varstring = "%s_%s_%d_%s/%s"%(genomeBuild,fields[fchrom],int(ppos)+pfix,pref[pfix:],palt[pfix:])
                if varstring not in VEP:
                  varstring = "%s_%s_%d_%s/%s"%(genomeBuild,fields[fchrom],int(ppos)+pfix,"-",palt[pfix:])
                if varstring not in VEP: 
                  varstring = "%s_%s_%d_%s/%s"%(genomeBuild,fields[fchrom],int(ppos)+pfix,pref[pfix:],"-")
                if varstring not in VEP:
                  sys.stderr.write("Can not retrieve variant from VEP output %s_%s_%s/%s\n"%(genomeBuild,fields[fchrom],fields[fpos],pref,palt))
                  continue
                
                if varObs != None and float(varObs)/int(values['DP']) < 0.05: continue
                
                if (('GQ' in values) and (float(values['GQ']) < 30)) or (('DP' in values) and (int(values['DP']) < 8)):
                  callQual = "LowQual"
                else:
                  if callQual == ".": callQual = "OK"
                  var_stats[individual]+=1
                
                variants[individual].append((varstring,values['GT'],values['GQ'],values['AD'],values['DP'],callQual))
          else:
            sample_GTs.append(0)
          if ('DP' in values) and (int(values['DP']) >= 3):
            sample_coverages.append(int(values['DP']))
          else:
            sample_coverages.append(0)
        else:
          sample_GTs.append(0)
          sample_coverages.append(0)
        
      if count_GT > len(individuals)//2 and (chrom,pos) not in sites:
        sites.add((chrom,pos))
        sites_gene[gene_name].add((chrom,pos))
        coverage_stats_by_region[region_name]['Count']+=1
        for ind,individual in enumerate(individuals):
          GT_stats[individual]+=sample_GTs[ind]
          GT_stats_gene[(gene_name,individual)]+=sample_GTs[ind]
          coverage_stats[individual]+=sample_coverages[ind]
          coverage_stats_by_region[region_name][individual+':Cov']+=sample_coverages[ind]
          coverage_stats_by_region[region_name][individual+':GT']+=sample_GTs[ind]
          if sample_GTs[ind] == 0:
            if len(coverage_holes[individual]) > 0:
              last = coverage_holes[individual][-1]
              if last[0] == chrom and last[2]+1 == pos:
                coverage_holes[individual][-1]=(chrom,last[1],pos)
              else:
                coverage_holes[individual].append((chrom,pos,pos))
            else:
              coverage_holes[individual].append((chrom,pos,pos))
infile.close()

total_sites = len(sites)
total_sites_genes = dict(map(lambda x: (x,len(sites_gene[x])),sites_gene.keys()))

try:
  os.makedirs('report')
except:
  pass

outfile3 = open('report/report.html','w')
outfile3.write("<html>\n <head>\n  <title>Sample Summary Report</title>\n </head>\n<body>\n")
outfile3.write("""<p align="center"><h1>Sample Summary Report</h1></p>\n""")
outfile3.write("""<table cellpadding="5" border="3">\n""")
outfile3.write("""<tr><th>Individual</th><th>Sex</th><th>Short variants</th><th>Incomplete coverage</th><th>Deletions<br>(<50% covered)</th><th>INT1</th><th>INT22</th><th>Status</th></tr>\n""")

outfile = open('report/summary.html','w')
outfile.write("<html>\n <head>\n  <title>Summary Report</title>\n </head>\n<body>\n")
outfile.write("""<p align="center"><h1>Summary Report</h1></p>\n""")
outfile.write("""<p>Total number of samples: <strong>%d</strong></p>\n"""%(len(individuals)))
outfile.write("""<p>Total number of sites considered [cov in >50%% samples]: <strong>%d</strong></p>\n"""%(total_sites))
outfile.write("""<p></p>\n""")
outfile.write("""<p><h2><a href="report.html">Sample summary</a></h2></p>\n""")
outfile.write("""<table cellpadding="5" border="3">\n""")
outfile.write("""<tr><th>SampleID</th><th>Sex</th><th>GTs</th><th>%GTs</th><th>Ave.Cov</th><th>Hets</th><th>Variants</th><th>VariantList (incl. low quality)</th></tr>\n""")

#sorted_regions = coverage_stats_by_region.keys()
#sorted_regions.sort()

variantStats = {}
coveragRegionStats = {}
inversionMipStats = {}

tbl_out1 = open('report/ind_status.csv','w')
tbl1_header = ["TubeID", "SampleID", "Sex", "SRY/Total", "Number of over-/under-performing MIPs", "Performance outlier MIPs", "Number of Sites with GTs", "Percent sites with GT"]
for gene in sorted(genes):
  tbl1_header.append("%s-Number of Sites with GTs"%gene)
  tbl1_header.append("%s-Percent sites with GT"%gene)
tbl1_header = tbl1_header+["Average Coverage", "Average Coverage of GT Calls", "Number Hets called", "Number short variants called", "Incomplete Coverage", "Deletions (<50% Covered)", "Status", "Notes"] #  "Well", "Plate"
tbl_out1.write('"%s"\r\n'%('","'.join(tbl1_header)))

tbl_out2 = open('report/inversion_calls.csv','w')
tbl2_header = ["ID", "SampleID", "Inversion MIP Reads", "Inversion Results", "Status"]
tbl_out2.write('"%s"\r\n'%('","'.join(tbl2_header)))

tbl_out3 = open('report/variant_calls.csv','w')
tbl3_header = ["ID", "SampleID", "Location", "GT", "GQ", "AD", "DP", "Status"]
tbl_out3.write('"%s"\r\n'%('","'.join(tbl3_header)))

tbl_out3_benign = open('report/variant_calls_benign.csv','w')
tbl_out3_benign.write('"%s"\r\n'%('","'.join(tbl3_header)))

tbl_out4 = open('report/variant_annotation.csv','w')
tbl4_header = ["Location", "Build", "Chrom", "Pos", "Ref", "Alt", "Gene", "Region", "cDNA", "CDS", "Protein", "HGVS Transcript", "HGVS Protein", "rsID", "AF1000G", "Notes"]
tbl_out4.write('"%s"\r\n'%('","'.join(tbl4_header)))

for location,(vepline,isBenign) in VEP.iteritems():
  build,chrom,pos,ref_alt = location.split("_")
  ref,alt = ref_alt.split('/')
  tbl4_fields = [location,build,chrom,pos,ref,alt]
  
  helper = dict(map(lambda x:splitFields(x),vepline[-1].split(";")))
  tbl4_fields.append(helper["SYMBOL"])
  
  if "EXON" in helper:
    tbl4_fields.append("Exon %s"%(helper["EXON"]))
  elif "INTRON" in helper:
    tbl4_fields.append("Intron %s"%(helper["INTRON"]))
  else:
    tbl4_fields.append("unknown")

  if vepline[10] != "-": tbl4_fields.append(vepline[10])
  else: tbl4_fields.append('')
  
  if vepline[11] != "-": tbl4_fields.append(vepline[11])
  else: tbl4_fields.append('')

  if vepline[12] != "-": tbl4_fields.append(vepline[12])
  else: tbl4_fields.append('')

  if ("HGVSc" in helper): tbl4_fields.append(":".join(helper["HGVSc"].split(":")[1:]))
  else: tbl4_fields.append("")
  if ("HGVSp" in helper): tbl4_fields.append(":".join(helper["HGVSp"].split(":")[1:]))
  else: tbl4_fields.append("")
  
  rsIDs = []
  for cID in vepline[-2].split(','):
    if cID.startswith('rs'): rsIDs.append(cID)
  tbl4_fields.append(','.join(rsIDs))
  
  if ("1000G_AF" in helper): tbl4_fields.append(helper["1000G_AF"])
  else: tbl4_fields.append("")
  
  tbl4_fields.append("")
  tbl_out4.write('"%s"\r\n'%('","'.join(tbl4_fields)))
tbl_out4.close()

for individual in individuals:
  status = None # None - OK, True - check, False - failed
  status_flags = set()

  tbl1_fields = ['']
  outfile2 = open('report/ind_%s.html'%(individual),'w')
  tbl1_fields.append(individual)
  #well = individual.split("_")[2].split(".")[0]
  #well = well[-1]+well[:-1]
  #tbl1_fields.append(well)
  #tbl1_fields.append(str(int(individual.split("_")[1])))
  outfile2.write("<html>\n <head>\n  <title>Report for %s</title>\n </head>\n<body>\n"%individual)
  outfile2.write("""<p align="center"><h1>Report for %s</h1></p>\n"""%individual)
  long_sex = "%s (%s)"%(sex_check[individual][0],sex_check[individual][1])
  tbl1_fields.append(sex_check[individual][0])
  tbl1_fields.append(sex_check[individual][1].split()[1])
  outfile2.write("""<p>Sex determined from SRY MIPs: <strong>%s</strong></p>\n"""%long_sex)
  if sex_check[individual][0].endswith("?"): 
    status = True
    status_flags.add("sex")
  
  tobs = 0 if individual not in failedMIPs else len(failedMIPs[individual])
  if tobs > 1:
    if status == None: 
      status = True
      status_flags.add("MIPs")
    outfile2.write("""<p>Number of target region performance outlier MIPs: <font color="#FF6600"><strong>%d</strong></font></p>\n"""%(tobs))
  else:
    outfile2.write("""<p>Number of target region performance outlier MIPs: <strong>%d</strong></p>\n"""%(tobs))
  tbl1_fields.append(str(tobs))
  if len(failedMIPs[individual]) > 10:
    tbl1_fields.append(",".join(failedMIPs[individual][:10])+",...")
  else:
    tbl1_fields.append(",".join(failedMIPs[individual]))
 
  outfile2.write("""<p>Number of sites with GTs: <strong>%d (%0.2f%%)</strong></p>\n"""%(0 if total_sites == 0 else GT_stats[individual],0 if total_sites == 0 else GT_stats[individual]/float(total_sites)*100))
  tbl1_fields.append(str(0 if total_sites == 0 else GT_stats[individual]))
  tbl1_fields.append(str(0 if total_sites == 0 else GT_stats[individual]/float(total_sites)*100))
  for gene in sorted(genes):
    tbl1_fields.append(str(0 if total_sites_genes[gene] == 0 else GT_stats_gene[(gene,individual)]))
    tbl1_fields.append(str(0 if total_sites_genes[gene] == 0 else GT_stats_gene[(gene,individual)]/float(total_sites_genes[gene])*100))
    outfile2.write("""<p>Number of sites with GTs [%s]: <strong>%d (%0.2f%%)</strong></p>\n"""%(gene,0 if total_sites_genes[gene] == 0 else GT_stats_gene[(gene,individual)],0 if total_sites_genes[gene] == 0 else GT_stats_gene[(gene,individual)]/float(total_sites_genes[gene])*100))
    
  outfile2.write("""<p>Average coverage: <strong>%0.2f</strong></p>\n"""%(0 if total_sites == 0 else coverage_stats[individual]/float(total_sites)))
  tbl1_fields.append("%0.2f"%(0 if total_sites == 0 else coverage_stats[individual]/float(total_sites)))
  outfile2.write("""<p>Average coverage of GT calls: <strong>%0.2f</strong></p>\n"""%(0 if GT_stats[individual] == 0 else coverage_stats[individual]/float(GT_stats[individual])))
  tbl1_fields.append("%0.2f"%(0 if GT_stats[individual] == 0 else coverage_stats[individual]/float(GT_stats[individual])))
  if het_stats[individual] > 0 and sex_check[individual][0].startswith('M'):
    outfile2.write("""<p>Number of hets called: <font color="#FF6600"><strong>%d</strong></font></p>\n"""%(het_stats[individual]))
  else:
    outfile2.write("""<p>Number of hets called: <strong>%d</strong></p>\n"""%(het_stats[individual]))
  tbl1_fields.append(str(het_stats[individual]))
  outfile2.write("""<p>Number of short variants called: <strong>%d</strong></p>\n"""%(var_stats[individual]))
  tbl1_fields.append(str(var_stats[individual]))
  if len(variants[individual])-var_stats[individual] > 0:
    outfile2.write("""<p>Number of low quality variants not counted above: <font color="#FF6600"><strong>%d</strong></font></p>\n"""%(len(variants[individual])-var_stats[individual]))
  else:
    outfile2.write("""<p>Number of low quality variants not counted above: <strong>%d</strong></p>\n"""%(len(variants[individual])-var_stats[individual]))
  #tbl1_fields.append(str(len(variants[individual])-var_stats[individual]))
  outfile2.write("""<p></p>\n""")
  outfile2.write("""<p><h2>Variants identified in target region:</h2></p>\n""")

  if het_stats[individual] > 0 and sex_check[individual][0].startswith('M'):
    status = True
    status_flags.add("variants")
  
  outfile3.write("""<tr><td><a href="ind_%s.html">%s</a></td><td>%s</td>\n"""%(individual,individual,sex_check[individual][0]))
  
  o3_variants = []
  
  if len(variants[individual]) > 0:
    outfile2.write("""<table cellpadding="5" border="3">\n""")
    outfile2.write("<tr><th>"+"</th><th>".join(["GT","GQ","AD","DP","Status"]+VEPheader_html[:3]+VEPheader_html[5:])+"</th></tr>\n")
    for variant,gt,gq,ad,dp,cstatus in variants[individual]:
      tbl3_fields = ['']
      tbl3_fields.append(individual)
      tbl3_fields.append(variant)
      tbl3_fields.append(gt)
      tbl3_fields.append(gq)
      tbl3_fields.append(ad)
      tbl3_fields.append(dp)
      tbl3_fields.append(cstatus)

      if cstatus == "LowQual":
        cstatus = """<font color="#FF6600">%s</font>"""%cstatus
      helper,isBenign = list(VEP[variant][0]),VEP[variant][1]
      VEP[variant][1]
      helper[-1] = helper[-1].replace(";","<br>")
      helper[-2] = helper[-2].replace(",","<br>")


      bgcolor = ""
      if isBenign: 
        bgcolor = 'bgcolor="#f0f0f0"'
        tbl_out3_benign.write('"%s"\r\n'%('","'.join(tbl3_fields)))
      else:
        tbl_out3.write('"%s"\r\n'%('","'.join(tbl3_fields)))

      outfile2.write("<tr><td %s>"%(bgcolor)+("</td><td %s>"%(bgcolor)).join([gt,gq,ad,dp,cstatus]+helper[:3]+helper[5:])+"</td></tr>\n")

      if variant not in variantStats:
        helper = dict(map(lambda x:splitFields(x),VEP[variant][0][-1].split(";")))
        variantStats[variant] = [helper["SYMBOL"],VEP[variant][0][9],1]
      else:
        variantStats[variant][-1]+=1

      if float(gq) >= 30 and int(dp) >= 8:
        helper = dict(map(lambda x:splitFields(x),VEP[variant][0][-1].split(";")))
        if isBenign:
          ovariant_text = '<font color="#737373"><strong>'+helper["SYMBOL"]
        else:
          ovariant_text = "<strong>"+helper["SYMBOL"]
        
        if "EXON" in helper:
          ovariant_text +=" (E:%s)"%(helper["EXON"])
        elif "INTRON" in helper:
          ovariant_text +=" (I:%s)"%(helper["INTRON"])
              
        if ("HGVSc" in helper) and ("HGVSp" in helper):
          ovariant_text +=" %s / %s"%(":".join(helper["HGVSc"].split(":")[1:]),":".join(helper["HGVSp"].split(":")[1:]))
        elif "HGVSp" in helper:
          ovariant_text +=" %s"%(":".join(helper["HGVSp"].split(":")[1:]))
        elif "HGVSc" in helper:
          ovariant_text +=" %s"%(":".join(helper["HGVSc"].split(":")[1:]))
        if isBenign:
          ovariant_text += "</strong></font> [%s %s:%s %s]"%tuple(variant.split("_"))
        else:
          ovariant_text += "</strong> [%s %s:%s %s]"%tuple(variant.split("_"))
        o3_variants.append(ovariant_text)
      else:
        if status == None: 
          status = True
        status_flags.add("variants")
    outfile2.write("</table>\n")
    outfile2.write("""<p>GT - Genotype call encoded as allele values separated by "/". The allele values are 0 for the reference allele, 1 for the first alternative allele, 2 for the second allele.<br>
GQ - Conditional genotype quality, encoded as a phred quality -10*log_10(p) (genotype call is wrong, conditioned on the site's being variant); we only report GQ >= 30 on the summary page.<br>
AD - Read depth for each variant at this position (first reference, followed by alternative alleles).<br>
DP - Read depth at this position for this sample; we only report DP >= 8 on the summary page.</p>\n""")
  else:
    outfile2.write("<p>None</p>\n")
  
  outfile3.write("<td>"+("&nbsp;" if len(o3_variants)==0 else "<br>".join(o3_variants))+"</td>")

  o3_variants = []
  o3_variants_del = []
  tbl_variants_del = ''
  outfile2.write("""<p></p>\n""")
  outfile2.write("""<p><h2>Coverage by target region</h2></p>\n""")
  outfile2.write("""<table cellpadding="5" border="3">\n""")
  outfile2.write("<tr><th>"+"</th><th>".join(["Region","Length","Sites","Called","Fraction","Ave.Cov"])+"</th></tr>\n")
  lind, lind_del = None, None
  for ind,region in enumerate(sorted_regions):
    cvalue = coverage_stats_by_region[region][individual+":Cov"]
    gvalue = coverage_stats_by_region[region][individual+":GT"]
    tvalue = coverage_stats_by_region[region]['Count']
    chrom,start,end = name2region[region]
    length = end-start
    if tvalue > 0:
      outfile2.write("<tr><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%.2f%%</td><td>%.2f</td></tr>\n"%(region,length,tvalue,gvalue,gvalue/float(tvalue)*100,cvalue/float(tvalue)))
    else:
      outfile2.write("<tr><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%.2f%%</td><td>%.2f</td></tr>\n"%(region,length,tvalue,gvalue,0,0))
    
    if gvalue < tvalue:
      if lind == ind-1:
        start = o3_variants[-1].split(" - ")[0]
        o3_variants[-1] = start+" - "+region
      else:
        o3_variants.append(region)
      lind = ind
    if gvalue < tvalue*0.5:
      if lind_del == ind-1:
        start = o3_variants_del[-1].split(" - ")[0]
        o3_variants_del[-1] = start+" - "+region
      else:
        o3_variants_del.append(region)
      lind_del = ind
    
    if region not in coveragRegionStats:
      coveragRegionStats[region] = [tvalue,0,0,0]
    coveragRegionStats[region][1]+=cvalue
    coveragRegionStats[region][2]+=gvalue
    coveragRegionStats[region][3]+=1

  outfile2.write("</table>\n")
  outfile2.write("""<p></p>\n""")
  outfile2.write("""<p><h2>Regions with missing genotype calls</h2></p>\n""")
  covholes = []
  
  tbl1_fields.append('')
  if len(coverage_holes[individual]) > 0:
    outfile2.write("""<table cellpadding="5" border="3">\n""")
    outfile2.write("<tr><th>"+"</th><th>".join(["Chrom","Start","End","Region"])+"</th></tr>\n")
    for chrom,start,end in coverage_holes[individual]:
      regionstr = None
      if bedanno != None: 
        if chrom not in bedanno: continue
        for cinterval in bedanno[chrom].find(start,start+1):
          regionstr = cinterval.value
        for cinterval in bedanno[chrom].find(end,end+1):
          if regionstr == None or regionstr == cinterval.value:
            regionstr = cinterval.value
          else:
            regionstr = regionstr+" - "+cinterval.value
      if regionstr == None: regionstr = "NA"
      covholes.append("%s:%d-%d (%s)"%(chrom,start,end,regionstr))
      outfile2.write("<tr><td>%s</td><td>%d</td><td>%d</td><td>%s</td><tr>\n"%(chrom,start,end,regionstr))
      tbl1_fields[-1]+="%s(%s:%d-%d),"%(regionstr,chrom,start,end)
    outfile2.write("</table>\n")
  else:
    outfile2.write("""<p>None</p>\n""")
  tbl1_fields[-1]=tbl1_fields[-1].rstrip(',')
  tbl1_fields.append(",".join(map(lambda x: x.replace(" - ","-"),o3_variants_del)))
  
  if len(covholes) < 10:
    outfile3.write("<td>"+("&nbsp;" if len(covholes)==0 else "<br>".join(covholes))+"</td>")
  else:
    outfile3.write("<td>"+("&nbsp;" if len(o3_variants)==0 else "<br>".join(o3_variants))+"</td>")
  outfile3.write("<td>"+("&nbsp;" if len(o3_variants_del)==0 else "<br>".join(o3_variants_del))+"</td>")
  
  if individual in inversion_obs:
    tbl2_fields_inv1 = ['',individual,'','','failed']
    tbl2_fields_inv22 = ['',individual,'','','failed']
    
    outfile2.write("""<p></p>\n""")
    outfile2.write("""<p><h2>Inversion MIP results</h2></p>\n""")
    tobs = sum(inversion_obs[individual].values())
    outfile2.write("""<p>Total inversion MIP reads: <strong>%d</strong></p>\n"""%(tobs))
    minObsInv1,minObsInv22 = None,None
    if tobs > 0:
      outfile2.write("""<p>INV1 MIPs:</p>\n<table cellpadding="5" border="3">\n""")
      outfile2.write("""<tr><th>MIP name</th><th>Count</th></tr>\n""")
      tbl2helper = ""
      for invkey in ["inv1_1IU+1ID","inv1_1IU+1ED"]:
        outfile2.write("<tr><td>%s</td><td>%d</td><tr>\n"%(invkey,inversion_obs[individual][invkey]))
        if inversion_obs[individual][invkey] > 0:
          tbl2helper += "%s:%d,"%(invkey,inversion_obs[individual][invkey])
          if minObsInv1 == None: minObsInv1 = inversion_obs[individual][invkey]
          minObsInv1 = min(minObsInv1, inversion_obs[individual][invkey])
          
      outfile2.write("</table>\n")
      outfile2.write("<p></p>\n")
      tbl2helper=tbl2helper.strip(',')
      if tbl2helper == "":
        tbl2_fields_inv1[-1]='Failed'
      else:
        tbl2_fields_inv1[2]=tbl2helper
        if minObsInv1 < 8:
          tbl2_fields_inv1[-1]='Low'
        else:
          tbl2_fields_inv1[-1]='OK'
      
      outfile2.write("""<p>INV22 MIPs:</p>\n<table cellpadding="5" border="3">\n""")
      outfile2.write("""<tr><th>MIP name</th><th>Count</th></tr>\n""")
      tbl2helper = ""
      for invkey in ["inv22_ID+IU","inv22_ED+2U","inv22_ED+3U","inv22_ID+2U","inv22_ID+3U","inv22_ED+IU"]:
        outfile2.write("<tr><td>%s</td><td>%d</td><tr>\n"%(invkey,inversion_obs[individual][invkey]))
        if inversion_obs[individual][invkey] > 0:
          tbl2helper += "%s:%d,"%(invkey,inversion_obs[individual][invkey])
          if minObsInv22 == None: minObsInv22 = inversion_obs[individual][invkey]
          minObsInv22 = min(minObsInv22, inversion_obs[individual][invkey])
          
      outfile2.write("</table>\n")
      tbl2helper=tbl2helper.strip(',')
      if tbl2helper == "":
        tbl2_fields_inv22[-1]='Failed'
      else:
        tbl2_fields_inv22[2]=tbl2helper
        if minObsInv22 < 8:
          tbl2_fields_inv22[-1]='Low'
        else:
          tbl2_fields_inv22[-1]='OK'

    inv1found = set()
    #Check INT1
    for invname,flags in INT1_inversion_types:
      check = True
      for key,value in zip(inversion_names,flags):
        if (value and inversion_obs[individual][key] == 0) or (value == False and inversion_obs[individual][key] > 0):
          check = False
          break
      if check: 
        inv1found.add(invname)
        break

    #Check INT22
    found = False
    inv22found = set()
    for invname,flags in INT22_inversion_types+noINT22_inversion_types+INT22failed_inversion_types:
      check = True
      for key,value in zip(inversion_names,flags):
        if (value and inversion_obs[individual][key] == 0) or (value == False and inversion_obs[individual][key] > 0):
          check = False
          break
      if check: 
        inv22found.add(invname.split("#")[0])
        found = True
    if not found:
      inv22found.add("Conflict: INT22")
      #tbl2_fields_inv22[2]=''
      tbl2_fields_inv22[-1]='Conflict'

    outfile2.write("""<p><h3>Resulting inversion calls</h3></p>\n""")
    for invname in (inv1found | inv22found):
      invfields = invname.split("Conflict: ")
      conflict = len(invfields) > 1
      cinvname = invfields[-1]
      if conflict and status == None: 
        status = True
        status_flags.add("inversions")
      if "unknown" in cinvname and status == None: 
        status = True 
        status_flags.add("inversions")
      if "FAILED" in cinvname: 
        status = False
        status_flags.add("inversions")

      if cinvname not in inversionMipStats:
        inversionMipStats[cinvname] = [0,0]
      inversionMipStats[cinvname][0] += 1
      if conflict:
        inversionMipStats[cinvname][1] += 1
      outfile2.write("<p><em>%s</em></p>\n"""%(invname))
    
    tbl2_fields_inv1[-2]=",".join(inv1found)
    tbl2_fields_inv22[-2]=",".join(inv22found)

    tbl_out2.write('"%s"\r\n'%('","'.join(tbl2_fields_inv1)))
    tbl_out2.write('"%s"\r\n'%('","'.join(tbl2_fields_inv22)))

  outfile2.write("""<p></p>\n""")
  outfile2.write("""<p><h2>Under-/over-performing target-region MIPs</h2></p>\n""")
  tobs = 0 if individual not in failedMIPs else len(failedMIPs[individual])
  if tobs > 0:
    outfile2.write("""<table cellpadding="5" border="3">\n""")
    outfile2.write("<tr><th>MIP</th><th>+/-</th><th>Coordinates</th><th>Region</th><th>Count/Total</th><tr>\n")
    sind = sample2ind[individual]
    for mip in failedMIPs[individual]:
      direction = mip[0]
      mip = mip[2:]
      if mip in MIPcoords: 
        chrom,start,end = MIPcoords[mip]
        region_name = None
        if chrom in bedanno:
          for cinterval in bedanno[chrom].find(start,end):
            region_name = cinterval.value
          if region_name == None:
            for cinterval in bedanno[chrom].find(start-50,end+50):
              region_name = cinterval.value
        if region_name == None: region_name = ""
        gene_name = region_name.split("/")[0]
        outfile2.write("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%d/%d</td><tr>\n"%(mip,direction,"%s:%d-%d"%MIPcoords[mip],region_name,MIPcounts[mip][sind],TotalSample[sind]))
    outfile2.write("""</table>\n""")
  else:
    outfile2.write("""<p>None</p>\n""")
  outfile2.write("""<p>&nbsp;</p>\n""")
  outfile2.write("""<p>Go to <a href="report.html">Sample Summary Report</a></p>\n""")
  outfile2.write("""<p>Go to <a href="summary.html">Summary Report</a></p>\n""")
  outfile2.write("</body>\n</html>\n")
  outfile2.close()

  outfile3.write("<td>"+("&nbsp;" if len(inv1found)==0 else "<br>".join(map(lambda x: x if x.startswith("no") else "<strong>%s</strong>"%x,list(inv1found))))+"</td>")
  outfile3.write("<td>"+("&nbsp;" if len(inv22found)==0 else "<br>".join(map(lambda x: x if x.startswith("no") else "<strong>%s</strong>"%x,list(inv22found))))+"</td>")
  
  if status == None:
    outfile3.write("""<td><font color="#008000">OK</font></td>""")
    tbl1_fields.append('OK')
  elif status == True:
    outfile3.write("""<td><font color="#FF6600">CHECK<br>%s</font></td>"""%(",".join(status_flags)))
    tbl1_fields.append('CHECK: %s'%(",".join(status_flags)))
  else:
    outfile3.write("""<td><font color="#FF0000">FAILED<br>%s</font></td>"""%(",".join(status_flags)))
    tbl1_fields.append('FAILED: %s'%(",".join(status_flags)))
  tbl1_fields.append('')
  tbl_out1.write('"%s"\r\n'%('","'.join(tbl1_fields)))
  
  outfile3.write("</tr>\n")

  outfile.write("""<tr><td><a href="ind_%s.html">%s</a></td><td>%s</td><td>%d</td><td>%0.2f%%</td><td>%0.4f</td><td>%d</td><td>%d</td><td>%s</td></tr>\n"""%(
    individual,
    individual,
    sex_check[individual][0],
    GT_stats[individual],
    0 if total_sites == 0 else GT_stats[individual]/float(total_sites)*100,
    0 if total_sites == 0 else coverage_stats[individual]/float(total_sites),
    het_stats[individual],
    var_stats[individual],", ".join(map(lambda x: "%s:%s %s"%(tuple(x[0].split("_")[1:])),variants[individual]))))

outfile.write("</table>\n")
outfile.write("""<p></p>\n""")
outfile.write("""<p><h2>Variant summary</h2></p>\n""")
outfile.write("""<table cellpadding="5" border="3">\n""")
outfile.write("<tr><th>"+"</th><th>".join(["Variant","Gene","Effect","Count"])+"</th></tr>\n")
to_sort = map(lambda (var,(gene,effect,count)):(count,var,gene,effect),variantStats.iteritems())
to_sort.sort()
for count,var,gene,effect in to_sort[::-1]:
  varfields = var.split("_")
  var = "%s:%s %s"%(varfields[1],varfields[2],varfields[3])
  outfile.write('<tr><td>%s</td><td>%s</td><td>%s</td><td>%d</td></tr>\n'%(var,gene,effect,count))
outfile.write("</table>\n")

outfile.write("""<p></p>\n""")
outfile.write("""<p><h2>Per region coverage summary</h2></p>\n""")
outfile.write("""<table cellpadding="5" border="3">\n""")
outfile.write("<tr><th>"+"</th><th>".join(["Region","Length","Sites","Called","Fraction","Ave.Cov"])+"</th></tr>\n")
for region in sorted_regions:
  tvalue,cvalue,gvalue,count = coveragRegionStats[region]
  chrom,start,end = name2region[region]
  length = end-start
  if tvalue > 0 and count > 0:
    outfile.write("<tr><td>%s</td><td>%d</td><td>%d</td><td>%.2f</td><td>%.2f%%</td><td>%.2f</td></tr>\n"%(region,length,tvalue,gvalue/float(count),gvalue/float(tvalue*count)*100,cvalue/float(tvalue*count)))
  else:
    outfile.write("<tr><td>%s</td><td>%d</td><td>%d</td><td>%.2f</td><td>%.2f%%</td><td>%.2f</td></tr>\n"%(region,length,tvalue,0,0,0))
outfile.write("</table>\n")
  
outfile.write("""<p></p>\n""")
outfile.write("""<p><h2>Inversion MIP summary</h2></p>\n""")
outfile.write("""<table cellpadding="5" border="3">\n""")
outfile.write("<tr><th>"+"</th><th>".join(["Inversion type","Count","Conflict","%Conflict"])+"</th></tr>\n")
to_sort = inversionMipStats.keys()
to_sort.sort()
for invname in to_sort:
  values = inversionMipStats[invname]
  if values[1] > 0:
    outfile.write("<tr><td>%s</td><td>%d</td><td>%d</td><td>%.2f%%</td></tr>\n"%(invname,values[0],values[1],values[1]/float(values[0])*100))
  else:
    outfile.write("<tr><td>%s</td><td>%d</td></tr>\n"%(invname,values[0]))
outfile.write("</table>\n")

outfile.write("""<p></p>\n""")
outfile.write("""<p><h2>MIP performance summary</h2></p>\n""")
outfile.write("""<table cellpadding="5" border="3">\n""")
outfile.write("<tr><th>"+"</th><th>".join(["MIP","+/-","Coordinates","Region","Observed","Individuals"])+"</th></tr>\n")
to_sort = map(lambda (x,y):(len(y),x,y),failedMIPs_summary.iteritems())
to_sort.sort()
for count,mip,failedIndividuals in to_sort[::-1]:
  direction = mip[0]
  mip = mip[2:]
  #print mip,direction,MIPcoords
  if mip in MIPcoords: 
    chrom,start,end = MIPcoords[mip]
    region_name = None
    if chrom in bedanno:
      for cinterval in bedanno[chrom].find(start,end):
        region_name = cinterval.value
      if region_name == None:
        for cinterval in bedanno[chrom].find(start-50,end+50):
          region_name = cinterval.value
    if region_name == None: region_name = ""
    gene_name = region_name.split("/")[0]
    outfile.write("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%d</td><td>%s</td><tr>\n"%(mip,direction,"%s:%d-%d"%MIPcoords[mip],region_name,count,", ".join(failedIndividuals)))
outfile.write("</table>\n")

outfile.write("""<p>&nbsp;</p>""")
outfile.write("""<p>Go to <a href="report.html">Sample Summary Report</a></p>""")
outfile.write("""</body>\n</html>\n""")
outfile.close()

outfile3.write("""</table>\n""")
outfile3.write("""<p>&nbsp;</p>""")
outfile3.write("""<p>Go to <a href="summary.html">Summary Report</a></p>""")
outfile3.write("</body>\n</html>\n")
outfile3.close()

tbl_out1.close()
tbl_out2.close()
tbl_out3.close()
tbl_out3_benign.close()
