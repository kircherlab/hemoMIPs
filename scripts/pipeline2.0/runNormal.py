#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
Run the Illumina Pipeline

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *04.01.2012
:Type: MAIN
:Input: -
:Output: -

"""

import sys
import os
import time
import subprocess
import math

import xml.dom.minidom
from xml.dom.minidom import Node

from optparse import OptionParser
from optparse import OptionGroup

root_bin = '/mnt/solexa/bin/'

IbisInstallFolder = root_bin+"Ibis/"
RunBaseCalling = IbisInstallFolder+"runBaseCalling.py"
MergeReads = root_bin+"pipeline2.0/MergeTrimReadsBAM.py"
BAMFilter = root_bin+"pipeline2.0/FilterBAM.py"
DoubleIndexSplit = root_bin+"pipeline2.0/SplitFastQdoubleIndexBAM.py"
BwaBAM = root_bin+"pipeline2.0/BAM2BWAdirect.py"
FastQ2BAM = root_bin+"pipeline2.0/FastQ2BAM.py"
BAM2FastQ = root_bin+"pipeline2.0/BAM2FastQ.py"
BCL2FastQ = root_bin+"pipeline2.0/BCL2FastQ.py"
RTAreport = root_bin+"RTAreport/generate_report.py"
FastQCreport = root_bin+"FastQC/fastqc"
soap = root_bin+"soap_1.11/soap"

if not os.path.exists(RunBaseCalling):
  print "FIX",RunBaseCalling
  sys.exit()
if not os.path.exists(MergeReads):
  print "FIX",MergeReads
  sys.exit()
if not os.path.exists(BAMFilter):
  print "FIX",BAMFilter
  sys.exit()
if not os.path.exists(FastQ2BAM):
  print "FIX",FastQ2BAM
  sys.exit()
if not os.path.exists(BAM2FastQ):
  print "FIX",BAM2FastQ
  sys.exit()
if not os.path.exists(BCL2FastQ):
  print "FIX",BCL2FastQ
  sys.exit()
if not os.path.exists(DoubleIndexSplit):
  print "FIX",DoubleIndexSplit
  sys.exit()
if not os.path.exists(BwaBAM):
  print "FIX",BwaBAM
  sys.exit()
if not os.path.exists(RTAreport):
  print "FIX",RTAreport
  sys.exit()
if not os.path.exists(FastQCreport):
  print "FIX",FastQCreport
  sys.exit()
if not os.path.exists(soap):
  print "FIX",soap
  sys.exit()

#########################
### USER PARAMETER ...
#########################

parser = OptionParser(usage="usage: %prog [options]")
group = OptionGroup(parser, "General","General options")
group.add_option("-c", "--cores", dest="cores", help="Maximum number of CPU cores to be used (default 8)",default=8,type="int")
group.add_option("-r", "--run_folder", dest="runfolder", help="Path to RUN folder (Data/Intensities/BaseCalls subfolder with BCL files required)")
group.add_option("-o", "--outpath", dest="outpath", help="Path for output files (default /mnt/ngs_data/)",default="/mnt/ngs_data/")
group.add_option("--temp", dest="tmp", help="Path to temporary folder (default /mnt/solexa/tmp)",default="/mnt/solexa/tmp")
group.add_option("--lanes", help="Where possible consider only a subset of lanes (default '1-8')",default="1-8")
group.add_option("--skipBustard", dest="skipBustard", help="Skip any processing of Bustard results (default '1-8')",default="1-8")
group.add_option("--skipIbis", dest="skipIbis", help="Skip any processing of Ibis results (default '')",default='')
group.add_option("--skipFinal", dest="skipFinal", help="Skip merging/trimming and quality filters (default '')",default='')
group.add_option("--skipBWA", dest="skipBWA", help="Don't run BWA on last Ibis output file for lanes (default '')",default='')
group.add_option("--FastQ", dest="FastQ", help="Generate filtered FastQ files for downstream processing (default '')",default='')
group.add_option("--mock", dest="mock", help="Don't do anything, do a mock run",action="store_true",default=False)
parser.add_option_group(group)

group = OptionGroup(parser, "Run specific parameters")
group.add_option("--keys",dest="keys",help="Keys in first and second read (default '1-8:')",default="1-8:")
group.add_option("--readlengths",dest="readlengths",help="Length of each read (default '76')",default="76")
group.add_option("--indexfile",dest="indexfile",help="Per lane file(s) with index sequences used (first column seq, second column name; default: '1-8:')",default='1-8:')
group.add_option("--no_index_dist", dest="noindexdist", help="Do not consider distance for indexes for lanes (default '')",default='')
group.add_option("--index_quality", dest="index_quality", help="Quality score cutoff for index read (default '1-8:0')",default='1-8:0')
group.add_option("--indexreadlength", dest="indexreadlength", help="Index read length (if longer than defined by index sequences)",default=0,type='int')
group.add_option("--2nd_indexreadlength", dest="indexreadlength2", help="Index 2 read length (if longer than defined by index sequences)",default=0,type='int')
group.add_option("--adapter",dest="adapter",help="Adapter sequences seen in first and second read (default '1-8:')",default="1-8:")
group.add_option("--chimera",dest="chimera",help="Chimera sequence(s) observed for this experiment (default '1-8:')",default="1-8:")
parser.add_option_group(group)

group = OptionGroup(parser, "Ibis training parameters (passed to "+RunBaseCalling+")")
group.add_option("--training_index",dest="train_index",help="Index sequence matched by sequences used for training (default 'TTGCCGC')",default="TTGCCGC")
group.add_option("--training_lanes",dest="train_lanes",help="Lanes used for training the Ibis (default '' = all)",default="")
group.add_option("--training_tiles",dest="train_tiles",help="Tiles used for training the Ibis (default '' = all)",default="")
group.add_option("--training_coordtype", dest="train_coordtype", help="Type of cluster coordinate transformation: round = Pipeline <= 1.3: round(x), floor = Pipeline 1.5.*: floor(abs(x)), shift_round (default) = Pipeline 1.6.*: round(x*10+1000)",choices=["round","floor","shift_round"],default="shift_round")
group.add_option("--training_reference",dest="reference",help="Reference sequence used for training the Ibis (default /mnt/solexa/Genomes/phiX/control/whole_genome.fa)",default="/mnt/solexa/Genomes/phiX/control/whole_genome.fa")
group.add_option("--ConsiderKey",dest="considerKey",help="Consider key in Ibis training",default=False,action="store_true")
group.add_option("--skipOneCycleAfterKey",dest="skipOneCycleAfterKey",help="For sequences with key, skip 1st cycle after key when training",default=False,action="store_true")
group.add_option("--maskMM", dest="maskMM", help="Mask mismatches from training data",default=False,action="store_true")
parser.add_option_group(group)

group = OptionGroup(parser, "Trimming/merging parameters for (passed to "+MergeReads+")")
group.add_option("--oneErrorKey", dest="oneErrorKey", help="Allow base of key to be wrong/missing (default '', on: '1-8')",default='')
group.add_option("--trimCutoff", dest="adapterTrim", help="Lowest number of adapter bases to be observed before trimming (default '1-8:1')",default='1-8:1')
group.add_option("--mergeOverlap", dest="mergeOverlap", help="Also merge overlapping reads without adapter bases at the read ends (default '', on: '1-8')",default='')
parser.add_option_group(group)

group = OptionGroup(parser, "Filter parameter (passed to "+BAMFilter+")")
group.add_option("--trimLastN", dest="trimLastN", help="Remove last N bases of reads (default '1-8:0'  0 = off)",default='1-8:0')
group.add_option("--filter_qval", dest="QvalCutoff", help="Quality cutoff value (default = '1-8:15', -1 = off).",default='1-8:15')
group.add_option("--filter_number", dest="QvalNumber", help="Maximum number of bases below cutoff (default = '1-8:5', -1 for average).",default='1-8:5')
group.add_option("--qualityTrim", dest="qualityTrim", help="Trim reads back by quality. Remove reads shorter N (default '1-8:0', 0 = off)",default='1-8:0')
group.add_option("--complex_val", dest="CompCutoff", help="Complexity cutoff value (default = '1-8:-1', -1 = off).",default='1-8:-1')
group.add_option("--complex_method", dest="cMethod", help="Use fraction of most abundant base instead of entropy as measure (default '').",default='')
group.add_option("--length_val", dest="LengthCutoff", help="Length cutoff value (default = '1-8:0', 0 = off).",default='1-8:0')
group.add_option("--length_method", dest="lMethod", help="Keep sequences below cutoff instead above or equal to cutoff (default '').",default='')
parser.add_option_group(group)

group = OptionGroup(parser, "Run BWA on results (passed to "+BwaBAM+")")
group.add_option("--bwa_genomes", dest="bwa_genomes", help="Path to genomes (folder containing bwa-0.4.9* and whole_genome.fa/.fa.fai) (default '1-8:/mnt/solexa/Genomes/hg19_1000g/')",default="1-8:/mnt/solexa/Genomes/hg19_1000g/")
group.add_option("--bwa_params", dest="bwa_params", help="Parameters passed to bwa (default '1-8:')",default="1-8:")
parser.add_option_group(group)

(options, args) = parser.parse_args()


##############################
### RANGE OPTION PARSERS
##############################

def parse_rangestr(rangestr):
  res = []
  fields = rangestr.split(',')
  for elem in fields:
    if "-" in elem:
      se = elem.split('-')
      if len(se) == 2:
        try:
          start = int(se[0])
          end = int(se[1])
          res.extend(range(start,end+1))
        except: return None
      else: return None
    else:
      try:
        pos = int(elem)
        res.append(pos)
      except: return None
  if len(res) == 0: return None
  else:
    res = list(set(res))
    res.sort()
    return map(str,res)

def parse_range_set(rangestr):
  res = parse_rangestr(rangestr)
  if res == None:
    res = set()
  else:
    res = set(res)
  return res

def parse_range_value(value,lanes,name,ctype="string"):
  res = {}
  cases = value.split(';')
  for case in cases:
    fields = case.split(":")
    if len(fields) == 2:
      clanes = parse_rangestr(fields[0])
      value = fields[1]
      if ctype == "int":
        value = 0
        try: value = int(fields[1])
        except:
          sys.stderr.write('Error parsing int from %s parameter string.\n'%(name))
          sys.exit()
      elif ctype == "float":
        value = 0.0
        try: value = float(fields[1])
        except:
          sys.stderr.write('Error parsing float from %s parameter string.\n'%(name))
          sys.exit()
      for elem in clanes:
        res[elem] = value
  if len(set(res.keys()).intersection(lanes)) != len(lanes):
    sys.stderr.write('%s %s'%(set(res.keys()).intersection(lanes),lanes))
    sys.stderr.write('Incomplete definition of %s parameter.\n'%(name))
    sys.exit()
  return res


##############################
### HANDLE MULTIPLE JOBS
##############################

def total_free_memory():
  proc = subprocess.Popen('top -n 1 | grep "^[S|M]"',shell=True,stdout=subprocess.PIPE)
  stdout_list = (proc.communicate())[0].split()
  helper = []
  for elem in stdout_list:
    if elem[0].isdigit() and elem[-1]=="k":
      helper.append(int(elem[:-1]))
    elif elem[0].isdigit() and elem[-1]=="m":
      helper.append(int(elem[:-1])*1024)
    elif elem[0].isdigit() and elem[-1]=="g":
      helper.append(int(elem[:-1])*1024*1024)
  if len(helper) == 8:
    return helper[2]+helper[7]
  else:
    return None

jobs = []
dev_null = open('/dev/null','w')
def wait_jobs():
  global jobs
  while len(jobs) > 0:
    proc = jobs.pop(0)
    proc.wait()

def handle_jobs(cjob,silent=False):
  global options
  global jobs
  global dev_null
  if options.mock:
    print cjob
    return None
  else: print cjob
  if len(jobs) < options.cores:
    if not silent: jobs.append(subprocess.Popen(cjob,shell=True))
    else: jobs.append(subprocess.Popen(cjob,shell=True,stderr=dev_null))
  else:
    iteration = 0
    while len(jobs) >= options.cores:
      njobs = []
      for elem in jobs:
        if None == elem.poll():
          njobs.append(elem)
      jobs = njobs
      iteration += 1
      if len(jobs) >= options.cores:
        #DEAD LOCK?
        if iteration >= 200:
          wait_jobs()
        else: time.sleep(30)
    if not silent: jobs.append(subprocess.Popen(cjob,shell=True))
    else: jobs.append(subprocess.Popen(cjob,shell=True,stderr=dev_null))

##############################
### HELPER FUNCTIONS
##############################

def makedirs(path):
  global options
  if os.path.isdir(path): return False
  if not options.mock:
    try: os.makedirs(path)
    except: pass
    if not os.path.isdir(path):
      print "Error: Could not create",path
      return False
    else:
      return True
  else:
    print "os.makedirs(%s)"%path
    return True

def read_runinfo(filename):
  res = None
  if os.path.exists(filename):
    cycle = 1
    res = {}
    doc = xml.dom.minidom.parse(filename)
    for node in doc.getElementsByTagName("Run"):
      res['runid'] = node.getAttribute("Id")
      for node2 in node.getElementsByTagName("Instrument"):
        h = ""
        for node3 in node2.childNodes:
          if node3.nodeType == Node.TEXT_NODE:
            h += node3.data
        res['instrument'] = h
      for node2 in node.getElementsByTagName("Cycles"):
        res['cycles'] = int(node2.getAttribute("Incorporation"))
      for node2 in node.getElementsByTagName("Tiles"):
        res['lanes'] = 0
        res['tiles'] = 0
        for node3 in node2.getElementsByTagName("Lane"):
          res['tiles'] = max(res['tiles'],int(node3.getAttribute("Incorporation")))
          res['lanes'] += 1
      for node2 in node.getElementsByTagName("FlowcellLayout"):
        res['lanes'] = int(node2.getAttribute("LaneCount"))
        res['tiles'] = int(node2.getAttribute("TileCount"))
        try:
          factor = int(node2.getAttribute("SwathCount"))*int(node2.getAttribute("SurfaceCount"))
        except:
          factor = 1
        res['tiles'] = res['tiles']*factor
      res['read_ranges'] = []
      for node2 in node.getElementsByTagName("Reads"):
        for node3 in node2.getElementsByTagName("Read"):
          try: 
            res['read_ranges'].append((int(node3.getAttribute("FirstCycle")),int(node3.getAttribute("LastCycle"))))
          except:
            NumCycles = int(node3.getAttribute("NumCycles"))
            res['read_ranges'].append((cycle,cycle+NumCycles-1))
            cycle+=NumCycles
      res['fix']=False
      if "_PEDI" in res['runid'].upper() and len(res['read_ranges']) == 3:
        lindex = res['read_ranges'][1][1]-res['read_ranges'][1][0]+1
        res['read_ranges'][-1]=(res['read_ranges'][-1][0],res['read_ranges'][-1][1]-lindex)
        res['read_ranges'].append((res['read_ranges'][-1][1]+1,res['read_ranges'][-1][1]+lindex))
        res['fix'] = True
    if 'cycles' not in res and cycle > 0: res['cycles'] = cycle-1
    del doc
  return res


##############################
### EVALUALTE SCRIPT PARAMETERS
##############################

lanes = parse_rangestr(options.lanes)
if lanes == None:
  lanes = map(str,range(1,9))
if len(lanes) != 8:
  print "\nWill consider the following lanes:",lanes
lanes_set = set(lanes)


skipBustard = parse_range_set(options.skipBustard)
if len(lanes_set - skipBustard) == 0: options.skipBustard = True
else: options.skipBustard = False

skipIbis = parse_range_set(options.skipIbis)
if len(lanes_set - skipIbis) == 0: options.skipIbis = True
else: options.skipIbis = False

skipFinal = parse_range_set(options.skipFinal)
if len(lanes_set - skipFinal) == 0: options.skipFinal = True
else: options.skipFinal = False

skipBWA = parse_range_set(options.skipBWA)
if len(lanes_set - skipBWA) == 0: options.skipBWA = True
else: options.skipBWA = False

FastQ = parse_range_set(options.FastQ)
if len(lanes_set & FastQ) == 0: options.FastQ = False
else: options.FastQ = True

noIndexDist = parse_range_set(options.noindexdist)
indexQualCutoff = parse_range_value(options.index_quality,lanes_set,"index read quality cutoff","int")

oneErrorKey = parse_range_set(options.oneErrorKey)
adapterTrim = parse_range_value(options.adapterTrim,lanes_set,"adapter trim cutoff","int")
mergeOverlap = parse_range_set(options.mergeOverlap)

qualTrim = parse_range_value(options.qualityTrim,lanes_set,"length cutoff of quality filter","int")
trimLastN = parse_range_value(options.trimLastN,lanes_set,"remove last N bases cutofff","int")
qualCutoff = parse_range_value(options.QvalCutoff,lanes_set,"quality score cutoff","int")
qualNumberCutoff = parse_range_value(options.QvalNumber,lanes_set,"number of bases below quality score cutoff","int")
compCutoff = parse_range_value(options.CompCutoff,lanes_set,"complexity cutoff","float")
cMethod = parse_range_set(options.cMethod)
lengthCutoff = parse_range_value(options.LengthCutoff,lanes_set,"length cutoff","int")
lMethod = parse_range_set(options.lMethod)

bwa_params = parse_range_value(options.bwa_params,lanes_set-skipBWA,"BWA parameters")

options.train_index = options.train_index.strip().upper()

if not(os.path.isdir(options.outpath)):
  if not makedirs(options.outpath):
    sys.exit()
options.outpath=options.outpath.rstrip("/")+"/"

if not os.path.isdir(options.tmp):
  print "Error: Need valid path for temporary files"
  sys.exit()
options.tmp=options.tmp.rstrip("/")+"/"

if (options.runfolder == None) or (not os.path.isdir(options.runfolder)) or (not os.path.isdir(options.runfolder+"/Data/Intensities/BaseCalls/")):
  print "Error: Need valid path to run folder"
  sys.exit()
options.runfolder=options.runfolder.rstrip("/")+"/"

runinfo = None
if os.path.exists(options.runfolder+"/RunInfo.xml"):
  print "Evaluating RunInfo file..."
  runinfo = read_runinfo(options.runfolder+"/RunInfo.xml")
  if options.train_lanes.strip() == "" and 'lanes' in runinfo:
    if runinfo['lanes'] > 1:
      options.train_lanes = "1-%d"%(runinfo['lanes'])
    else:
      options.train_lanes = "1"
  if options.train_tiles.strip() == "" and 'tiles' in runinfo:
    if runinfo['tiles'] == 100: # GAII
      options.train_tiles = "1-100"
    elif runinfo['tiles'] == 120: # GAIIx
      options.train_tiles = "1-120"
    elif runinfo['tiles'] == 48: # HiSeq 3x8x2
      options.train_tiles = "1101-1108,1201-1208,1301-1308,2101-2108,2201-2208,2301-2308"
    elif runinfo['tiles'] == 36: # HiSeq 2x8x2
      options.train_tiles = "1101-1108,1201-1208,2101-2108,2201-2208"
    elif runinfo['tiles'] == 24: # HiScan 3x8
      options.train_tiles = "1101-1108,1201-1208,1301-1308"
    elif runinfo['tiles'] == 16: # HiScan 2x8
      options.train_tiles = "1101-1108,1201-1208"
    elif runinfo['tiles'] == 12: # MiSeq
      options.train_tiles = "1-12"

if not os.path.exists(options.runfolder+'NoIbis.txt') and (options.train_lanes.strip("") == "" or options.train_tiles.strip("") == ""):
  print "Error: Need valid Ibis training lanes and training tiles."
  sys.exit()

options.expID = None
options.folderID = None
if runinfo != None:
  options.expID = "_".join(runinfo['runid'].split("_")[1:])
  options.folderID = runinfo['runid']
  print "Experiment name obtained from RunInfo.xml:",options.expID
else:
  fields = options.runfolder.split("/")
  options.expID="_".join(fields[-1].split("_")[1:])
  options.folderID = fields[-1]
  print "Extracted experiment name from run folder name:",options.expID


# CHECK PROVIDED INDEX FILES...
print ""
lengthindex = {}
lengthindex2 = {}
indexnames = {}
isIndexRun = False
isDoubleIndexRun = False
hindexfiles = parse_range_value(options.indexfile,lanes_set,"index files")
for lane in lanes:
  lengthindex[lane]=0
  lengthindex2[lane]=0
  indexnames[lane] = set()
  indexnames[lane].add("conflict")
  indexnames[lane].add("unknown")
  indexnames[lane].add("wrong")
  if (lane in hindexfiles) and (hindexfiles[lane]!=""):
    if os.path.exists(hindexfiles[lane].replace("_index","_dindex")):
      print "Found double index file, replacing %s by %s"%(hindexfiles[lane],hindexfiles[lane].replace("_index","_dindex"))
      hindexfiles[lane] = hindexfiles[lane].replace("_index","_dindex")
    if os.path.exists(hindexfiles[lane]):
      infile = open(hindexfiles[lane])
      nucleotides=["A","C","G","T"]
      count_index = 0
      for line in infile:
        if (len(line.strip()) > 0) and (line[0] != "#"):
          fields = line.split()
          if len(fields) >= 3:
            cindex = fields[0].upper()
            # SKIP NON VALID LINES
            valid = True
            for elem in cindex:
              if elem not in nucleotides:
                valid=False
                break
            if not valid: continue
            if lengthindex[lane] < len(cindex): lengthindex[lane] = len(cindex)
            count_index += 1
            cindex = fields[1].upper()
            # CHECK VALID SECOND INDEX
            valid = True
            for elem in cindex:
              if elem not in nucleotides:
                valid=False
                break
            if not valid:
              indexnames[lane].add(fields[1].split("-")[-1])
            else:
              if lengthindex2[lane] < len(cindex): lengthindex2[lane] = len(cindex)
              indexnames[lane].add(fields[2].split("-")[-1])
          elif len(fields) >= 2:
            cindex = fields[0].upper()
            # SKIP NON VALID LINES
            valid = True
            for elem in cindex:
              if elem not in nucleotides:
                valid=False
                break
            if not valid: continue
            if lengthindex[lane] < len(cindex): lengthindex[lane] = len(cindex)
            count_index += 1
            indexnames[lane].add(fields[1].split("-")[-1])
      infile.close()
      if (count_index != 0) and (lengthindex[lane] == 0):
        print "Error: Given index file is not valid:",hindexfiles[lane]
        sys.exit()
      elif (count_index == 0):
        del hindexfiles[lane]
        isIndexRun = True
      else:
        print "%s defines %d index sequences with a length of %dnt for lane %s"%(hindexfiles[lane],count_index,lengthindex[lane],lane)
        isIndexRun = True
      if lengthindex2[lane] != 0: 
        if not isDoubleIndexRun: print "Run was identified as double index run."
        isIndexRun = True
        isDoubleIndexRun = True
    else:
      print "Error: Given index file is not available:",hindexfiles[lane]
      sys.exit()


readlength = None
ispaired = None
maxlengthindex,maxlengthindex2 = None,None
if runinfo != None:
  for ind,readrange in enumerate(runinfo['read_ranges']):
    creadlength = readrange[1]-readrange[0]+1
    if ind == 0:
      readlength = creadlength
    elif (creadlength < 10) and ind == 1:
      maxlengthindex = creadlength
    elif (creadlength >= 10) and ind == 1:
      readlength = (readlength,creadlength)
    elif (creadlength < 10) and ind > 1:
      maxlengthindex2 = creadlength
    elif (creadlength >= 10) and ind > 1:
      readlength = (readlength,creadlength)
if (maxlengthindex != None or maxlengthindex2 != None) and not isIndexRun:
  isIndexRun = True
  print "Run was identified as multiplex run from RunInfo.xml."
if maxlengthindex != None and maxlengthindex2 != None and not isDoubleIndexRun:
  isDoubleIndexRun = True
  print "Run was identified as double index run from RunInfo.xml."
if readlength != None and type(readlength) == type((1,1)): ispaired = True


# CHECK PROVIDED KEY SEQUENCES...
lane_keys = {}
options.keys = options.keys.strip()
if ((options.keys == "") or (options.keys.count(':') == 1 and options.keys.split(':')[1] == "")) and (ispaired == None):
  lane_keys = {'1':'','2':'','3':'','4':'','5':'','6':'','7':'','8':''}
  ispaired=False
elif ((options.keys == ",") or (options.keys.count(':') == 1 and options.keys.split(':')[1] == ",")) and (ispaired == None):
  lane_keys = {'1':('',''),'2':('',''),'3':('',''),'4':('',''),'5':('',''),'6':('',''),'7':('',''),'8':('','')}
  ispaired=True
elif ((options.keys == "") or (options.keys == ",") or (options.keys.count(':') == 1 and (options.keys.split(':')[1] == "," or options.keys.split(':')[1] == ""))):
  if ispaired: lane_keys = {'1':('',''),'2':('',''),'3':('',''),'4':('',''),'5':('',''),'6':('',''),'7':('',''),'8':('','')}
  else: lane_keys = {'1':'','2':'','3':'','4':'','5':'','6':'','7':'','8':''}
else:
  cases = options.keys.split(";")
  for case in cases:
    fields = case.split(":")
    if len(fields) == 2:
      clanes = parse_rangestr(fields[0])
      creads = fields[1].split(",")
      if len(creads) == 2:
        if (ispaired in [None,True]): ispaired = True
        else:
          print "Error: Cannot handle a mixed single read/paired end run:",clanes,creads
          sys.exit()
        for elem in clanes:
          lane_keys[elem]=(creads[0],creads[1])
      elif len(creads) == 1:
        if (ispaired in [None,False]): ispaired = False
        else:
          print "Error: Cannot handle a mixed single read/paired end run:",clanes,creads
          sys.exit()
        for elem in clanes:
          lane_keys[elem]=creads[0]
      else:
        print "Error: Unexpected pattern for keys",options.keys,"->",case,"->",creads
        sys.exit()
    else:
      print "Error: Unexpected pattern for keys",options.keys,"->",case
      sys.exit()
if (len(lane_keys) < len(lanes)):
  print "Expected",len(lanes),"lanes in keys argument got",len(lane_keys)
  sys.exit()

print ""
if ispaired: print "Parameters define a paired end run."
else: print "Parameters define a single read run."

print ""
print "Got the following keys from commandline options:"
print lane_keys

clen = None
clen2 = None
if options.considerKey:
  wrong = False
  for elem in parse_rangestr(options.train_lanes):
    if (elem in lane_keys):
      if ispaired:
        if (clen == None): clen=len(lane_keys[elem][0])
        elif (len(lane_keys[elem][0]) > clen):
          clen = len(lane_keys[elem][0])
          wrong = True
          break
        if (clen2 == None): clen2=len(lane_keys[elem][1])
        elif (len(lane_keys[elem][1]) > clen2):
          clen2 = len(lane_keys[elem][1])
          wrong = True
          break
      else:
        if (clen == None): clen=len(lane_keys[elem])
        elif (len(lane_keys[elem]) > clen):
          clen = len(lane_keys[elem])
          wrong = True
          break
  if wrong:
    print "Warning: Lanes used for training should have same key lengths. Considering maximum."
if clen == None: clen = 0
if clen2 == None: clen2 = 0

if readlength == None:
  fields = options.readlengths.split(",")
  if (len(fields) == 2) and ispaired:
    try:
      readlength = (int(fields[0]),int(fields[1]))
    except:
      print "Unexpected readlength parameter:",fields
      sys.exit()
  elif (len(fields) == 1) and not ispaired:
    try:
      readlength = int(fields[0])
    except:
      print "Unexpected readlength parameter:",fields
      sys.exit()
  else:
    print "Error: Read length parameter does not fit run type."
    sys.exit()
elif options.readlengths.strip() != "":
  try:
    fields = map(int,options.readlengths.split(","))
    if (len(fields) == 2 and tuple(fields) != readlength) or (len(fields) == 1 and fields[0] != readlength):
      print "Warning: Defined read length does not match determined read length! Continue with determined."
  except:
    print "Warning: Unexpected readlength parameter: %s. Continue with determined read length values."%(options.readlengths)

if maxlengthindex == None:
  maxlengthindex=max(options.indexreadlength,max(lengthindex.values()))
  if maxlengthindex > 0: isIndexRun = True
if maxlengthindex2 == None:
  maxlengthindex2=max(options.indexreadlength2,max(lengthindex2.values()))
  if maxlengthindex2 > 0 and ispaired: isDoubleIndexRun = True
if not isDoubleIndexRun: lengthindex2 = {}

print ""
print "Determined the following readlength parameters:"

print readlength,"Index:",maxlengthindex,maxlengthindex2

lane_adapter = {}
options.adapter = options.adapter.strip()
if (options.adapter == "") or (options.adapter == ",") or (options.adapter.count(':') == 1 and (options.adapter.split(':')[1] == "" or options.adapter.split(':')[1] == ",")):
  if ispaired: lane_adapter = {'1':('',''),'2':('',''),'3':('',''),'4':('',''),'5':('',''),'6':('',''),'7':('',''),'8':('','')}
  else: lane_adapter = {'1':'','2':'','3':'','4':'','5':'','6':'','7':'','8':''}
else:
  cases = options.adapter.split(";")
  for case in cases:
    fields = case.split(":")
    if len(fields) == 2:
      clanes = parse_rangestr(fields[0])
      creads = fields[1].split(",")
      if len(creads) == 2 and ispaired:
        for elem in clanes:
          lane_adapter[elem]=(creads[0],creads[1])
      elif len(creads) == 1 and not ispaired:
        for elem in clanes:
          lane_adapter[elem]=creads[0]
      else:
        print "Unexpected pattern for adapters",options.adapter,"->",case,"->",creads
        sys.exit()
    else:
      print "Unexpected pattern for adapters",options.adapter,"->",case
      sys.exit()
if (len(lane_adapter) == 0) and ispaired:
  lane_keys = {'1':('',''),'2':('',''),'3':('',''),'4':('',''),'5':('',''),'6':('',''),'7':('',''),'8':('','')}
if (len(lane_adapter) == 0) and not ispaired:
  lane_keys = {'1':'','2':'','3':'','4':'','5':'','6':'','7':'','8':''}
elif (len(lane_adapter) < len(lanes)):
  print "Expected",len(lanes),"lanes in adapter argument got",len(lane_adapter)
  sys.exit()

print ""
print "Got the following adapter sequences from commandline options:"
print lane_adapter

lane_chimera = {}
options.chimera = options.chimera.strip()
if (options.chimera == "") or (options.chimera.split(':')[1] == ""):
  lane_chimera = {'1':[],'2':[],'3':[],'4':[],'5':[],'6':[],'7':[],'8':[]}
else:
  cases = options.chimera.split(";")
  for case in cases:
    fields = case.split(":")
    if len(fields) == 2:
      clanes = parse_rangestr(fields[0])
      creads = fields[1].split(",")
      for elem in clanes:
        lane_chimera[elem]=creads
    else:
      print "Unexpected pattern for chimeras",options.chimera,"->",case
      sys.exit()
if (len(lane_chimera) == 0):
  lane_chimera = {'1':[],'2':[],'3':[],'4':[],'5':[],'6':[],'7':[],'8':[]}
elif (len(lane_chimera) < len(lanes)):
  print "Expected",len(lanes),"lanes in chimera argument got",len(lane_chimera)
  sys.exit()

print ""
print "Got the following chimera sequences from commandline options:"
print lane_chimera

bwa_genomes = {}
if not options.skipBWA:
  bwa_genomes = parse_range_value(options.bwa_genomes,lanes_set-skipBWA,"BWA genome")
  for lane,path in bwa_genomes.iteritems():
    if not os.path.exists(path+"/bwa-0.4.9.bwt") or not os.path.exists(path+"/whole_genome.fa") and os.path.exists(path+"/whole_genome.fa.fai"):
      print "Path to bwa genome is not valid",path
      sys.exit()
print ""
print "Got the following bwa genome parameters from default/commandline options:"
print bwa_genomes

firecrest = options.runfolder+'Data/Intensities/'
bustard = firecrest+'BaseCalls/'
ibisfolder = firecrest+"Ibis_Basecall_"+time.strftime("%d-%m-%Y", time.localtime())+"/"
for elem in os.listdir(firecrest):
  if (elem.startswith("SVM_Basecall_") and os.path.isdir(firecrest+elem)) or (elem.startswith("Ibis_Basecall_") and os.path.isdir(firecrest+elem)):
    ibisfolder = firecrest+elem+"/"
    print "Found SVM/Ibis folder:",ibisfolder
    break

##############################
### DOING ACTUAL PROCESSING
##############################

step = 1
print ""
print "Seems we have everything to start..."

print str(step)+". Creating output folders."
makedirs(options.outpath+options.folderID+"/Bustard/Report/")
if not os.path.exists(options.runfolder+'NoIbis.txt'): makedirs(options.outpath+options.folderID+"/Ibis/FastQC/")
if (not options.skipBustard):
  makedirs(options.outpath+options.folderID+"/Bustard/Raw_Sequences/")
  if (not options.skipFinal): makedirs(options.outpath+options.folderID+"/Bustard/Final_Sequences/")
  if (not options.skipBWA): makedirs(options.outpath+options.folderID+"/Bustard/BWA/")
  if (options.FastQ): makedirs(options.outpath+options.folderID+"/Bustard/FastQ/")
if (not options.skipIbis):
  makedirs(options.outpath+options.folderID+"/Ibis/Raw_Sequences/")
  if (not options.skipFinal): makedirs(options.outpath+options.folderID+"/Ibis/Final_Sequences/")
  if (not options.skipBWA): makedirs(options.outpath+options.folderID+"/Ibis/BWA/")
  if (options.FastQ): makedirs(options.outpath+options.folderID+"/Ibis/FastQ/")

step+=1
print ""
print str(step)+". Save parameter file",options.outpath+options.folderID+"/pipeline_params"
outfile = sys.stdout
if (not options.mock):
  if not os.path.exists(options.outpath+options.folderID+"/pipeline_params"):
    outfile = open(options.outpath+options.folderID+"/pipeline_params",'w')
  else:
    outfile = open(options.outpath+options.folderID+"/pipeline_params",'a')
    outfile.write('\n')
outfile.write("#### Pipeline 2.0 run: "+time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())+" ####\n")
outfile.write(" ".join(map(lambda x: x if ("=" not in x) or (len(x.split('=')) == 2 and x.split('=')[1].replace('.','').isdigit()) else x.replace("=","='")+"'",sys.argv))+"\n\n")
for elem in dir(options):
  if (elem not in ['__cmp__', '__doc__', '__init__', '__module__', '__repr__', '__str__', '_update', '_update_careful', '_update_loose','ensure_value','read_file', 'read_module']):
    outfile.write(elem+"\t"+str(eval("options."+elem))+"\n")
outfile.write("GLOBAL:soap\t"+str(soap)+"\n")
outfile.write("GLOBAL:IbisInstallFolder\t"+str(IbisInstallFolder)+"\n")
outfile.write("GLOBAL:RunBaseCalling\t"+str(RunBaseCalling)+"\n")
outfile.write("GLOBAL:MergeReads\t"+str(MergeReads)+"\n")
outfile.write("GLOBAL:BAMFilter\t"+str(BAMFilter)+"\n")
outfile.write("GLOBAL:DoubleIndexSplit\t"+str(DoubleIndexSplit)+"\n")
outfile.write("GLOBAL:BwaBAM\t"+str(BwaBAM)+"\n")
outfile.write("GLOBAL:FastQ2BAM\t"+str(FastQ2BAM)+"\n")
outfile.write("GLOBAL:BAM2FastQ\t"+str(BAM2FastQ)+"\n")
outfile.write("GLOBAL:BCL2FastQ\t"+str(BCL2FastQ)+"\n")
outfile.write("PARSED:skipIbis\t"+str(skipIbis)+"\n")
outfile.write("PARSED:skipBustard\t"+str(skipBustard)+"\n")
outfile.write("PARSED:skipFinal\t"+str(skipFinal)+"\n")
outfile.write("PARSED:skipBWA\t"+str(skipBWA)+"\n")
outfile.write("PARSED:FastQ\t"+str(FastQ)+"\n")
if (not options.mock): outfile.close()


step+=1
print ""
print str(step)+". Copy sample sheet file..."
import re
samplesheet = re.compile("^([F]{0,1}[C]{0,1}[0-9]{2,3}[a-zA-Z0-9]{5,7}|[A-Z0-9]{9,10}).xml$")
found = False
for elem in os.listdir(options.runfolder):
  if samplesheet.match(elem) != None:
    handle_jobs("cp "+options.runfolder+elem+" "+options.outpath+options.folderID+"/")
    found = True
if not found and os.path.exists(options.runfolder+"SampleSheet.csv"):
  handle_jobs("cp "+options.runfolder+"SampleSheet.csv "+options.outpath+options.folderID+"/")


if options.mock or (len(os.listdir(options.outpath+options.folderID+"/Bustard/Report/")) == 0):
  step+=1
  print ""
  print str(step)+". Generating/Copying RTA or GERALD Report..."
  if os.path.isdir(bustard+"/RTAreport/"):
    print "Copying available RTA report..."
    handle_jobs("cp -R "+bustard+"/RTAreport/* "+options.outpath+options.folderID+"/Bustard/Report/")
  else:
    geraldfolder = filter(lambda x:x.startswith("GERALD") and os.path.isdir(bustard+"/"+x),os.listdir(bustard))
    if len(geraldfolder) >= 1: 
      geraldfolder = geraldfolder[0]
      print "Copying available GERALD report..."
      if os.path.isdir(geraldfolder+"/Plots"):
        handle_jobs("cp -R "+geraldfolder+"/Plots "+options.outpath+options.folderID+"/Bustard/Report/")
        handle_jobs("cp "+geraldfolder+"*.xml "+geraldfolder+"*.xsl "+geraldfolder+"*.htm "+geraldfolder+"config.txt "+options.outpath+options.folderID+"/Bustard/Report/")
        handle_jobs("cp -R "+bustard+"/Plots "+options.outpath+options.folderID+"/Bustard/")
        handle_jobs("cp -R "+bustard+"/*.htm "+bustard+"/BustardSummary.x* "+options.outpath+options.folderID+"/Bustard/")
      else:
        handle_jobs("cp "+geraldfolder+"*.png "+geraldfolder+"*.htm "+options.outpath+options.folderID+"/Bustard/Report/")
    else:
      print "Generating RTA report and copying it..."
      handle_jobs(RTAreport+" "+options.runfolder)
      wait_jobs()
      handle_jobs("cp -R "+bustard+"/RTAreport/* "+options.outpath+options.folderID+"/Bustard/Report/")


if not os.path.exists(options.runfolder+'NoIbis.txt') and ((not os.path.isdir(ibisfolder) or (len(os.listdir(ibisfolder)) < 4))):
  step+=1
  print ""
  print str(step)+". Running Ibis training and prediction (handing over to "+RunBaseCalling+")..."
  wait_jobs() # RunBaseCalling will all cores...
  if options.skipOneCycleAfterKey and ((clen > 0) or (clen2 > 0)):
    if ispaired:
      if (clen > 0) and (clen2 == 0): start=str(clen+2)+","+str(readlength[0]+maxlengthindex+1+clen2)
      elif (clen == 0) and (clen2 > 0): start=str(clen+1)+","+str(readlength[0]+maxlengthindex+2+clen2)
      else: start=str(clen+2)+","+str(readlength[0]+maxlengthindex+2+clen2)
      end=str(readlength[0])+","+str(readlength[0]+maxlengthindex+readlength[1])
    else:
      if clen > 0: start=str(clen+2)
      else: start=str(clen+1)
      end=str(readlength)
  else:
    if ispaired:
      start=str(clen+1)+","+str(readlength[0]+maxlengthindex+1+clen2)
      end=str(readlength[0])+","+str(readlength[0]+maxlengthindex+readlength[1])
    else:
      start=str(clen+1)
      end=str(readlength)
  params = ""
  if options.maskMM:
    params += " --maskMM"
  if isIndexRun and options.train_index != '':
    params += " --control_index="+options.train_index
  handle_jobs(RunBaseCalling+" --coordianteType="+options.train_coordtype+" --NoFinishCheck -c "+str(options.cores)+" -e "+options.expID+" -o "+ibisfolder+" -b "+bustard+" --temp="+options.tmp+" --start="+start+" --end="+end+" --indexlength="+str(maxlengthindex)+" --2nd_indexlength="+str(maxlengthindex2)+" -l "+options.train_lanes+" -t "+options.train_tiles+" -a "+soap+" -r "+options.reference+params)
  wait_jobs()  # RunBaseCalling will take all cores...
if not os.path.exists(options.runfolder+'NoIbis.txt') and not os.path.exists(options.outpath+options.folderID+"/Ibis/error_profile.pdf"):
  handle_jobs("cd "+ibisfolder+"Models; R --vanilla --quiet < "+IbisInstallFolder+"plot_error.R; cp error_profile.pdf "+options.outpath+options.folderID+"/Ibis/")


if (options.mock and not os.path.exists(options.runfolder+'NoIbis.txt')) or os.path.isdir(ibisfolder):
  step+=1
  print ""
  print str(step)+". Generate FastQC reports from Ibis raw sequences..."
  for lane in lanes_set:
    if not os.path.exists(options.outpath+options.folderID+"/Ibis/FastQC/s_"+lane+"_sequence_fastqc.zip") and (os.path.exists("%s/s_%s_sequence.txt.gz"%(ibisfolder,lane)) or options.mock):
      conversion_str = FastQCreport+" -q -o %s/Ibis/FastQC/ -f fastq %s/s_%s_sequence.txt.gz"%(options.outpath+options.folderID,ibisfolder,lane)
      handle_jobs(conversion_str)


if (not options.skipBustard):
  step+=1
  print ""
  print str(step)+". Extracting Bustard raw sequences and converting to raw sequence BAMs..."
  for lane in (lanes_set - skipBustard):
    if not os.path.exists(options.outpath+options.folderID+"/Bustard/Raw_Sequences/s_"+lane+"_sequence.bam"):
      if ispaired:
        max_cycle = readlength[0]+maxlengthindex+readlength[1]+maxlengthindex2
      else:
        max_cycle = readlength+maxlengthindex+maxlengthindex2
      conversion_str = BCL2FastQ+" -e %s -p %s -l %s --max_cycle=%d --PIPE | "%(options.expID,bustard,lane,max_cycle)
      if isIndexRun:
        if isDoubleIndexRun:
          conversion_str += DoubleIndexSplit+" --summary -l %d -m %d --PIPE "%(maxlengthindex,maxlengthindex2)
        else:
          conversion_str += DoubleIndexSplit+" --summary -l %d -m 0 --PIPE "%(maxlengthindex)
        if (lengthindex[lane] > 0):
          conversion_str += "-i %s "%hindexfiles[lane]
          if lane in options.noindexdist:
            conversion_str += "--no_skip_first_base --no_mutants --no_Ns "
        if ispaired:
          conversion_str += " --start=%d "%(readlength[0]+maxlengthindex+1)
        if lane in indexQualCutoff:
          conversion_str += " --quality=%i "%indexQualCutoff[lane]
      else:
        conversion_str += FastQ2BAM+" --PIPE "
        if ispaired:
          conversion_str += "--start=%d "%(readlength[0]+maxlengthindex+1)
      conversion_str += "> %s/Bustard/Raw_Sequences/s_%s_sequence.bam"%(options.outpath+options.folderID,lane)
      handle_jobs(conversion_str)


if (not options.skipIbis):
  wait_jobs()
  step+=1
  print ""
  print str(step)+". Converting compressed Ibis sequences to raw sequence BAMs."
  for lane in (lanes_set - skipIbis):
    if not os.path.exists(options.outpath+options.folderID+"/Ibis/Raw_Sequences/s_"+lane+"_sequence.bam") and (os.path.exists("%s/s_%s_sequence.txt.gz"%(ibisfolder,lane)) or options.mock):
      conversion_str = "zcat %s/s_%s_sequence.txt.gz | "%(ibisfolder,lane)
      if isIndexRun:
        if isDoubleIndexRun:
          conversion_str += DoubleIndexSplit+" --summary -l %d -m %d --PIPE "%(maxlengthindex,maxlengthindex2)
        else:
          conversion_str += DoubleIndexSplit+" --summary -l %d -m 0 --PIPE "%(maxlengthindex)
        if (lengthindex[lane] > 0):
          conversion_str += "-i %s "%hindexfiles[lane]
          if lane in options.noindexdist:
            conversion_str += "--no_skip_first_base --no_mutants --no_Ns "
        if ispaired:
          conversion_str += "--start=%d "%(readlength[0]+maxlengthindex+1)
        if lane in indexQualCutoff:
          conversion_str += "--quality=%i "%indexQualCutoff[lane]
      else:
        conversion_str += FastQ2BAM+" --PIPE "
        if ispaired:
          conversion_str += "--start=%d "%(readlength[0]+maxlengthindex+1)
      conversion_str += "> %s/Ibis/Raw_Sequences/s_%s_sequence.bam"%(options.outpath+options.folderID,lane)
      handle_jobs(conversion_str)


if (not options.skipBustard) and (not options.skipFinal):
  wait_jobs()
  step+=1
  print ""
  print str(step)+". Adapter trim/merge and filter sequences from Bustard."
  for lane in (lanes_set - skipBustard - skipFinal):
    if not os.path.exists(options.outpath+options.folderID+"/Bustard/Final_Sequences/s_"+lane+"_sequence.bam") and (os.path.exists(options.outpath+options.folderID+"/Bustard/Raw_Sequences/s_"+lane+"_sequence.bam") or options.mock):
      conversion_str = "cat "+options.outpath+options.folderID+"/Bustard/Raw_Sequences/s_"+lane+"_sequence.bam | "+MergeReads+" --PIPE "
      if ispaired:
        conversion_str += "-k '%s,%s' -f '%s' -s '%s' -c '%s' "%(lane_keys[lane][0],lane_keys[lane][1],lane_adapter[lane][0],lane_adapter[lane][1],",".join(lane_chimera[lane]))
        if lane in mergeOverlap: conversion_str += "--mergeoverlap "
      else:
        conversion_str += "-k '%s' -f '%s' -c '%s' "%(lane_keys[lane],lane_adapter[lane],",".join(lane_chimera[lane]))
      if lane in oneErrorKey: conversion_str += "--allowMissing "
      conversion_str += "-t %d "%(adapterTrim[lane])
      if not (qualTrim[lane] == 0 and trimLastN[lane] == 0 and qualCutoff[lane] == -1 and compCutoff[lane] == -1 and lengthCutoff[lane] == 0):
        conversion_str += "| "+BAMFilter+" --PIPE "
        if compCutoff[lane] != -1:
          if lane in cMethod: conversion_str += "--frequency --comp_cutoff=%.4f "%(compCutoff[lane])
          else: conversion_str += "--entropy --comp_cutoff=%.4f "%(compCutoff[lane])
        if lengthCutoff[lane] != 0:
          if lane in qualTrim and qualTrim[lane] != 0:
            if lane in lMethod: conversion_str += "-l %d -m %d "%(qualTrim[lane],lengthCutoff[lane])
            else: conversion_str += "-l %d -m -1 "%(max(lengthCutoff[lane],qualTrim[lane]))
          else:
            if lane in lMethod: conversion_str += "-l 0 -m %d "%(lengthCutoff[lane])
            else: conversion_str += "-l %d -m -1 "%(lengthCutoff[lane])
        if qualCutoff[lane] != -1:
          if qualNumberCutoff[lane] == -1: conversion_str += "--average --qual_cutoff=%d "%(qualCutoff[lane])
          elif lane in qualTrim and qualTrim[lane] != 0: conversion_str += "--trim --qual_cutoff=%d --qual_number=%d "%(qualCutoff[lane],qualNumberCutoff[lane])
          else: conversion_str += "--quality --qual_cutoff=%d --qual_number=%d "%(qualCutoff[lane],qualNumberCutoff[lane])
        if trimLastN[lane] != 0:
          conversion_str += "--clip=-%d "%abs(trimLastN[lane])
      conversion_str += "> "+options.outpath+options.folderID+"/Bustard/Final_Sequences/s_"+lane+"_sequence.bam"
      handle_jobs(conversion_str)


if (not options.skipIbis) and (not options.skipFinal):
  wait_jobs()
  step+=1
  print ""
  print str(step)+". Adapter trim/merge and filter sequences from Ibis."
  for lane in (lanes_set - skipIbis - skipFinal):
    if not os.path.exists(options.outpath+options.folderID+"/Ibis/Final_Sequences/s_"+lane+"_sequence.bam") and (os.path.exists(options.outpath+options.folderID+"/Ibis/Raw_Sequences/s_"+lane+"_sequence.bam") or options.mock):
      conversion_str = "cat "+options.outpath+options.folderID+"/Ibis/Raw_Sequences/s_"+lane+"_sequence.bam | "+MergeReads+" --PIPE "
      if ispaired:
        conversion_str += "-k '%s,%s' -f '%s' -s '%s' -c '%s' "%(lane_keys[lane][0],lane_keys[lane][1],lane_adapter[lane][0],lane_adapter[lane][1],",".join(lane_chimera[lane]))
        if lane in mergeOverlap: conversion_str += "--mergeoverlap "
      else:
        conversion_str += "-k '%s' -f '%s' -c '%s' "%(lane_keys[lane],lane_adapter[lane],",".join(lane_chimera[lane]))
      if lane in oneErrorKey: conversion_str += "--allowMissing "
      conversion_str += "-t %d "%(adapterTrim[lane])
      if not (qualTrim[lane] == 0 and trimLastN[lane] == 0 and qualCutoff[lane] == -1 and compCutoff[lane] == -1 and lengthCutoff[lane] == 0):
        conversion_str += "| "+BAMFilter+" --PIPE "
        if compCutoff[lane] != -1:
          if lane in cMethod: conversion_str += "--frequency --comp_cutoff=%.4f "%(compCutoff[lane])
          else: conversion_str += "--entropy --comp_cutoff=%.4f "%(compCutoff[lane])
        if lengthCutoff[lane] != 0:
          if lane in qualTrim and qualTrim[lane] != 0:
            if lane in lMethod: conversion_str += "-l %d -m %d "%(qualTrim[lane],lengthCutoff[lane])
            else: conversion_str += "-l %d -m -1 "%(max(lengthCutoff[lane],qualTrim[lane]))
          else:
            if lane in lMethod: conversion_str += "-l 0 -m %d "%(lengthCutoff[lane])
            else: conversion_str += "-l %d -m -1 "%(lengthCutoff[lane])
        if qualCutoff[lane] != -1:
          if qualNumberCutoff[lane] == -1: conversion_str += "--average --qual_cutoff=%d "%(qualCutoff[lane])
          elif lane in qualTrim and qualTrim[lane] != 0: conversion_str += "--trim --qual_cutoff=%d --qual_number=%d "%(qualCutoff[lane],qualNumberCutoff[lane])
          else: conversion_str += "--quality --qual_cutoff=%d --qual_number=%d "%(qualCutoff[lane],qualNumberCutoff[lane])
        if trimLastN[lane] != 0:
          conversion_str += "--clip=-%d "%abs(trimLastN[lane])
      conversion_str += "> "+options.outpath+options.folderID+"/Ibis/Final_Sequences/s_"+lane+"_sequence.bam"
      handle_jobs(conversion_str)


tmemfree = total_free_memory()
mem_instance = 3.5*1024*1024
if ispaired:
  mem_instance = mem_instance*3
try: 
  max_instance = min(math.floor(tmemfree/mem_instance),len(lanes_set - skipBWA))
  max_threads = int(math.floor(options.cores/float(max_instance)))
  if max_threads == 0: max_threads = 1
except:
  max_instance = 1
  max_threads = 1
print ""
print "Max BWA instances: %d, Number of cores per instance: %d"%(max_instance,max_threads)

if (not options.skipBustard) and (not options.skipBWA):
  wait_jobs()
  step+=1
  count_instance = 0
  print ""
  print str(step)+". Bwa mappings of Bustard reads..."
  for lane in (lanes_set - skipBustard - skipBWA):
    if lane not in bwa_genomes: continue
    infilename = options.outpath+options.folderID+"/Bustard/Final_Sequences/s_"+lane+"_sequence.bam"
    if lane in skipFinal: infilename = options.outpath+options.folderID+"/Bustard/Raw_Sequences/s_"+lane+"_sequence.bam"
    if not os.path.exists(infilename) and not options.mock:
      print "Error: No input file for BWA alignment!"
      continue
    if os.path.exists(options.outpath+options.folderID+"/Bustard/BWA/s_"+lane+"_sequence"+("_"+bwa_params[lane].replace(" ","_")+"_"+bwa_genomes[lane].strip('/').split('/')[-1]).replace("__","_")+".bam"): continue
    conversion_str = BwaBAM+" --cores=%.0f --outprefix=%s -o %s --bwa_genome='%s' --bwa_params='%s' %s"%(max_threads,"s_"+lane+"_sequence",options.outpath+options.folderID+"/Bustard/BWA/",bwa_genomes[lane],bwa_params[lane],infilename)
    if count_instance >= max_instance:
      wait_jobs()
      count_instance = 0
    handle_jobs(conversion_str)
    count_instance += 1
  wait_jobs()
  for lane in (lanes_set - skipIbis - skipBWA):
    if lane not in bwa_genomes: continue
    infilename = options.outpath+options.folderID+"/Ibis/Final_Sequences/s_"+lane+"_sequence.bam"
    if lane in skipFinal: infilename = options.outpath+options.folderID+"/Ibis/Raw_Sequences/s_"+lane+"_sequence.bam"
    if not os.path.exists(infilename) and not options.mock: continue
    outfilename = options.outpath+options.folderID+"/Ibis/BWA/s_"+lane+"_sequence"+("_"+bwa_params[lane].replace(" ","_")+"_"+bwa_genomes[lane].strip('/').split('/')[-1]).replace("__","_")+".bam"
    if os.path.exists(outfilename):
      if not options.mock: os.remove(infilename)
      handle_jobs('ln -s %s %s'%(outfilename,infilename))


if (not options.skipIbis) and (not options.skipBWA):
  wait_jobs()
  step+=1
  count_instance = 0
  print ""
  print str(step)+". Bwa mappings of Ibis reads..."
  for lane in (lanes_set - skipIbis - skipBWA):
    if lane not in bwa_genomes: continue
    infilename = options.outpath+options.folderID+"/Ibis/Final_Sequences/s_"+lane+"_sequence.bam"
    if lane in skipFinal: infilename = options.outpath+options.folderID+"/Ibis/Raw_Sequences/s_"+lane+"_sequence.bam"
    if not os.path.exists(infilename) and not options.mock:
      print "Error: No input file for BWA alignment!"
      continue
    if os.path.exists(options.outpath+options.folderID+"/Ibis/BWA/s_"+lane+"_sequence"+("_"+bwa_params[lane].replace(" ","_")+"_"+bwa_genomes[lane].strip('/').split('/')[-1]).replace("__","_")+".bam"): continue
    conversion_str = BwaBAM+" --cores=%.0f --outprefix=%s -o %s --bwa_genome='%s' --bwa_params='%s' %s"%(max_threads,"s_"+lane+"_sequence",options.outpath+options.folderID+"/Ibis/BWA/",bwa_genomes[lane],bwa_params[lane],infilename)
    if count_instance >= max_instance:
      wait_jobs()
      count_instance = 0
    handle_jobs(conversion_str)
    count_instance += 1
  wait_jobs()
  for lane in (lanes_set - skipIbis - skipBWA):
    if lane not in bwa_genomes: continue
    infilename = options.outpath+options.folderID+"/Ibis/Final_Sequences/s_"+lane+"_sequence.bam"
    if lane in skipFinal: infilename = options.outpath+options.folderID+"/Ibis/Raw_Sequences/s_"+lane+"_sequence.bam"
    if not os.path.exists(infilename) and not options.mock: continue
    outfilename = options.outpath+options.folderID+"/Ibis/BWA/s_"+lane+"_sequence"+("_"+bwa_params[lane].replace(" ","_")+"_"+bwa_genomes[lane].strip('/').split('/')[-1]).replace("__","_")+".bam"
    if os.path.exists(outfilename):
      if not options.mock: os.remove(infilename)
      handle_jobs('ln -s %s %s'%(outfilename,infilename))


if (not options.skipBustard) and options.FastQ:
  wait_jobs()
  step+=1
  print ""
  print str(step)+". FastQ conversion of Ibis reads..."
  for lane in (FastQ - skipBustard):
    infilename = options.outpath+options.folderID+"/Bustard/Final_Sequences/s_"+lane+"_sequence.bam"
    if lane in skipFinal: infilename = options.outpath+options.folderID+"/Bustard/Raw_Sequences/s_"+lane+"_sequence.bam"
    if not os.path.exists(infilename) and not options.mock:
      print "Error: No input files for FastQ conversion!"
      continue
    conversion_str = BAM2FastQ
    if isIndexRun: conversion_str += " --RG"
    conversion_str += " --outprefix=%s -o %s %s"%("s_"+lane+"_sequence",options.outpath+options.folderID+"/Bustard/FastQ/",infilename)
    handle_jobs(conversion_str)
  wait_jobs()


if (not options.skipIbis) and options.FastQ:
  wait_jobs()
  step+=1
  print ""
  print str(step)+". FastQ conversion of Ibis reads..."
  for lane in (FastQ - skipIbis):
    infilename = options.outpath+options.folderID+"/Ibis/Final_Sequences/s_"+lane+"_sequence.bam"
    if lane in skipFinal: infilename = options.outpath+options.folderID+"/Ibis/Raw_Sequences/s_"+lane+"_sequence.bam"
    if not os.path.exists(infilename) and not options.mock:
      print "Error: No input files for FastQ conversion!"
      continue
    conversion_str = BAM2FastQ
    if isIndexRun: conversion_str += " --RG" 
    conversion_str += " --outprefix=%s -o %s %s"%("s_"+lane+"_sequence",options.outpath+options.folderID+"/Ibis/FastQ/",infilename)
    handle_jobs(conversion_str)
  wait_jobs()


print ""
print "Waiting for last jobs to finish..."
wait_jobs()

print ""
print "Making sure everyone has read/write/execute rights on output folder..."
handle_jobs("chmod a+rwX -R -f "+options.outpath+options.folderID)
wait_jobs()

dev_null.close()
