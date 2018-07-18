#!/usr/bin/env python
# -*- coding: ASCII -*-

import sys,os
import numpy
from optparse import OptionParser
import gzip
import StringIO

sys.path.append('/net/shendure/vol1/home/mkircher/bin/')
from library import parse_rangestr

parser = OptionParser()
parser.add_option("-p", "--path", dest="path", help="Path to base caller results (default .)", default=".")
parser.add_option("-e", "--expID", dest="expID", help="Name of experiment to use for sequence names")
parser.add_option("-l", "--lanes", dest="lanes", help="Lanes, example: 1,3,4-8 (default 1-8)",default="1-8")
parser.add_option("-t", "--tiles", dest="tiles", help="Tiles, example: 1-13,100 (default all)")
parser.add_option("-o", "--outpath", dest="outdir", help="Path for output files (default .)",default=".")
parser.add_option("-s", "--PIPE", dest="PIPE", help="Write FastQ output to stdout and not to file (default off)",default=False,action="store_true")
parser.add_option("-c", "--coordianteType", dest="coordtype", help="Type of cluster coordinate transformation: round = Pipeline <= 1.4: round(x), floor = Pipeline 1.5.*: floor(abs(x)), shift_floor (default) = Pipeline 1.6.*: round(x*10+1000)",choices=["round","floor","shift_round"],default="shift_round")
parser.add_option("--min_cycle", dest="min_cycle", help="Fist cycle to include in conversion (default 1)",default=1,type="int")
parser.add_option("--max_cycle", dest="max_cycle", help="Last cycle to include in conversion (default MAX)",default=None,type="int")
parser.add_option("--skip_cycles", dest="skip_cycles", help="Skip cycles (default '')",default="")
parser.add_option("-b", "--allow_broken_cycles",dest="allow_broken_cycles", help="Replace cycles with broken BCL files by N bases/! qualities. (Def: skip whole tile)",default=False,action="store_true")
parser.add_option("-v", "--verbose",dest="verbose", help="Activate verbose output of messages on STDERR",default=False,action="store_true")
(options, args) = parser.parse_args()

pskipCycles = parse_rangestr(options.skip_cycles)
skipCycles = set() if pskipCycles == None else set(pskipCycles)

if options.verbose and len(skipCycles) != 0:
  sys.stderr.write('Excluding cycles: %s\n'%(str(skipCycles)))

lendian = 1 ## use -1 for big endian systems
bases = 'ACGT'

def coord_round(x): return "%d"%(round(float(x)))
def coord_floor(x): return "%d"%(math.floor(abs(float(x))))
def coord_shift_round(x): return "%d"%(round(float(x)*10+1000))
coordconv = coord_shift_round

def to_int(s):
  global lendian
  power = 1
  total = 0
  for elem in s[::lendian]:
    total += ord(elem)*power
    power *= 256
  return total

def read_position_txt(filename):
  global coordconv
  if filename.endswith('.gz'):
    infile = gzip.open(filename)
  else:
    infile = open(filename)
  for line in infile:
    fields = map(float,line.split())
    yield coordconv(fields[0]),coordconv(fields[1])
  infile.close()
  raise StopIteration

def read_locs(filename):
  global coordconv
  data = []
  if filename.endswith('.gz'):
    infile = gzip.open(filename,'rb')
    infile.read(8) # First 8 Byte are unused
    clusters = to_int(infile.read(4))
    shape = (clusters,2)
    data = numpy.fromstring(infile.read(), dtype=numpy.float32).reshape(shape)
  else:
    infile = open(filename,'rb')
    infile.read(8) # First 8 Byte are unused
    clusters = to_int(infile.read(4))
    shape = (clusters,2)
    data = numpy.fromfile(file=infile, dtype=numpy.float32).reshape(shape)
  for x,y in data:
    yield coordconv(x),coordconv(y)
  infile.close()
  raise StopIteration

def read_clocs(filename):
  EXPECTED_CLOCS_VERSION = 1
  BLOCK_SIZE = 25
  IMAGE_WIDTH = 2048
  BLOCKS_PER_LINE = (IMAGE_WIDTH + BLOCK_SIZE - 1) / BLOCK_SIZE
  totalBlocks = 0
  currentBlock = 0
  currentBlockUnreadClusters = 0

  if filename.endswith('.gz'):
    infile = gzip.open(filename,'rb')
  else:
    infile = open(filename,'rb')
  clocsVersion = ord(infile.read(1))
  totalBlocks = to_int(infile.read(4))
  currentBlockUnreadClusters = ord(infile.read(1))
  currentBlock+=1

  while (currentBlock < totalBlocks or ( currentBlock == totalBlocks and currentBlockUnreadClusters > 0)):
     while (currentBlockUnreadClusters == 0 and currentBlock < totalBlocks):
        currentBlockUnreadClusters = ord(infile.read(1))
        currentBlock += 1
     dx = ord(infile.read(1))
     dy = ord(infile.read(1))
     x = 10 * BLOCK_SIZE * ((currentBlock - 1) % BLOCKS_PER_LINE) + dx + 1000;
     y = 10 * BLOCK_SIZE * ((currentBlock - 1) / BLOCKS_PER_LINE) + dy + 1000;
     yield x,y
     currentBlockUnreadClusters -= 1
  infile.close()
  raise StopIteration


def read_bcl(filename,clusteridx=None):
  global bases,options
  if filename.endswith('.gz') or filename.endswith(".bgzf"):
    infile = gzip.open(filename,'rb')
  else:
    infile = open(filename,'rb')
  if options.verbose:
    sys.stderr.write("Reading BCL file: %s\n"%(filename))
  try:
    nrclusters = to_int(infile.read(4))
    data = infile.read()
    if len(data) == nrclusters:
      iterclusters = xrange(nrclusters)
      if clusteridx != None:
        iterclusters = clusteridx
      for cluster in iterclusters:
        base = (ord(data[cluster]) & 3)
        quality = ord(data[cluster]) >> 2
        if quality != 0: 
          yield bases[base],quality
        else:
          yield "N",quality
    else:
      sys.stderr.write("File content does not reflect number of clusters (%d vs %d).\n"%(len(data),nrclusters))
      if len(data) < nrclusters:
        sys.stderr.write("Trying to continue...\n")
        nrclusters = len(data)
        iterclusters = xrange(nrclusters)
        if clusteridx != None:
          iterclusters = clusteridx
        for cluster in iterclusters:
          base = (ord(data[cluster]) & 3)
          quality = ord(data[cluster]) >> 2
          if quality != 0: 
            yield bases[base],quality
          else:
            yield "N",quality
        
  except IOError:
    sys.stderr.write("Error reading: %s\n"%filename)
  infile.close()
  raise StopIteration

def get_bcl_cycles(root_path,lane,tile,start,end,clusteridx=None,qualities=True):
  global skipCycles
  lanes = map(lambda x: "L%03d"%(x+1),range(8))
  res_bases = None
  res_quals = None
  prefix_bases,prefix_qualities = "",""
  if lane > 0 and lane <= len(lanes):
    cycles = map(lambda x: ("C%d.1"%x,x),range(start,end+1))
    if os.path.isdir(root_path+"/"+lanes[lane-1]):
      res_bases = []
      if qualities: res_quals = []
      for ind,(cycle,cycleVal) in enumerate(cycles):
        if cycleVal in skipCycles: continue
        fullname = root_path+"/"+lanes[lane-1]+"/"+cycle+"/s_%d_%d.bcl"%(lane,tile)
        if os.path.isfile(root_path+"/"+lanes[lane-1]+"/"+cycle+"/s_%d_%d.bcl.gz"%(lane,tile)):
          fullname = root_path+"/"+lanes[lane-1]+"/"+cycle+"/s_%d_%d.bcl.gz"%(lane,tile)
        elif os.path.isfile(root_path+"/"+lanes[lane-1]+"/%04d.bcl.bgzf"%(cycleVal)):
          fullname = root_path+"/"+lanes[lane-1]+"/%04d.bcl.bgzf"%(cycleVal)                    
        if os.path.isfile(fullname):
          count = 0
          for cluster in read_bcl(fullname,clusteridx):
            if ind == 0:
              res_bases.append(prefix_bases+cluster[0])
              if qualities: res_quals.append(prefix_qualities+chr(cluster[1]+33))
            else:
              res_bases[count]+=cluster[0]
              if qualities: res_quals[count]+=chr(cluster[1]+33)
            count+=1
          if count == 0:
            sys.stderr.write("Error reading %s\n"%(fullname))
            if options.allow_broken_cycles:
              if ind == 0:
                prefix_bases = 'N'
                prefix_qualities = '!'
              elif len(res_bases) == 0:
                prefix_bases += 'N'
                prefix_qualities += '!'
              else:
                for ccount,value in enumerate(res_bases):
                  res_bases[ccount]+='N'
                  if qualities: res_quals[ccount]+='!'
            else:
              return None,None
        else:
          sys.stderr.write("Error: Could not find %s\n"%(root_path+"/"+lanes[lane-1]+"/"+cycle))
    else:
      sys.stderr.write("Error: Could not find %s\n"%(root_path+"/"+lanes[lane-1]))
      return None,None
  return res_bases,res_quals

if options.expID == None:
  sys.stderr.write("Need a name for the experiment.\n")
  sys.exit()

path = options.path.rstrip('/')+'/'
outpath = options.outdir.rstrip('/')+'/'
expID = options.expID
start = options.min_cycle
end = options.max_cycle

lanes = parse_rangestr(options.lanes)
if lanes == None:
  sys.stderr.write("Need a valid range of lanes.\n")
  sys.exit()

if options.tiles <> None:
  tiles = parse_rangestr(options.tiles)
  if tiles:
      sys.stderr.write("Will use only tiles: %s\n"%str(tiles))
else:
  tiles = None

if options.coordtype == "round":  coordconv = coord_round
elif options.coordtype == "floor": coordconv = coord_floor
else: coordconv = coord_shift_round

for lane in lanes:
  if not options.PIPE:
    outfile = open(outpath+"s_%s_sequence.txt"%(lane),'w')
  else:
    outfile = sys.stdout
  if end == None:
    cycles = []
    try:
      cycles = map(lambda x: int(x[1:].split(".")[0]),filter(lambda x: x.startswith('C') and os.path.isdir(path+"L00%s/"%(lane)+x), os.listdir(path+"L00%s/"%lane)))
    except:
      sys.stderr.write("Error: No valid BaseCall input folder provided!\n")
      sys.exit()
    if len(cycles) == 0:
      try:
        cycles = map(lambda x: int(x.split(".")[0]),filter(lambda x: x.endswith('.bcl.bgzf') and not os.path.isdir(path+"L00%s/"%(lane)+x), os.listdir(path+"L00%s/"%lane)))
      except:
        sys.stderr.write("Error: No valid BaseCall input folder provided!\n")
        sys.exit()
    cycles.sort()
    end = cycles[-1]
    if options.verbose: sys.stderr.write("Determined run length with %d cycles.\n"%(end))
  if tiles == None:
    try:
      tiles = map(lambda x: int(x.split('.')[0].split('_')[2]),filter(lambda x:'.bcl' in x,os.listdir(path+"L00%s/C%d.1/"%(lane,options.min_cycle))))
    except:
      sys.stderr.write("Error: No valid BaseCall input folder provided!\n")
      sys.exit()
    tiles.sort()
    tiles = map(str,tiles)
    if options.verbose: sys.stderr.write("Found a total of %d tiles.\n"%(len(tiles)))
  for tile in tiles:
    if os.path.exists(path+"/../s_%s_%04d_pos.txt"%(lane,int(tile))):
      if options.verbose: sys.stderr.write("Opening position file: %s\n"%(path+"/../s_%s_%04d_pos.txt"%(lane,int(tile))))
      positionfile = read_position_txt(path+"/../s_%s_%04d_pos.txt"%(lane,int(tile)))
    elif os.path.exists(path+"/../s_%s_%04d_pos.txt.gz"%(lane,int(tile))):
      if options.verbose: sys.stderr.write("Opening position file: %s\n"%(path+"/../s_%s_%04d_pos.txt.gz"%(lane,int(tile))))
      positionfile = read_position_txt(path+"/../s_%s_%04d_pos.txt.gz"%(lane,int(tile)))
    elif os.path.exists(path+"/../L00%s/s_%s_%d.locs"%(lane,lane,int(tile))):
      if options.verbose: sys.stderr.write("Opening position file: %s\n"%(path+"/../L00%s/s_%s_%d.locs"%(lane,lane,int(tile))))
      positionfile = read_locs(path+"/../L00%s/s_%s_%d.locs"%(lane,lane,int(tile)))
    elif os.path.exists(path+"/../L00%s/s_%s.locs"%(lane,lane)):
      if options.verbose: sys.stderr.write("Opening position file: %s\n"%(path+"/../L00%s/s_%s.locs"%(lane,lane)))
      positionfile = read_locs(path+"/../L00%s/s_%s_%d.locs"%(lane,lane,int(tile)))
    elif os.path.exists(path+"/../L00%s/s_%s_%d.locs.gz"%(lane,lane,int(tile))):
      if options.verbose: sys.stderr.write("Opening position file: %s\n"%(path+"/../L00%s/s_%s_%d.locs.gz"%(lane,lane,int(tile))))
      positionfile = read_locs(path+"/../L00%s/s_%s_%d.locs.gz"%(lane,lane,int(tile)))
    elif os.path.exists(path+"/../L00%s/s_%s_%d.clocs"%(lane,lane,int(tile))):
      if options.verbose: sys.stderr.write("Opening position file: %s\n"%(path+"/../L00%s/s_%s_%d.clocs"%(lane,lane,int(tile))))
      positionfile = read_clocs(path+"/../L00%s/s_%s_%d.clocs"%(lane,lane,int(tile)))
    elif os.path.exists(path+"/../L00%s/s_%s_%d.clocs.gz"%(lane,lane,int(tile))):
      if options.verbose: sys.stderr.write("Opening position file: %s\n"%(path+"/../L00%s/s_%s_%d.clocs.gz"%(lane,lane,int(tile))))
      positionfile = read_clocs(path+"/../L00%s/s_%s_%d.clocs.gz"%(lane,lane,int(tile)))
    else:
      sys.stderr.write("Did not find a cluster position file\n")
      positionfile = None
    #count = 0
    #for coords in positionfile:
      #count += 1
    #print count
    if options.verbose: sys.stderr.write("Reading BCL files for this tile...\n")
    reads,qualities = get_bcl_cycles(path,int(lane),int(tile),start,end,None,True)
    if options.verbose: sys.stderr.write("Generating FastQ output for this tile...\n")
    if reads != None:
      lcoords = None,None
      for ind,seq in enumerate(reads):
        if positionfile != None:
          try:
            coords = positionfile.next()
          except:
            sys.stderr.write("Error: More clusters than positions. Last coords were lane %d tile %s %s!\n"%(lane,tile,str(lcoords)))
            sys.stderr.write("Continuing with truncated sequence file!\n")
            break
        else:
          coords = ind+1,ind+1
        lcoords = coords
        outfile.write("@%s:%s:%s:%s:%s\n%s\n+\n%s\n"%(options.expID,lane,tile,coords[0],coords[1],seq,qualities[ind]))
  outfile.close()
