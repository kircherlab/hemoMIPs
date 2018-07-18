#!/usr/bin/env python
# -*- coding: ASCII -*-

import sys,os
import pysam
from optparse import OptionParser
import string
import math
import gzip
import subprocess

import xml.dom.minidom
from xml.dom.minidom import Node

parser = OptionParser(usage="usage: %prog RUNFOLDER")
parser.add_option("-o","--outfolder", dest="outfolder", help="Generate Report Folder in X (def RUNFOLDER/Data/Intensities/BaseCalls/RTAreport/)",default="RUNFOLDER/Data/Intensities/BaseCalls/RTAreport/")
(options, args) = parser.parse_args()

nr2base_int = string.maketrans('0123','ACGT')
nr2base_focus = string.maketrans('1234','ACGT')
index_length = 10

#---------------------------
# RTA OUTPUT FILE READERS
#---------------------------

def sum_none(vlist):
  res = 0
  for elem in vlist:
    if elem != None: res += elem
  return res

def ave_none(vlist):
  res = 0
  count = 0
  for elem in vlist:
    if elem != None: 
      res += elem
      count += 1
  if count > 0:
    return res/float(count)
  else:
    return res

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

def read_status(filename):
  res = None
  if os.path.exists(filename):
    res = {}
    doc = xml.dom.minidom.parse(filename)
    for node in doc.getElementsByTagName("Software"):
      h = ""
      for node2 in node.childNodes:
        if node2.nodeType == Node.TEXT_NODE:
          h += node2.data
      res['version'] = h
    for node in doc.getElementsByTagName("Configuration"):
      for i in ['CopyAllFiles','CopyImages','DeleteImages','DeleteIntensity','IsPairedEndRun']:
        h = ""
        for node2 in node.getElementsByTagName(i):
          for node3 in node2.childNodes:
            if node3.nodeType == Node.TEXT_NODE:
              h += node3.data
        h=h.strip()
        if len(h) > 0: res[i] = bool(eval(h.replace("true","True").replace("false","False")))
  return res

def read_stats_read(filename):
  res = None
  if os.path.exists(filename):
    res = {}
    doc = xml.dom.minidom.parse(filename)
    for node in doc.getElementsByTagName("Summary"):
      res['Tile2SqMM'] = float(node.getAttribute("densityRatio"))
    for node in doc.getElementsByTagName("Lane"):
      lane = int(node.getAttribute('key'))
      res[lane]={}
      for i in ['ClustersRaw','ClustersRawSD','ClustersPF','ClustersPFSD','PrcPFClusters','PrcPFClustersSD','Phasing','Prephasing','PrcAlign','PrcAlignSD','FirstCycleIntPF','FirstCycleIntPFSD','PrcIntensityAfter20CyclesPF','PrcIntensityAfter20CyclesPFSD']:
        res[lane][i] = float(node.getAttribute(i))
  return res

def read_error_chart_xml(filename):
  res = None
  if os.path.exists(filename):
    res = {}
    doc = xml.dom.minidom.parse(filename)
    for node in doc.getElementsByTagName("TL"):
      lane,tile = map(int,node.getAttribute("Key").split('_'))
      value = float(node.getAttribute("Val"))
      if lane not in res: res[lane] = {}
      res[lane][tile] = value
  return res

def mean(values,lane):
  if lane in values and len(values[lane]) > 0:
    return sum(values[lane].values())/float(len(values[lane]))
  else:
    return None

def sd(values,lane,mean):
  if lane in values and len(values[lane]) > 0:
    return math.sqrt(sum(map(lambda x:(x-mean)**2,values[lane].values()))/float(len(values[lane])))
  else:
    return None

def det_mean_cycle_error(reads,lanes,tiles,root):
  res = []
  for (start,end) in reads:
    if (end-start+1) > index_length:
      res.append({"Start":start,'End':end,"Errors":{},"ErrorsSD":{},"AveError":{},"AveErrorSD":{}})
      for lane in range(lanes):
        res[-1]["Errors"][lane+1]=[]
        res[-1]["ErrorsSD"][lane+1]=[]
        res[-1]["AveError"][lane+1]=None
        res[-1]["AveErrorSD"][lane+1]=None
      errsum = [0]*lanes
      sdsum = [0]*lanes
      counts = [0]*lanes
      for cycle in range(start,end+1):
        helper = read_error_chart_xml(root+"%d.xml"%cycle)
        if helper != None:
          for lane in range(lanes):
            val = mean(helper,lane+1)
            val2 = sd(helper,lane+1,val)
            if val != None:
              errsum[lane] += val
              sdsum[lane] += val2
              counts[lane] += 1
            res[-1]["Errors"][lane+1].append(val)
            res[-1]["ErrorsSD"][lane+1].append(val2)
        else:
          for lane in range(lanes): res[-1]["Errors"][lane+1].append(None)
      for lane in range(lanes):
        res[-1]["AveError"][lane+1] = None if counts[lane] == 0 else errsum[lane]/counts[lane]
        res[-1]["AveErrorSD"][lane+1] = None if counts[lane] == 0 else sdsum[lane]/counts[lane]
  return res

def read_intensities(filename):
  res = None
  if os.path.exists(filename):
    res = {}
    if filename.endswith('.gz'):
      infile = gzip.open(filename)
    else:
      infile = open(filename)
    first = False
    for line in infile:
      if line.startswith("Color\tCycle\tLane\tTile\t"): 
        first = True
        continue
      elif not first: continue
      fields = line.split()
      if len(fields) == 5:
        base = fields[0].translate(nr2base_int)
        cycle = int(fields[1])
        lane = int(fields[2])+1
        tile = int(fields[3])
        if fields[4].lower() == "nan":
          intensity = None
        else:
          intensity = float(fields[4])
        if not lane in res: res[lane] = {}
        if tile not in res[lane]: res[lane][tile] = { "A":[],"C":[],"G":[],"T":[] }
        res[lane][tile][base].append(intensity)
    infile.close()
  return res

def read_focus(filename):
  res = None
  if os.path.exists(filename):
    res = {}
    if filename.endswith('.gz'):
      infile = gzip.open(filename)
    else:
      infile = open(filename)
    first = False
    for line in infile:
      if line.startswith("Color\tCycle\tLane\tTile\t"): 
        first = True
        continue
      elif not first: continue
      fields = line.split()
      if len(fields) == 5:
        base = fields[0].translate(nr2base_focus)
        cycle = int(fields[1])
        lane = int(fields[2])
        tile = int(fields[3])
        if fields[4].lower() == "nan":
          intensity = None
        else:
          intensity = float(fields[4])
        if not lane in res: res[lane] = {}
        if tile not in res[lane]: res[lane][tile] = { "A":[],"C":[],"G":[],"T":[] }
        res[lane][tile][base].append(intensity)
    infile.close()
  return res

def read_cluster_counts(filename):
  res = None
  if os.path.exists(filename):
    res = {}
    if filename.endswith('.gz'):
      infile = gzip.open(filename)
    else:
      infile = open(filename)
    first = False
    for line in infile:
      if line.startswith("Lane\tTile\tCluster"):
        first = True
        continue
      elif not first: continue
      fields = line.split()
      if len(fields) == 3:
        lane = int(fields[0])+1
        tile = int(fields[1])
        clusters = float(fields[2])
        if not lane in res: res[lane] = []
        res[lane].append(clusters)
    infile.close()
  return res

def read_stats_html(read_stats,info=''):
  if read_stats == None:
    return ""
  else:
    outstr  = "<h3>%s Read:</h3><p>"%info
    outstr += "<table border=1 align='center'>\n<tr><th>Lane</th>"
    headers  = [('Clusters [Tile]','ClustersRaw','ClustersRawSD'),
                ('Clusters PF','ClustersPF','ClustersPFSD'),
                ('PF %','PrcPFClusters','PrcPFClustersSD'),
                ('1st Cycle Int (PF)','FirstCycleIntPF','FirstCycleIntPFSD'),
                ('% Intensity 20cycles','PrcIntensityAfter20CyclesPF','PrcIntensityAfter20CyclesPFSD'),
                ('Phasing','Phasing'),
                ('Prephasing','Prephasing')]
    for cols in headers:
      outstr += "<th>%s</th>"%cols[0]
    outstr += "</tr>\n"
    lanes = filter(lambda x: type(x) == type(1),read_stats.keys())
    for lane in range(min(lanes),max(lanes)+1):
      outstr += "<tr><td>%d</td>"%lane
      for cols in headers:
        if len(cols) == 3:
          if read_stats[lane][cols[1]] > 100:
            outstr += "<td align='right'>%.0f &plusmn;%.2f</td>"%(read_stats[lane][cols[1]],read_stats[lane][cols[2]])
          else:
            outstr += "<td align='right'>%.2f &plusmn;%.2f</td>"%(read_stats[lane][cols[1]],read_stats[lane][cols[2]])
        elif len(cols) == 2:
          outstr += "<td align='right'>%.4f</td>"%(read_stats[lane][cols[1]])
      outstr += "</tr>\n"
    outstr += "</table></p>"
    return outstr

#---------------------------
# REPORT WRITER FUNCTIONS
#---------------------------

def error_stats_html(errors,read_stats,info=''):
  if read_stats == None:
    return ""
  else:
    rrange = (errors['Start'],errors['End'])
    outstr  = "<h3>Error of controls for %s Read (%d-%d):</h3>"%(info,rrange[0],rrange[1])
    outstr += "<p><table border=1 align='center'>\n<tr><th>Lane</th>"
    headers  = [('% Align (PF)','PrcAlign','PrcAlignSD')]
    for cols in headers:
      outstr += "<th>%s</th>"%cols[0]
    middle = (rrange[1]-rrange[0]+1)/2 - 1
    for i in ["Error 1st base","%dth base"%(middle+1),"Last base","Average Error"]:
      outstr += "<th>%s</th>"%i
    outstr += "</tr>\n"

    lanes = filter(lambda x: type(x) == type(1),read_stats.keys())
    for lane in range(min(lanes),max(lanes)+1):
      outstr += "<tr><td>%d</td>"%lane

      for cols in headers:
        if len(cols) == 3:
          if read_stats[lane][cols[1]] > 100:
            outstr += "<td align='right'>%.0f &plusmn;%.2f</td>"%(read_stats[lane][cols[1]],read_stats[lane][cols[2]])
          else:
            outstr += "<td align='right'>%.2f &plusmn;%.2f</td>"%(read_stats[lane][cols[1]],read_stats[lane][cols[2]])
        elif len(cols) == 2:
          outstr += "<td align='right'>%.4f</td>"%(read_stats[lane][cols[1]])

      for i in [0,middle,-1]:
        value = errors['Errors'][lane][i]
        if value != None: outstr += "<td align='right'>%.2f%% &plusmn;%.2f</td>"%(value,errors['ErrorsSD'][lane][i])
        else: outstr += "<td align='center'>NA</td>"

      if errors['AveError'][lane] != None: outstr += "<td align='right'>%.2f%% &plusmn;%.2f</td>"%(errors['AveError'][lane],errors['AveErrorSD'][lane])
      else: "<td align='center'>NA</td>"

      outstr += "</tr>\n"
    outstr += "</table></p>"
    return outstr

def generate_image_error(errors):
  rcmd  = "x <- c(%s)\n"%(",".join(map(str,reduce(lambda x,y:x+y,map(lambda x: range(x['Start'],x['End']+1),errors)))))
  rcmd += "y <- c()\n"
  maxlanes = len(errors[0]['Errors'].keys())
  for lane in errors[0]['Errors'].keys():
    rcmd += "y <- cbind(y,c(%s))\n"%(",".join(map(lambda x: 'NA' if x == None else str(x),reduce(lambda x,y:x+y,map(lambda x:x['Errors'][lane],errors)))))
  rcmd += "png('"+options.outfolder+"/images/error_lanes.png',width=800,height=400)\n"
  rcmd += "matplot(x,y,xlab='Cycle',main='Per cycle error rate of control reads',ylab='Error rate [%%]',type='l',lty=1,lwd=2,col=1:%d)\n"%(maxlanes)
  rcmd += "legend('topleft',sprintf('Lane %s',1:%d),fill=1:%d)\n"%("%d",maxlanes,maxlanes)
  rcmd += "dev.off()\n"
  R = subprocess.Popen("R --vanilla --quiet",stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
  R.communicate(rcmd)
  #print rcmd

def cluster_int_focus_tbl(clusters,PFclusters,intensities,focus,cfactor):
  lanes = clusters.keys()
  lanes.sort()
  outstr = ""
  for lane in lanes:
    outstr += "<h3>Lane %d</h3>\n"%lane
    outstr += "<p><table border=1 align='center'>\n"
    outstr += "<tr><th>Lane</th><th>Tile</th><th>Clusters</th><th>Clusters PF</th><th>% PF</th><th>1st Cycle Int</th><th>20th Cycle Int</th><th>% Intensity</th><th>Min Focus</th><th>Median Focus</th><th>Max Focus</th></tr>\n"
    nrtiles = len(clusters[lane])
    tiles = intensities[lane].keys()
    tiles.sort()
    for ind,tile in enumerate(tiles):
      outstr += "<tr><td>%d</td><td>%d</td><td align='right'>%.0f</td><td align='right'>%.0f</td><td align='right'>%.2f%%</td>"%(lane,tile,clusters[lane][ind]/cfactor,PFclusters[lane][ind]/cfactor,(PFclusters[lane][ind]*100.0)/clusters[lane][ind])
      first_int = ave_none([intensities[lane][tile]['A'][0],intensities[lane][tile]['C'][0],intensities[lane][tile]['G'][0],intensities[lane][tile]['T'][0]])
      twenty_int = ave_none([intensities[lane][tile]['A'][0],intensities[lane][tile]['C'][0],intensities[lane][tile]['G'][20],intensities[lane][tile]['T'][20]])
      helper = focus[lane][tile]['A']+focus[lane][tile]['C']+focus[lane][tile]['G']+focus[lane][tile]['T']
      helper.sort()
      if len(helper) % 2 == 0 and helper[len(helper) / 2] != None: focusval = helper[len(helper) / 2]
      else: focusval = ave_none([helper[len(helper) / 2],helper[(len(helper)+1) / 2]])
      outstr += "<td align='right'>%.2f</td><td align='right'>%.2f</td><td align='right'>%.2f</td><td align='right'>%.2f</td><td align='right'>%.2f</td><td align='right'>%.2f</td></tr>\n"%(first_int,twenty_int,(twenty_int*100.0)/first_int if first_int > 0 else 0,0 if helper[0] == None else helper[0],focusval,0 if helper[-1] == None else helper[-1])
    outstr += "</table></p>"
  return outstr

def generate_images_int_focus(ctype,values):
  lanes = values.keys()
  lanes.sort()
  tiles = values[lanes[0]].keys()
  tiles.sort()
  max_cycle = len(values[lanes[0]][tiles[0]]['A'])
  rcmd  = "x <- c(%s)\n"%(",".join(map(str,range(1,max_cycle+1))))
  outstr = "<h2>%s</h2>"%(ctype[0].upper()+ctype[1:])
  for lane in lanes:
    outstr += "<h3>Lane %d</h3>\n"%(lane)
    outstr += "<p align='center'><a href='images/%s_lane%d.png'><img src='images/%s_lane%d.png'></a></p>\n"%(ctype,lane,ctype,lane)
    for base in "ACGT":
      clist = []
      for cycle in range(max_cycle):
        val,count = 0.0,0
        for tile in tiles:
          if values[lane][tile][base][cycle] != None:
            val += values[lane][tile][base][cycle]
            count += 1
        if count > 0: clist.append(str(val/float(count)))
        else: clist.append("NA")
      rcmd += "%s <- c(%s)\n"%(base,",".join(clist))
    rcmd += "png('"+options.outfolder+"/images/%s_lane%d.png',width=800,height=400)\n"%(ctype,lane)
    rcmd += "matplot(x,cbind(A,C,G,T),xlab='Cycle',main='Per cycle %s values (Lane %d)',ylab='%s',type='l',lty=1,lwd=2,col=c('green','blue','black','red'))\n"%(ctype,lane,ctype[0].upper()+ctype[1:])
    rcmd += "legend('topright',c('A','C','G','T'),fill=c('green','blue','black','red'))\n"
    rcmd += "dev.off()\n"
  #if ctype.startswith("focus"): print rcmd
  R = subprocess.Popen("R --vanilla --quiet",stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
  R.communicate(rcmd)
  #print rcmd
  return outstr

#---------------------------
# MAIN
#---------------------------

if (len(args) == 1) and os.path.isdir(args[0]):
  print "Evaluating RunInfo file..."
  run_info = read_runinfo(args[0]+"/RunInfo.xml")
  if run_info == None:
    sys.stderr.write("Minimum input files missing for report.\n")
    sys.exit()

  print "Evaluating Status file..."
  run_status = read_status(args[0]+"/Data/reports/Status.xml")
  if run_status == None:
    sys.stderr.write("Minimum input files missing for report.\n")
    sys.exit()

  print "Checking for output folder..."
  options.outfolder = options.outfolder.replace("RUNFOLDER",args[0])
  if not os.path.isdir(options.outfolder): os.makedirs(options.outfolder)
  if not os.path.isdir(options.outfolder+'/images'): os.mkdir(options.outfolder+'/images')

  print "Evaluating Read Summary files..."
  read_stats_forward = read_stats_read(args[0]+"/Data/reports/Summary/read1.xml")
  if 'IsPairedEndRun' in run_status and run_status['IsPairedEndRun']:
    if os.path.exists(args[0]+"/Data/reports/Summary/read3.xml"):
      read_stats_reverse = read_stats_read(args[0]+"/Data/reports/Summary/read3.xml")
    else:
      read_stats_reverse = read_stats_read(args[0]+"/Data/reports/Summary/read2.xml")
  else:
    read_stats_reverse = None
  #print str(read_stats_forward)[:100]
  #print str(read_stats_reverse)[:100]
  forward_stats_text = read_stats_html(read_stats_forward,'Forward')
  reverse_stats_text = read_stats_html(read_stats_reverse,'Reverse')

  print "Copying some images from RTA Summary..."
  pf_image = ""
  if os.path.exists(args[0]+"/Data/reports/NumClusters By Lane.png"):
    outfile = open(options.outfolder+'/images/pf_stats.png','w')
    infile = open(args[0]+"/Data/reports/NumClusters By Lane.png")
    outfile.write(infile.read())
    outfile.close()
    infile.close()
    pf_image = "<p align='center'><a href='images/pf_stats.png'><img src='images/pf_stats.png'></a></p>\n"
  else:
    print "Missing cluster/cluster PF overview figure. Skipping part of the report."
    print args[0]+"/Data/reports/NumClusters By Lane.png",os.path.exists(args[0]+"/Data/reports/NumClusters By Lane.png")

  flowcell_overview_images = ""
  if (os.path.exists(args[0]+"/Data/reports/NumClusters_Chart.png") and 
      os.path.exists(args[0]+"/Data/reports/NumPassedFilter25_Chart.png") and 
      os.path.exists(args[0]+"/Data/reports/PassedFilter25_Chart.png")):
    outfile = open(options.outfolder+'/images/num_cl_fc.png','w')
    infile = open(args[0]+"/Data/reports/NumClusters_Chart.png")
    outfile.write(infile.read())
    outfile.close()
    infile.close()
    outfile = open(options.outfolder+'/images/num_cl_pf_fc.png','w')
    infile = open(args[0]+"/Data/reports/NumPassedFilter25_Chart.png")
    outfile.write(infile.read())
    outfile.close()
    infile.close()
    outfile = open(options.outfolder+'/images/frac_cl_pf_fc.png','w')
    infile = open(args[0]+"/Data/reports/PassedFilter25_Chart.png")
    outfile.write(infile.read())
    outfile.close()
    infile.close()
    flowcell_overview_images = """<p align='center'>
    <a href='images/num_cl_fc.png'><img src='images/num_cl_fc.png'></a>
    <a href='images/num_cl_pf_fc.png'><img src='images/num_cl_pf_fc.png'></a>
    <a href='images/frac_cl_pf_fc.png'><img src='images/frac_cl_pf_fc.png'></a>
    </p>\n"""
  else:
    print "Missing at least one of the overview images. Skipping part of the report."
    print args[0]+"Data/reports/NumClusters_Chart.png",os.path.exists(args[0]+"/Data/reports/NumClusters_Chart.png")
    print args[0]+"Data/reports/NumPassedFilter25_Chart.png",os.path.exists(args[0]+"/Data/reports/NumPassedFilter25_Chart.png")
    print args[0]+"Data/reports/PassedFilter25_Chart.png",os.path.exists(args[0]+"/Data/reports/PassedFilter25_Chart.png")

  forward_error_text = ""
  reverse_error_text = ""
  error_image = ""
  print "Determining average error rates..."
  errors = det_mean_cycle_error(run_info['read_ranges'],run_info['lanes'],run_info['tiles'],args[0]+"/Data/reports/ErrorRate/Chart_")
  #print str(errors)[:100]
  forward_error_text = error_stats_html(errors[0],read_stats_forward,'Forward')
  reverse_error_text = error_stats_html(errors[-1],read_stats_reverse,'Reverse')

  generate_image_error(errors)
  error_image = ""
  if os.path.exists(options.outfolder+'/images/error_lanes.png'):
    error_image = "<p align='center'><a href='images/error_lanes.png'><img src='images/error_lanes.png'></a></p>"

  print "Reading tile intensity values..."
  intensities = read_intensities(args[0]+"/Data/reports/Intensity By Color And Cycle.txt")
  #print str(intensities)[:100]

  print "Reading tile focus values..."
  focus = read_focus(args[0]+"/Data/reports/FWHM By Color And Cycle.txt")
  #print str(focus)[:100]

  intensities_dev_images = generate_images_int_focus('intensity',intensities)
  #print intensities_dev_images
  focus_dev_images =generate_images_int_focus('focus',focus)
  #print focus_dev_images

  print "Reading tile cluster values..."
  clusters = read_cluster_counts(args[0]+"/Data/reports/NumClusters By Lane.txt")
  #print str(clusters)[:100]

  print "Reading tile cluster PF values..."
  PFclusters = read_cluster_counts(args[0]+"/Data/reports/NumClusters By Lane PF.txt")
  #print str(PFclusters)[:100]

  table_txt = cluster_int_focus_tbl(clusters,PFclusters,intensities,focus,read_stats_forward['Tile2SqMM'])

  print "Generating output HTML..."
  html = """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
  <html><body>
  <h1 align="center">%s Report<br>%s</h1>
  <p align="center">Lanes: %d, Tiles: %d, Tiles per mm&sup2;: %.4f, Total number of cycles: %d</p>
  <h2>Run Information:</h2><p>%s</p>
  <h3>Reads:</h3><p>%s</p>
  %s%s
  <h3>Sequencing error:</h3>
  %s
  %s%s
  </body></html>
  """%(run_status['version'],run_info['runid'],
  run_info['lanes'],run_info['tiles'],read_stats_forward['Tile2SqMM'],run_info['cycles'],
  ", ".join(map(lambda (x,y):"%s : %s"%(x,y),filter(lambda (x,y): x!='version',run_status.iteritems()))),
  "<br>".join(map(lambda (x,y): "Read %d%s: %d - %d"%(x+1," [Index]" if (y[1]-y[0]) < index_length else "",y[0],y[1]),enumerate(run_info['read_ranges']))),
  forward_stats_text,reverse_stats_text,
  error_image,
  forward_error_text,reverse_error_text)

  outfile = open(options.outfolder+"/Summary_short.htm",'w')
  outfile.write(html)
  outfile.close()

  html = """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
  <html><body>
  <h1 align="center">%s Report<br>%s</h1>
  <p align="center">Lanes: %d, Tiles: %d, Tiles per mm&sup2;: %.4f, Total number of cycles: %d</p>
  <h2>Run Information:</h2><p>%s</p>
  <h3>Reads:</h3><p>%s</p>
  %s
  %s%s
  %s
  <h3>Sequencing error:</h3>
  %s
  %s%s
  %s
  %s
  <h2>Tile statistics</h2>
  %s
  </body></html>
  """%(run_status['version'],run_info['runid'],
  run_info['lanes'],run_info['tiles'],read_stats_forward['Tile2SqMM'],run_info['cycles'],
  ", ".join(map(lambda (x,y):"%s : %s"%(x,y),filter(lambda (x,y): x!='version',run_status.iteritems()))),
  "<br>".join(map(lambda (x,y): "Read %d%s: %d - %d"%(x+1," [Index]" if (y[1]-y[0]) < index_length else "",y[0],y[1]),enumerate(run_info['read_ranges']))),
  pf_image,
  forward_stats_text,reverse_stats_text,
  flowcell_overview_images,
  error_image,
  forward_error_text,reverse_error_text,
  intensities_dev_images,
  focus_dev_images,
  table_txt)

  outfile = open(options.outfolder+"/Summary.htm",'w')
  outfile.write(html)
  outfile.close()

else:
  sys.stderr.write('Expect exactly one argument: valid path to RUNFOLDER.\n')
  sys.exit()
