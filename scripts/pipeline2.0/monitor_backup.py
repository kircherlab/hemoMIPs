#!/usr/bin/env python

"""
Periodically scans for old finished runs and backups them.
Runs as daemon

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *30.02.2009
:type: daemon
:input: -
:output: -

options:
    -b, --backup        Path to backup storage (default /raid01/backup/)
    -d, --data          Path to data storage (default /raid01/data/)
    -o, --old           Path to old run storage (default /raid01/old_runs/)
    -m, --mock          Make a test run

"""

import sys
import os
from optparse import OptionParser
import time
import subprocess
import smtplib
from email.mime.text import MIMEText

time_wait = 15*60
old_day_wait = 7.0*6 # Number of days before run is moved to old runs
image_day_wait = 7.0*3 # Number of days before images are removed after analysis

timeout_backup_file = 60 # Number of minuts before backup creation is considered dead

sender = 'sbsuser@bioaps01'
recipient = 'martin_kircher@eva.mpg.de,schultz@eva.mpg.de,sequencing@eva.mpg.de,kelso@eva.mpg.de,udo_stenzel@eva.mpg.de'

parser = OptionParser()
parser.add_option("-b", "--backup", dest="backup", help="Path to backup storage (default /raid01/backup/)",default="/raid01/backup/")
parser.add_option("-d", "--data", dest="data", help="Path to data storage (default /raid01/data/)",default="/raid01/data/")
parser.add_option("-o", "--old", dest="old", help="Path to old run storage (default /raid01/old_runs/)",default="/raid01/old_runs/")
parser.add_option("-m", "--mock", dest="mock", help="Make a test run",default=False,action="store_true")
(options, args) = parser.parse_args()

data_folder = options.data.rstrip('/')+'/'
backup_folder = options.backup.rstrip('/')+'/'
backup_info_file = "backup_history.txt.tmp"
old_folder = options.old.rstrip('/')+'/'

if not os.path.isdir(data_folder):
  print "Need valid path for data folder:",data_folder
  sys.exit()
if not os.path.isdir(backup_folder):
  print "Need valid path for backup folder:",backup_folder
  sys.exit()
if not os.path.isdir(old_folder):
  print "Need valid path for old run folder:",old_folder
  sys.exit()
if not os.path.exists(backup_folder+backup_info_file):
  print "Could not find backup info file",backup_folder+backup_info_file
  sys.exit()

def get_diskfree(path):
  proc = subprocess.Popen('/bin/df -h '+path,shell=True,stdout=subprocess.PIPE)
  proc.wait()
  lines = proc.stdout.readlines()
  fields = lines[-1].split()
  percent = 0
  free = '0'
  if fields[-2].endswith('%'):
    try:
      percent = 100-int(fields[-2][:-1])
      free = fields[-3]
    except:
      pass
  return percent,free

backup_runames = set()
while True:
  start = time.time()
  finished = set()
  finished_analysis_time = {}
  finished_analysis = set()
  not_finished = set()
  send_report = False

  text = "Status check on run folders"
  text2=""
  ### CHECK STATUS OF RUN FOLDERS
  for elem in os.listdir(data_folder):
    if os.path.isdir(data_folder+elem) and os.path.isdir(data_folder+elem+"/Data/") and os.path.isdir(data_folder+elem+"/Thumbnail_Images/"): # HAVE A RUN FOLDER with Thumbnail Images and Data folder...
      if (os.path.exists(data_folder+elem+"/Run.completed") and (os.path.exists(data_folder+elem+"/Basecalling_Netcopy_complete_READ2.txt") or os.path.exists(data_folder+elem+"/Basecalling_Netcopy_complete_SINGLEREAD.txt"))) or os.path.exists(data_folder+elem+"/RTAComplete.txt"):
        finished.add(elem)
        if not os.path.exists(data_folder+elem+"/Run.completed"):
          outfile = open(data_folder+elem+"/Run.completed",'w')
          outfile.close()

      runfolder_data = os.listdir(data_folder+elem+"/Data/")
      for dataelem in runfolder_data:
        if os.path.isdir(data_folder+elem+"/Data/"+dataelem):
          firecrest_data = os.listdir(data_folder+elem+"/Data/"+dataelem)
          for firecrestelem in firecrest_data:
            if firecrestelem.startswith('Bustard') or firecrestelem.startswith('Ibis'):
              time_analyzed = None
              if firecrestelem.startswith('Bustard'):
                time_analyzed = firecrestelem.split("_")[-2].replace("-",".")
              elif firecrestelem.startswith('Ibis'):
                time_analyzed = firecrestelem.split("_")[-1].replace("-",".")
              try:
                time_analyzed = (start-time.mktime(time.strptime(time_analyzed,"%d.%m.%Y")))/60/60/24
                finished_analysis_time[elem]=time_analyzed
                finished_analysis.add(elem)
              except:
                text2 += "%s:Problem identifying time analyzed, giving up for this iteration\n"%elem
                send_report = True
              break
              break

      if (elem in finished_analysis) and (elem not in finished) and (not os.path.exists(data_folder+elem+"/NoAutoFinish.txt")):
        text2 += "%s: -> Fixing finished flag\n"%elem
        send_report = True
        if not options.mock:
          if not os.path.exists(data_folder+elem+"/Run.completed"): 
            outfile = open(data_folder+elem+"/Run.completed",'w')
            outfile.close()
          if not os.path.exists(data_folder+elem+"/Basecalling_Netcopy_complete_READ2.txt") and not os.path.exists(data_folder+elem+"/Basecalling_Netcopy_complete_SINGLEREAD.txt") or os.path.exists(data_folder+elem+"/RTAComplete.txt"):
            outfile = open(data_folder+elem+"/RTAComplete.txt",'w')
            outfile.close()
          finished.add(elem)

      if (elem not in finished):
        not_finished.add(elem)

    elif os.path.isdir(data_folder+elem) and (os.path.isdir(data_folder+elem+"/Data/") or os.path.isdir(data_folder+elem+"/Images/")): # HAVE A RUN FOLDER with Images or Data folder...
      if "SYBR" in elem.upper(): continue ## ANOTHER SCRIPT SHOULD TAKE CARE OF THESE FOLDERS

      if os.path.exists(data_folder+elem+"/Run.completed") and (os.path.exists(data_folder+elem+"/Basecalling_Netcopy_complete_READ2.txt") or os.path.exists(data_folder+elem+"/Basecalling_Netcopy_complete_SINGLEREAD.txt") or os.path.exists(data_folder+elem+"/RTAComplete.txt") or os.path.exists(data_folder+elem+"/IPAR_Netcopy_Complete.txt") or os.path.exists(data_folder+elem+"/GA_Netcopy_Complete.txt")):
        finished.add(elem)

      runfolder_data = os.listdir(data_folder+elem+"/Data/")
      for dataelem in runfolder_data:
        if os.path.isdir(data_folder+elem+"/Data/"+dataelem):
          firecrest_data = os.listdir(data_folder+elem+"/Data/"+dataelem)
          for firecrestelem in firecrest_data:
            if firecrestelem.startswith('Bustard') or firecrestelem.startswith('SVM') or firecrestelem.startswith('Ibis'):
              time_analyzed = None
              if firecrestelem.startswith('Bustard'):
                time_analyzed = firecrestelem.split("_")[-2].replace("-",".")
              elif firecrestelem.startswith('Ibis') or firecrestelem.startswith('SVM'):
                time_analyzed = firecrestelem.split("_")[-1].replace("-",".")
              try:
                time_analyzed = (start-time.mktime(time.strptime(time_analyzed,"%d.%m.%Y")))/60/60/24
                finished_analysis_time[elem]=time_analyzed
                finished_analysis.add(elem)
              except:
                text2 += "%s:Problem identifying time analyzed, giving up for this iteration\n"%elem
                send_report = True
              break
              break

      if (elem in finished_analysis) and (elem not in finished) and (not os.path.exists(data_folder+elem+"/NoAutoFinish.txt")):
        text2 += "%s: -> Fixing finished flag\n"%elem
        send_report = True
        if not options.mock:
          if not os.path.exists(data_folder+elem+"/Run.completed"): 
            outfile = open(data_folder+elem+"/Run.completed",'w')
            outfile.close()
          if not os.path.exists(data_folder+elem+"/GA_Netcopy_Complete.txt") and not os.path.exists(data_folder+elem+"/Image_Netcopy_complete.txt") and not os.path.exists(data_folder+elem+"/RTAComplete.txt") and not os.path.exists(data_folder+elem+"/Basecalling_Netcopy_complete_READ2.txt") and not os.path.exists(data_folder+elem+"/Basecalling_Netcopy_complete_SINGLEREAD.txt"):
            outfile = open(data_folder+elem+"/RTAComplete.txt",'w')
            outfile.close()
          finished.add(elem)

      if (elem not in finished):
        not_finished.add(elem)
  if text2 != "": text+="\n"+text2
  text += " ... finished\n"

  ### CHECK STATUS OF BACKUP - SHOULD ONE HAVE A SMALL DATABASE IN CASE THE BACKUP SYSTEM GETS REPLACED?
  text += "\nChecking backup info file"
  text2=""
  infile = open(backup_folder+backup_info_file)
  valid_entry = False ## REQUIRE ONE VALID BACKUP ENTRY BEFORE RESET OF BACKUP SET
  for line in infile:
    fields = line.split()
    if len(fields) == 7: # OLD FILE FORMAT
      filename = fields[-1].split("/")[-1]
      size = fields[0]
      path = "/".join(fields[-1].split("/")[:-1])
      if size == "0": 
        print "Skipping zero entry: %s"%(filename)
        continue
      if os.path.isdir(path) and not filename.startswith(backup_info_file[:6]):
        if not valid_entry: # FIRST VALID ENTRY?
          valid_entry = True
          backup_runames = set()
        if os.path.exists(backup_folder+filename):
          text2 += "%s: File with finished backup\n"%filename
          send_report = True
          if not options.mock:
            try:
              os.remove(backup_folder+filename)
              text2 += "%s: Removed\n"%filename
            except:
              text2 += "%s: Removing failed\n"%filename
        runame = filename.split(".")[0]
        backup_runames.add(runame)
    elif len(fields) == 9: # NEW FILE FORMAT
      filename = fields[4].split("/")[-1]
      size = fields[0]
      path = "/".join(fields[4].split("/")[:-1])
      if size == "0": 
        print "Skipping zero entry: %s"%(filename)
        continue
      if os.path.isdir(path) and not filename.startswith(backup_info_file[:6]):
        if not valid_entry: # FIRST VALID ENTRY?
          valid_entry = True
          backup_runames = set()
        if os.path.exists(backup_folder+filename):
          text2 += "%s: File with finished backup\n"%filename
          send_report = True
          if not options.mock:
            try:
              os.remove(backup_folder+filename)
              text2 += "%s: Removed\n"%filename
            except:
              text2 += "%s: Removing failed\n"%filename
        runame = filename.split(".")[0]
        backup_runames.add(runame)
  infile.close()
  if text2 != "": text+="\n"+text2
  text += " ... finished evaluating %d run names in backup\n"%len(backup_runames)

  ### PRINT FOR RUNS NOT YET IN BACKUP THE DETAILED STATUS
  text += "\nChecking runs finished but not found in backup file"
  text2=""
  backup_todo = set()
  for elem in finished-backup_runames:
    if os.path.exists(backup_folder+elem+".tar.gz") or os.path.exists(backup_folder+elem+".tgz") or os.path.exists(backup_folder+elem+".tar.bz2"):
      text2 += "%s: Archive file is waiting for backup to finish\n"%elem
    else:
      # CHECK whether we are already creating a backup
      if os.path.exists(data_folder+elem+".tar.gz") or os.path.exists(data_folder+elem+".tgz") or os.path.exists(data_folder+elem+".tar.bz2"):
        text2 += "%s: Someone is creating backup for this run in data folder\n"%elem
      elif os.path.exists(data_folder+"/tmp/"+elem+".tar.gz") or os.path.exists(data_folder+"/tmp/"+elem+".tgz") or os.path.exists(data_folder+"/tmp/"+elem+".tar.bz2"):
        text2 += "%s: Someone is creating backup for this run in tmp folder\n"%elem
      elif os.path.exists(backup_folder+elem+".tar.gz.tmp") or os.path.exists(backup_folder+elem+".tar.bz2.tmp"):
        text2 += "%s: Someone is creating backup tmp for this run\n"%elem
      elif os.path.exists(backup_folder+elem+".tgz.tmp"):
        text2 += "%s: Already creating backup for this run\n"%elem
        try:
          filestats = os.stat(backup_folder+elem+".tgz.tmp")
          if ((start-filestats.st_mtime)/60) > timeout_backup_file:
            text2 += "%s: Considering backup process to be dead, restarting the backup\n"%elem
            send_report = True
            try:
              os.remove(backup_folder+elem+".tgz.tmp")
              backup_todo.add(elem)
            except:
              text2 += "%s: Deleting old backup file failed\n"%elem
        except:
          text2 += "%s: Error reading time stamp of archive file\n"%elem
          send_report = True
      else:
        backup_todo.add(elem)
  if text2 != "": text+="\n"+text2
  text += " ... finished\n"

  ## START BACKUP FOR RUNS NOT YET IN BACKUP OR UNDER WAY
  text += "\nChecking for folders waiting to be archived"
  text2=""
  for elem in backup_todo:
    text2 += "%s: Creating archive\n"%elem
    send_report = True
    if not options.mock:
      if os.path.isdir(data_folder+elem+"/Images"):
        compress = subprocess.Popen("cd "+data_folder+"; tar c --exclude-tag-all=config.xml "+elem+" | ahazip > "+backup_folder+elem+".tgz.tmp; mv "+backup_folder+elem+".tgz.tmp "+backup_folder+elem+".tgz",shell=True)
      else:
        compress = subprocess.Popen("cd "+data_folder+"; tar c --exclude-tag-all=Makefile --exclude-tag-all=s_4_finished.txt "+elem+" | ahazip > "+backup_folder+elem+".tgz.tmp; mv "+backup_folder+elem+".tgz.tmp "+backup_folder+elem+".tgz",shell=True)
  if text2 != "": text+="\n"+text2
  text += " ... finished\n"

  text += "\nCleaning old runs"
  text2=""
  ### FOR RUNS IN OR CLOSE TO A FINISHED BACKUP REMOVE IMAGES AND MOVE TO OLD RUNS WITH OLD ENOUGH
  saved_runs = set(finished_analysis.intersection(backup_runames))
  for elem in finished_analysis-backup_runames:
    if os.path.exists(backup_folder+elem+".tar.gz") or os.path.exists(backup_folder+elem+".tgz") or os.path.exists(backup_folder+elem+".tar.bz2"):
      saved_runs.add(elem)

  for elem in saved_runs:
    if (elem in finished_analysis_time) and (finished_analysis_time[elem] >= image_day_wait) and (os.path.isdir(data_folder+elem+"/Images")) and (not os.path.exists(data_folder+elem+"/NoImageDelete.txt")):
      text2 += "%s: Removing images\n"%elem
      send_report = True
      if not options.mock:
        delete = subprocess.Popen("/bin/rm -fR "+data_folder+elem+"/Images",shell=True)
        delete.wait()
    if (elem in finished_analysis_time) and (finished_analysis_time[elem] >= image_day_wait) and os.path.isdir(data_folder+elem+"/Thumbnail_Images/") and (not os.path.exists(data_folder+elem+"/NoImageDelete.txt")):
      text2 += "%s: Removing thumbnail images\n"%elem
      send_report = True
      if not options.mock:
        delete = subprocess.Popen("/bin/rm -fR "+data_folder+elem+"/Thumbnail_Images",shell=True)
        delete.wait()
    if (elem in finished_analysis_time) and (finished_analysis_time[elem] >= old_day_wait) and (not os.path.exists(data_folder+elem+"/NoImageDelete.txt")):
      text2 += "%s: Moving run to old runs folder\n"%elem
      send_report = True
      if not options.mock:
        move = subprocess.Popen("mv "+data_folder+elem+" "+old_folder,shell=True)
        move.wait()
  if text2 != "": text+="\n"+text2
  text += " ... finished\n"

  percent, free = get_diskfree(data_folder)
  text += "\nThe free disk space is: %s %d%%\n"%(free,percent)

  if send_report:
    text = "THIS IS AN AUTOMATICALLY CREATED MESSAGE PLEASE DO NOT REPLY\n\n"+text
    msg = MIMEText(text)
    msg['Subject'] = 'Sequencing backup script report'
    msg['From'] = sender
    msg['To'] = recipient
    if options.mock:
      print msg.as_string()
    else:
      print "Sending out email regarding backup."
      s = smtplib.SMTP('smtp.eva.mpg.de')
      s.set_debuglevel(1)
      s.sendmail(sender, recipient.split(','), msg.as_string())
      s.quit()

  end = time.time()
  print "Going to sleep till next iteration in",max(0,time_wait-(end-start))/60.0,"mins"
  if options.mock: sys.exit()
  time.sleep(max(0,time_wait-(end-start)))


