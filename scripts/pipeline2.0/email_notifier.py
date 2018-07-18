#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
Periodically scans for new or finished solexa processing jobs. Sends email notifications.
Runs as daemon.

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *16.11.2008
:Type: daemon
:Input: -
:Output: -

options:
    -t, --target            Path to final storage (default /mnt/solexa)
    -m, --mock              Make a test run

"""

import sys
import os
import time
import subprocess
from optparse import OptionParser
import smtplib
from email.mime.text import MIMEText

time_wait = 15*60
storage = 20
time_storage = 60*60*24
last_storage = time.time()-time_storage
sender = 'sbsuser@bioaps01'
recipient = 'martin_kircher@eva.mpg.de,schultz@eva.mpg.de,sequencing@eva.mpg.de,kelso@eva.mpg.de,udo_stenzel@eva.mpg.de'

def get_diskfree(path):
  proc = subprocess.Popen('df -h '+path,shell=True,stdout=subprocess.PIPE)
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

parser = OptionParser()
parser.add_option("-t", "--target", dest="target", help="Path to final storage (default /mnt/solexa)",default="/mnt/solexa/")
parser.add_option("-m", "--mock", dest="mock", help="Make a test run",default=False,action="store_true")
(options, args) = parser.parse_args()

data_folder = options.target.rstrip('/')+'/'
old_set = set()
miseq_notified = set()

while True:
  start = time.time()
  if ((start - last_storage) > time_storage) or options.mock:
    percent, free = get_diskfree(data_folder)
    if percent <= storage: # or options.mock:
      text = "THIS IS AN AUTOMATICLY CREATED MESSAGE PLEASE DON'T REPLY\n\n"
      text += "The disk space for storing Solexa run data is running low.\n"
      text += "The free disk space is less or equal to "+str(storage)+"%:\n\n"
      text += "     "+free+" ("+str(percent)+"%)\n\n"
      text += "Please consider archiving runs or removing images.\n "
      msg = MIMEText(text)
      msg['Subject'] = 'Disk space running low'
      msg['From'] = sender
      msg['To'] = recipient
      if options.mock:
        print msg.as_string()
        last_storage=start
      else:
        try:
          print "Sending out email regarding disk space."
          s = smtplib.SMTP('smtp.eva.mpg.de')
          s.set_debuglevel(1)
          s.sendmail(sender, recipient.split(','), msg.as_string())
          s.quit()
          last_storage=start
        except:
          print "Error connecting to SMTP server"

  infolder = os.listdir(data_folder)
  text = "THIS IS AN AUTOMATICLY CREATED MESSAGE PLEASE DON'T REPLY\n\n"
  finished = []
  not_finished = []
  notify = False

  for elem in infolder:
    if os.path.isdir(data_folder+elem) and (os.path.isdir(data_folder+elem+"/Data/") or os.path.isdir(data_folder+elem+"/Images/")): # HAVE A RUN FOLDER...
      if not os.path.exists(data_folder+elem+"/RTAComplete.txt") and os.path.exists(data_folder+elem+'/SampleSheet.csv'):
        if elem in miseq_notified: continue
        infile = open(data_folder+elem+'/SampleSheet.csv')
        reads_length = []
        ExperimentName, ExperimenterName, DescriptionText, ProjectName = '','','',''
        Email = None
        ReadsStarted = False
        for line in infile:
          line = line.strip()
          if line.startswith('Investigator Name,'):
            ExperimenterName = line.split(',')[1].strip()
          elif line.startswith('Project Name,'):
            ProjectName = line.split(',')[1].strip()
          elif line.startswith('Experiment Name,'):
            ExperimentName = line.split(',')[1].strip()
          elif line.startswith('Description,'):
            DescriptionText = line.split(',')[1].strip()
            if DescriptionText == "": DescriptionText = "<i>None</i>"
          elif line.startswith('Email'):
            Email = line.split(',')[1].strip().split(':')
          elif line.startswith('[Reads]'):
            ReadsStarted = True
          elif line.isdigit() and ReadsStarted:
            reads_length.append(int(line))
          elif not line.isdigit() and ReadsStarted:
            ReadsStarted = False
          elif line.startswith('Sample_ID'): 
            if "index," in line: reads_length.insert(1,7)
            if "index2," in line: reads_length.append(7)
        infile.close()
        if Email != None:
          cycle_text = "+".join(map(str,reads_length))
          cycle_time = round(sum(reads_length)/10.+2)

          html_text = """<p> ======= THIS IS AN AUTOMATICALLY CREATED MESSAGE ======= </p>
          <p>Dear all,</p>
          <p>We started a new experiment (%s) on %s. The run id is <b>%s</b> and it will finish within about <b>%d hours</b> (%s cycles).</p>
          <p>Please fill out an analysis request for this run before the run finishes using the internal form available <a href="http://biofs05/~sbsuser/index.py/form?runid=%s">here</a>. 
          <i>If this form is not available when the run finishes, your lanes will not be processed!</i> </p>
          <p>Operator/Contact:&nbsp;&nbsp;%s<br>
          Comment:&nbsp;&nbsp;%s</p>
          <p>Library: %s</p>
          <p>If you have any questions feel free to contact us.</p>
          <p>Best regards,<br>
          SeqGroup</p>"""%(ExperimentName,elem.split('_')[1],elem,cycle_time,cycle_text,elem,ExperimenterName,DescriptionText,ProjectName)

          msg = MIMEText("""<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd"><html><body>"""+html_text+"""</body></html>""")
          msg.set_type('text/html')
          msg['Subject'] = '%s started/analysis information request'%elem
          msg['From'] = 'sequencing@eva.mpg.de'
          msg['To'] = ",".join(Email)
          if not options.mock:
            print "Sending out email regarding a new MiSeq experiment."
            try:
              s = smtplib.SMTP('smtp.eva.mpg.de')
              s.set_debuglevel(1)
              s.sendmail('sequencing@eva.mpg.de', Email, msg.as_string())
              s.quit()
              miseq_notified.add(elem)
            except:
              print "Sending message failed... Will try again in next iteration."
          else:
            print msg.as_string()

      if (os.path.exists(data_folder+elem+"/Run.completed") and (os.path.exists(data_folder+elem+"/Basecalling_Netcopy_complete_SINGLEREAD.txt") or os.path.exists(data_folder+elem+"/Basecalling_Netcopy_complete_READ2.txt") or os.path.exists(data_folder+elem+"/IPAR_Netcopy_Complete.txt") or os.path.exists(data_folder+elem+"/GA_Netcopy_Complete.txt"))) or os.path.exists(data_folder+elem+"/RTAComplete.txt"):
        if elem in miseq_notified: miseq_notified.remove(elem)
        if "SYBR" in elem.upper(): continue ## ANOTHER SCRIPT SHOULD TAKE CARE OF THESE FOLDERS

        analyzed = False
        if os.path.isdir(data_folder+elem+"/Data/Intensities/BaseCalls/RTAreport/"): analyzed=True
        else:
          runfolder_data = os.listdir(data_folder+elem+"/Data/")
          for dataelem in runfolder_data:
            if os.path.isdir(data_folder+elem+"/Data/"+dataelem):
              firecrest_data = os.listdir(data_folder+elem+"/Data/"+dataelem)
              for firecrestelem in firecrest_data:
                if firecrestelem.startswith('Bustard') or firecrestelem.startswith('SVM') or firecrestelem.startswith('Ibis'):
                  analyzed=True
                  break
                  break

        if analyzed:
          finished.append(elem+" (analyzed)")
        else:
          if elem not in old_set:
            old_set.add(elem)
            notify=True
            text += "A new run finished and is waiting for analysis\n\n   "+elem+"\n\n"
          finished.append(elem)
      else:
        not_finished.append(elem)
  finished.sort()
  not_finished.sort()
  text += "The following run folders finished sequencing:\n\n"
  text += '\n'.join(finished)
  text += "\n\nThe following run folders seem not finished:\n\n"
  text += '\n'.join(not_finished)
  text += "\n"
  msg = MIMEText(text)
  msg['Subject'] = 'New solexa runs available'
  msg['From'] = sender
  msg['To'] = recipient
  if options.mock:
    print msg.as_string()
  elif notify and not options.mock:
    not_send = True
    while not_send:
      try:
        print "Sending out email regarding a new run."
        s = smtplib.SMTP('smtp.eva.mpg.de')
        s.set_debuglevel(1)
        s.sendmail(sender, recipient.split(','), msg.as_string())
        s.quit()
        not_send = False
      except:
        print "Error connecting to SMTP server"
        end = time.time()
        print "Going to sleep and trying to send email again in",max(0,time_wait-(end-start))/60.0,"mins"
        time.sleep(max(0,time_wait-(end-start)))

  end = time.time()
  print "Going to sleep till next iteration in",max(0,time_wait-(end-start))/60.0,"mins"
  time.sleep(max(0,time_wait-(end-start)))


