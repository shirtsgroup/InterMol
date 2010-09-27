#!/usr/bin/python

# Submit jobs to LSF queue, storing jobids in job directories to guard against accidental restarts.

#========================================================================
# IMPORTS
#========================================================================

import os
import commands
import re

#========================================================================
# PARAMETERS
#========================================================================

# job queue
queue = "SP" # serial queue, standard priority

# runtime limit
runtime_limit = "168:00" # one week -- the maximum

#========================================================================
# SUBROUTINES
#========================================================================

def read_file(filename):
   """Read the contents of a file into an array of lines.

   ARGUMENTS
       filename (string) - name of the file to be read

   RETURN VALUES
       lines (list of strings) - the contents of the file
   """

   # Open the file for writing
   file = open(filename, 'r')

   # read the file contents
   lines = file.readlines()

   # Close the file
   file.close()

   return lines

def write_file(filename, contents):
   """Write the given contents to a text file.

   ARGUMENTS
       filename (string) - name of the file to write to, creating if it doesn't exist
       contents (string) - contents of the file to be written
   """

   # Open the file for writing
   file = open(filename, 'w')

   # Write the file contents
   file.write(contents)

   # Close the file
   file.close()

   return

def lsf_submit(directory, jobfile, hold_jobid = None, shell = '/bin/tcsh', jobname = None, queue = None, runtime_limit = None):
  """Submit the given job file to the LSF queue.

  ARGUMENTS
      directory (string) - the directory to change to before submission
      jobfile (string) - the name of the job script to submit

  OPTIONAL ARGUMENTS
      hold_jobid (string or int) - job name or id to hold on as a prerequisite for execution
      jobname (string) - a name for the job
      queue (string) - job queue to submit to

  RETURNS
      jobid (integer) - the jobid
  """

  # Form command
  command = 'bsub'
  command += ' -L %s' % shell
  if jobname: command += ' -J %s' % jobname
  if hold_jobid: command += ' -w "done(%s)"' % hold_jobid 
  if queue: command += ' -q %s' % queue
  if runtime_limit: command += ' -W %s' % runtime_limit
  command += ' -o stdout -e stderr'
  command += ' %s' % jobfile

  # Change to directory
  import os
  current_directory = os.getcwd()
  os.chdir(directory)

  # Submit the job and capture output.
  import commands
  print "> " + command
  output = commands.getoutput(command)
  print output

  # Match job id
  import re
  matches = re.match('Job <(\d+)>', output)
  jobid = matches.group(1)

  # Restore working directory
  os.chdir(current_directory)

  return int(jobid)

#========================================================================
# MAIN
#========================================================================


#========================================================================
# get a list of directories in which jobs are to be submitted
#========================================================================

# read contents of file as newline-terminated lines
lines = read_file('directories')
# strip newlines and compile list of directories
directories = list()
for line in lines:
   directories.append(line.strip())

# read jobids, if file exists
jobids_filename = 'jobids'
if os.path.exists(jobids_filename):
   lines = read_file(jobids_filename)
else:
   lines = list()

# form list
jobids = set()
for line in lines:
   jobid = int(line.strip())
   jobids.add(jobid)

# submit job from each directory
for directory in directories:
   # if there is a 'jobid' file listed, skip if the process is still running
   jobid_filename = os.path.join(directory, 'jobid')
   if os.path.exists(jobid_filename):
      # get jobid
      lines = read_file(jobid_filename)
      jobid = lines[0].strip()
      
      # if job is running, skip this job
      command = 'bjobs ' + jobid
      output = commands.getoutput(command)
      if re.search('RUN', output) or re.search('WAIT', output):
         print 'Job ' + jobid + ' is still running or waiting -- refusing to start job again.'
         print output
         continue
   
   # construct jobname from directory name
   print directory
   jobid = lsf_submit(directory, 'tcsh run.sh', queue = queue, runtime_limit = runtime_limit)

   # write jobid
   write_file(jobid_filename, str(jobid) + '\n')

   # keep track of submitted jobids
   jobids.add(jobid)
   
print "Done."

# write submitted jobids
# form text file from existing job ids
jobids_text = ""
for jobid in jobids:
   jobids_text += "%d\n" % jobid

filename = 'jobids'
write_file(filename, jobids_text)

