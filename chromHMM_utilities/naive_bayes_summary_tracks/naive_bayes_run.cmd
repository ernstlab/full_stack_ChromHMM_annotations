#!/bin/csh -f
#  naive_bayes_run.cmd
#
#  UGE job for naive_bayes_run built Thu Feb 22 11:21:06 PST 2018
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/project/ernst/havu73/source/naive_bayes/naive_bayes_run.joblog.$JOB_ID
#$ -o /u/project/ernst/havu73/source/naive_bayes/naive_bayes_run.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/project/ernst/havu73/source/naive_bayes/naive_bayes_run
#  arguments       = 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Resources requested
#
#$ -l h_data=16G,h_rt=72:00:00,highp
# #
#  Name of application for log
#$ -v QQAPP=job
#  Email address to notify
#$ -M havu73@g.ucla.edu
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n
#
# Initialization for serial execution
#
  unalias *
  set qqversion = 
  set qqapp     = "job serial"
  set qqidir    = /u/project/ernst/havu73/source/naive_bayes
  set qqjob     = naive_bayes_run
  set qqodir    = /u/project/ernst/havu73/source/naive_bayes
  cd     /u/project/ernst/havu73/source/naive_bayes
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for naive_bayes_run built Thu Feb 22 11:21:06 PST 2018"
  echo ""
  echo "  naive_bayes_run directory:"
  echo "    "/u/project/ernst/havu73/source/naive_bayes
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  SCRATCH directory:"
  echo "    "$qqscratch
#
  echo ""
  echo "naive_bayes_run started on:   "` hostname -s `
  echo "naive_bayes_run started at:   "` date `
  echo ""
#
# Run the user program
#
  source /u/local/Modules/default/init/modules.csh
  module load intel/13.cs
#
  echo naive_bayes_run "" \>\& naive_bayes_run.output.$JOB_ID
  echo ""
  module load pypy
  module load python/2.7

  /usr/bin/time pypy naive_bayes_version2.py ../../naive_bayes/conditional_probabilities.txt ../../naive_bayes/prior_probabilities.txt ../../even_smaller_binary_data_80/ ../../ChromHMM/JASON_MODEL/segments/genome_100_segments.bed  ../../naive_bayes/output_cluster_02252018/trial.txt 80  >& /u/project/ernst/havu73/source/naive_bayes/naive_bayes_run.output.$JOB_ID
#
  echo ""
  echo "naive_bayes_run finished at:  "` date `
#
# Cleanup after serial execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/project/ernst/havu73/source/naive_bayes/naive_bayes_run.joblog.$JOB_ID --------" >> /u/local/apps/queue.logs/job.log.serial
  if (`wc -l /u/project/ernst/havu73/source/naive_bayes/naive_bayes_run.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
	head -50 /u/project/ernst/havu73/source/naive_bayes/naive_bayes_run.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
	echo " "  >> /u/local/apps/queue.logs/job.log.serial
	tail -10 /u/project/ernst/havu73/source/naive_bayes/naive_bayes_run.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  else
	cat /u/project/ernst/havu73/source/naive_bayes/naive_bayes_run.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  endif
  exit (0)
