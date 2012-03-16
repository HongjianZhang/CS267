#! /usr/bin/python
# Purpose: Collect performance data from the current build
import re
import sys


if(len(sys.argv) < 2):
	print "ERROR: Need build file to analyze";
	exit();

datafile = open(sys.argv[1], 'r');
datalines = datafile.readlines();
datafile.close();

print 'Aggregating performance data in file'

current_location = 0;
# 0 means not yet in any region of interest
# 1 means read output summary line
# 2 means reading IPM summary
# 3 means finished IPM summary
# 4 means reading compute
# 5 means compute %MPI read
# 6 means finished reading compute

found_n = -1;
found_nproc = -1;
found_time = -1;
found_wclock = -1;
found_flops = -1;
found_cclock = -1;
found_pMPI = -1;

print "n".ljust(12) + "nproc".ljust(12) + "time(s)".ljust(12) + "futil%".ljust(12) + "comm%".ljust(12) + "barr%".ljust(12) + "comp%".ljust(12)

for line in datalines:
        if current_location == 0 and line.startswith('n = '):
                found = re.search('n = ([0-9.+-]+), n_procs = ([0-9.+-]+), simulation time = ([0-9.+-]+)', line)
                found_n = int(found.group(1))
                found_nproc = int(found.group(2))
                found_time = float(found.group(3))
                current_location = 1
        if current_location == 1 and line.startswith('# command'):
                current_location = 2
        if current_location == 2 and line.startswith('# stop'):
                found = re.search('# stop[\s\w:]*wallclock\s:\s([0-9.]+)', line)
                found_wclock = float(found.group(1))
        if current_location == 2 and line.startswith('# mem'):
                found = re.search('# mem \[GB][\s\w.:]*gflop/sec\s*:\s*([0-9.]+)',line)
                found_flops = float(found.group(1))
                current_location = 3
        if current_location == 3 and line.startswith('# region    :\'compute\''):
                current_location = 4
        if current_location == 4 and line.startswith('# wallclock'):
                found = re.search('# wallclock :\s*[0-9.]*\s*([0-9.]+)', line)
                found_cclock = float(found.group(1))
        if current_location == 4 and line.startswith('#   MPI     :'):
                found = re.search('#   MPI     :\s*([0-9.]+)', line)
                found_pMPI = float(found.group(1))
                current_location = 5
        if current_location == 5 and line.startswith('# MPI_Barrier'):
                found = re.search('# MPI_Barrier\s*[0-9.]+\s*[0-9.]+\s*([0-9.]+)', line)
                found_pBarrier = float(found.group(1))
        if current_location == 5 and line.startswith('###'):
                current_location = 6
        if current_location == 6 and line.startswith('Application'):
		current_location = 0
		calc_flops = (found_flops/found_nproc)*(found_wclock/found_cclock)
		calc_futil = calc_flops/8.504*100
		calc_comm = found_pMPI - found_pBarrier
		calc_comp = 100 - calc_comm - found_pBarrier
		print str(found_n).ljust(12) + str(found_nproc).ljust(12) + str(found_time).ljust(12) + ("%.6f"%calc_futil).ljust(12) + str(calc_comm).ljust(12) + str(found_pBarrier).ljust(12) + str(calc_comp).ljust(12)
		found_n = -1;
		found_nproc = -1;
		found_time = -1;
		found_wclock = -1;
		found_flops = -1;
		found_cclock = -1;
		found_pMPI = -1;
		found_pBarrier = -1;
