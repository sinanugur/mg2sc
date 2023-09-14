#!/usr/bin/env python
import os
import logging
import argparse
import subprocess
from k2sc import mg2sc

############# Import and set up #############

logging.basicConfig(level=logging.DEBUG)
logging.info("Started the run.")

# Define the CLI defaults:
defaults = {'verbosity': 'debug', 'threads': 2, 'kraken': 'kraken2','confidence':0.05,'minimum-hit-groups':4,'complexity':0}

# Create the command line interface:
parser = argparse.ArgumentParser(description='Assign metagenomic assignment to single cells')

# required arguments
parser.add_argument('-i', '--input', dest = 'bamfile', help = "Input bam file", required = True)
parser.add_argument('-o', '--outdir', dest = 'outdir', help = "Directory for output", required = True)
parser.add_argument('-db', '--DBpath', dest ='dbfile', help ="Path to kraken database", required = True)

# optional arguments
parser.add_argument('-v', '--verbose', dest = 'verbosity', default = defaults['verbosity'], choices = ['error', 'warning', 'info', 'debug'], help = 'Set logging level (default {verbosity})'.format(**defaults))
parser.add_argument('-n', '--threads', dest = 'threads', default = defaults['threads'], help = "n cores")
parser.add_argument('-k', '--kraken', dest = 'kraken', default = defaults['kraken'], help = "path to kraken2 executable, if not in $PATH")
parser.add_argument('-prefix', '--prefix', dest = 'prefix', help = "Prefix file output filenames")
parser.add_argument('-confidence', '--confidence', dest = 'confidence',default=defaults['confidence'], help = "Confidence")
parser.add_argument('-minimum-hit-groups', '--minimum-hit-groups', dest = 'hitgroups',default=defaults['minimum-hit-groups'], help = "Minimum hit groups")
parser.add_argument('-complexity', '--complexity', dest = 'complexity',default=defaults['complexity'], help = "Complexity treshold for fastp")
parser.add_argument('-bdb', '--bowtie', dest = 'bowtie', help = "path to a bowtie2 db file to filter out host reads")
#parser.add_argument('-classified-out', '--classified-out', dest = 'classified', help = "Classified file name")


# Parse the CLI:
args = parser.parse_args()
args.threads = str(args.threads)

# Set up logging based on the verbosity level set by the command line arguments:
logging.basicConfig(format='%(levelname)s: %(message)s', level=args.verbosity.upper())
    
# extract name prefix from input bam
if args.prefix is not None:
    prefix = args.prefix
else:  
    prefix = os.path.splitext(args.bamfile)[0]

# Generate variables based on input
bamfile_out = os.path.join(args.outdir,prefix + "_unmapped.bam")
fqfile = os.path.join(args.outdir,prefix + "_unmapped.fq")
krakenoutfile = os.path.join(args.outdir,prefix + "_output.kraken")
reportf = os.path.join(args.outdir,prefix + "_krakenreport.txt")
classout = os.path.join(args.outdir,prefix + "_classified_sequences.txt")

# Make output directories and check that all files exist
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

# Give user feedback about the run
logging.info("Input bamfile: {}".format(args.bamfile))
logging.info("Output krakenfile: {}".format(krakenoutfile))
logging.info("Threads used: {}".format(args.threads))
    
############# Prepare files #############

# Extract unmapped reads from bam
# samtools view -b -f 4 starsoloinput.bam > output_unmapped.bam
cmd1 = "samtools view -@ " + args.threads + " -b -f 4 " + args.bamfile + " > " + bamfile_out
cmd11 = "samtools view -H " + bamfile_out  + " > " + bamfile_out + ".header"
proc1 = subprocess.Popen(cmd1, shell=True)
proc1.wait()
proc11 = subprocess.Popen(cmd11, shell=True)
proc11.wait()
logging.info("Unmapped reads were extracted and saved to {}".format(bamfile_out))

# Convert to fastq
# bedtools bamtofastq -i output_unmapped.bam -fq output_unmapped.fq
#cmd2 = "samtools fastq -@ " + args.threads + " -n " + bamfile_out + " > " + fqfile

cmd2 = "samtools fastq -@ " + args.threads + " -n " + bamfile_out + " | fastp --stdin --stdout --complexity_threshold " + args.complexity + " -Q -L -A -y " + " > " + fqfile #add fastp low complexity filtering


proc2 = subprocess.Popen(cmd2, shell=True)
proc2.wait()
logging.info("FASTQ generated and saved to {}".format(fqfile))



############# Run metagenomic tool on fastq and report summary #############
# run kraken
# kraken2 --threads 2 --db dbpath --report reportfilepath inputfastqpath > krakenoutputfilepath
cmd3 = args.kraken + " --report-minimizer-data " + " --threads " + args.threads + " --confidence " + str(args.confidence) + " --classified-out " + classout + " --minimum-hit-groups " + str(args.hitgroups) + " --db " + args.dbfile + " --report " + reportf + " " + fqfile + " > " + krakenoutfile
proc3 = subprocess.Popen(cmd3, shell=True)
proc3.wait()
logging.info("Kraken2 finished running, krakenoutput saved to {}".format(krakenoutfile))

cmd4= "cat " + fqfile +  " | grep ^@ | tr -d ^@  > " + fqfile + ".id"
proc4 = subprocess.Popen(cmd4, shell=True)
proc4.wait()


if args.bowtie is not None:
    
    ppp=subprocess.Popen("bowtie2 --no-unal --very-fast -x " + args.bowtie +  " -p " + args.threads + " -b  " + bamfile_out + " | cut -f1 | grep -v ^@ > " + fqfile + ".exclude", shell=True)
    ppp.wait()
    ppp=subprocess.Popen("bowtie2 --no-unal --align-paired-reads --very-fast -x " + args.bowtie +  " -p " + args.threads + " -b  " + bamfile_out + " | cut -f1 | grep -v ^@ >> " + fqfile + ".exclude", shell=True)
    ppp.wait()
    ppp=subprocess.Popen("cat " + fqfile + ".id | grep -F -w -v -f " + fqfile + ".exclude | sponge " + fqfile + ".id" , shell=True)
    ppp.wait()



cmd5 = "samtools view -@ " + args.threads + " " + bamfile_out + "  | fgrep -w -f " + fqfile + ".id | " + "awk 'NR==FNR { l[$1]=$0; next } $1 in l {print l[$1]}' - " + fqfile + ".id | samtools view -bS | samtools reheader " + bamfile_out + ".header" +  " - > " + bamfile_out + ".tmp"

proc5 = subprocess.Popen(cmd5, shell=True)
proc5.wait()


cmd6 = "cat " + krakenoutfile + "  | fgrep -w -f " + fqfile + ".id | " + "awk 'NR==FNR { l[$1]=$0; next } $1 in l {print l[$1]}' - " + fqfile + ".id | sponge " + krakenoutfile + ".tmp"

proc6 = subprocess.Popen(cmd6, shell=True)
proc6.wait()


# Split output into single cell level and report sparse matrix (mg2sc.py)
mg2sc(bamfile_out + ".tmp", krakenoutfile + ".tmp", args.dbfile, args.outdir)
logging.info("Sparse matrix with single cell information created")
logging.info("Run finished.")