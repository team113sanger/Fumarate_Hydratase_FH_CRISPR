import ruffus.cmdline as cmdline
from ruffus import *
import json
import re
import shutil
from subprocess import *
import drmaa
import pysam
import tabix
import sys
import os
from distutils.util import strtobool

from drmaa_wrapper_vvi import run_job, error_drmaa_job

#
import json
# Read in a json config file with --config <filename.json> 
#
parser = cmdline.get_argparse(description="fetch seqs for crispr stats")
parser.add_argument('--seq_file_info',type=str,help="flat file containing sample / file information for seq fetch",required=True)
parser.add_argument('--config',type=str,help="json config file containing locations",required=True)
args = parser.parse_args()

mylogger, logger_mutex = cmdline.setup_logging (__name__, args.log_file, args.verbose)

#start the drmaa business
drmaa_session = drmaa.Session()
drmaa_session.initialize()

# read the config file into an object
configfile = open(args.config,'r')
config = json.load(configfile)

# args: cgp_dir, sample_list, output_dir
test_locally = strtobool(config["test_locally"])

output_dir = config["output_dir"]
stdout_dir = config["stdout_dir"]

lsf_8M_option = "-P team113 -M 8000 -R 'select[mem>8000] rusage[mem=8000] span[hosts=1]'" 

cram_dir = output_dir + "/cram"
bam_dir = output_dir + "/bam"
fastq_dir = output_dir + "/fastq"
sample_fastq_dir = output_dir + "/sample_fastq"
stats_dir = output_dir + "/stats/"
job_output_dir = output_dir + "/job_stdout/"

cmd_output = check_output(['mkdir','-p',job_output_dir])
cmd_output = check_output(['mkdir','-p',cram_dir])
cmd_output = check_output(['mkdir','-p',bam_dir])
cmd_output = check_output(['mkdir','-p',fastq_dir])
cmd_output = check_output(['mkdir','-p',sample_fastq_dir])

# 
# Split to get an indicator first - it's easier to parallelise the pass-through 
# than this step
@split(None, cram_dir + "/*.get_cram")
def get_cram_indicators(input_file, output_file):
	seq_files = open(args.seq_file_info, 'r')
	for line in seq_files:
		line = line.rstrip()
		bits = line.split("\t")
		# sample = bits[5] 
		sample = bits[5] + "_" + bits[7] 
		#sample = bits[7] 
		run = bits[12]
		lane = bits[13]
		tag = bits[14]
# Check tehse positions carefully, they are changing
#		run = bits[9]
#		lane = bits[10]
#		tag = bits[11]
		indicator_file = cram_dir + "/" + sample + "__" + run + "_" + lane + "#" + tag + ".get_cram" 
		print(indicator_file)
		oo = open(indicator_file,"w") 

@transform( get_cram_indicators, suffix(".get_cram"), ".cram") 
def get_cram(input_file, output_file):
        # The file we fetch does not have the Sample appended and DOES have the full irods path appended
        fetch_pattern = re.compile("(.*)__(\d+)_(\d+#\d+).get_cram$")
        match = fetch_pattern.match(input_file) 
	print("Input file "+input_file)
	output_cram_file = match.group(1) + "__" + match.group(2) + "_" + match.group(3) + ".cram" 
        irods_fetch = "iget -f /seq/" + match.group(2) + "/" + match.group(2) + "_" + match.group(3) + ".cram  " + output_cram_file
	print(irods_fetch)

	job_other_options = lsf_8M_option 

	run_job_and_move_to_output_dir( cmd=irods_fetch, other_options=lsf_8M_option, step_name="get_cram", 
			named_input_file=input_file, named_output_file=output_file, thelogger=mylogger, thedrmaa_session=drmaa_session)

@transform( get_cram, formatter(), "{subpath[0][1]}/bam/{basename[0]}.bam") 
def make_bam(input_file, output_file):
	print(input_file + ":" + output_file)
	cramtobamcmd = "/software/vertres/bin-external/samtools-1.2/samtools view -b -o " + output_file + " " + input_file 
	
	job_other_options = lsf_8M_option 

	run_job_and_move_to_output_dir( cmd=cramtobamcmd, other_options=lsf_8M_option, step_name="make_bam", 
			named_input_file=input_file, named_output_file=output_file, thelogger=mylogger, thedrmaa_session=drmaa_session)
	
		
@transform( make_bam, formatter(), "{subpath[0][1]}/fastq/{basename[0]}.fastq" ) 
def make_fastq(input_file, output_file):
	print(input_file + ":" + output_file[0])
	# /software/hpag/bin/bamtofastq F=BB_D14_R2_v2__23290_2_2.R1.fastq F2=BB_D14_R2_v2__23290_2_2.R2.fastq filename=BB_D14_R2_v2__23290_2_2.bam
	bamtofastqcmd= "/software/hpag/bin/bamtofastq S="+output_file+" filename="+input_file	
	print(bamtofastqcmd)

	job_other_options = lsf_8M_option 

	run_job_and_move_to_output_dir( cmd=bamtofastqcmd, other_options=lsf_8M_option, step_name="make_fastq", 
			named_input_file=input_file, named_output_file=output_file[0]+":"+output_file[1], thelogger=mylogger, thedrmaa_session=drmaa_session)

@collate(make_fastq,formatter(".+/(?P<BASE>.+)__(.+)\.fastq"),"{subpath[0][1]}/sample_fastq/{BASE[0]}.fastq")
def collate_sample_fastq(input_file, output_file):
	r1_files = list() 
	input_files = input_file
	print("input files "+str(input_files))
	for thefile in input_files:
		print("input file " + thefile)
		file_pattern = re.compile("^(.*).fastq$") 
		match = file_pattern.match(thefile)
		if(match):
			r1_files.append(thefile)

	r1_file_string = " ".join(r1_files)
	output_name1 = output_file 
	#print("detangled output " + output_name1 + " ,,, " + output_name2)
	cat1_command = "cat " + r1_file_string + " > " + output_name1 
	run_job_and_move_to_output_dir( cmd=cat1_command, other_options=lsf_8M_option, step_name="make_sample_fastq", 
			named_input_file="-".join(r1_file_string), named_output_file=output_name1, thelogger=mylogger, thedrmaa_session=drmaa_session)
	print(cat1_command)

def run_job_and_move_to_output_dir(cmd, other_options, step_name, named_input_file, named_output_file, thelogger, thedrmaa_session):
	if(not test_locally):
		stdout_res, stderr_res = run_job(cmd_str = cmd, logger = thelogger, drmaa_session = thedrmaa_session, job_other_options = other_options)
	else:
		stdout_res, stderr_res = run_job(run_locally = test_locally, cmd_str = cmd, logger = thelogger, drmaa_session = thedrmaa_session, job_other_options = other_options)

	move_to_output_dir(step_name, named_input_file, named_output_file, stdout_res, "stdout")
	move_to_output_dir(step_name, named_input_file, named_output_file, stderr_res, "stderr")


def move_to_output_dir(step_name, input_file, output_file, file_to_move, suffix):
	file_leaf_pattern = re.compile(".*/(.*)$") 
	match = file_leaf_pattern.match(input_file)
	if(match):
		input_leaf = match.group(1)
	else:
		input_leaf = input_file


	match = file_leaf_pattern.match(output_file)
	if(match):
		output_leaf = match.group(1)
	else:
		output_leaf = output_file

	new_output_name = stdout_dir + "/" + step_name + "__" + input_leaf + "__" + output_leaf + "."+suffix
	if(not test_locally):
		move_out = check_output(['mv', file_to_move, new_output_name])
	else:
		oo = open(new_output_name,"w")
		for line in file_to_move:
			oo.write(line)

# pipeline_run(verbose = 4, multiprocess=10) multiprocess if running in fg, multithread if running with drmaa
pipeline_run(forcedtorun_tasks=[collate_sample_fastq], verbose = 4, multithread=170, logger = mylogger, multiprocess=10)
drmaa_session.exit()
