import sys
sys.path.insert(0, '../')

import argparse
import os
import pheflux


# Description:
# ------------
parser = argparse.ArgumentParser(description='')

# Arguments:
# ----------
parser.add_argument("-i","--input_file", help="Name of the csv file with required information", type=str)
parser.add_argument("-l","--prefix_log_file", help="Prefix to *_record_XXXX.log.csv output file. Example: Glucose_record_XXXX.log.csv", type=str, default="")
parser.add_argument("-o","--output_directory", help="Name of the directory where output files will be stored", type=str)
parser.add_argument("-v","--verbosity", help="Verbose mode", action="store_true")

args = parser.parse_args()

print("Input file: ", args.input_file, "\n")

# Create an output directory:
# ---------------------------
# output directory:
resultsDir = args.output_directory
if os.path.exists(resultsDir) == False:
	path = os.path.join(resultsDir)
	os.mkdir(path)

# Run Pheflux:
# ------------
inputFile = args.input_file
pheflux.getFluxes(inputFile,resultsDir,args.prefix_log_file,args.verbosity) #(-i, -o, -l, -v)

print("")
print("Output directory:", resultsDir)
