import subprocess
import os

# Define the path to your C++ executable
cpp_executable = "./bin/CreateTuples"
base_input_directory = "/home/truga/projects/phys565/project/"
output_directory = "/home/truga/projects/phys565/project/data/proccess_tuples/"
os.makedirs(output_directory, exist_ok=True)

datasets = ["dy", "qcd", "wjets", "single_top", "ww", "wz", "zz", "ttbar", "data"]


# List of input arguments
input_arguments = []

# ./bin/CreateHistograms /eos/uscms/store/user/csanmart/analyzer_HiggsMuMu/MC_background/TTto2L2Nu_Summer22/SumGenWeight.root ZZto4Lhist.root 2022 ZZto4L F
for dataset in datasets:
    input = base_input_directory + "data/tuples/"+ dataset + ".root"
    input_arguments.append([input, output_directory, dataset])

# Loop over each set of input arguments and execute the C++ program
for args in input_arguments:
    # Run the C++ executable with the current arguments
    result = subprocess.run([cpp_executable] + args, capture_output=True, text=True)

    # Print the output and any error messages
    print("Output:", result.stdout)
    print("Errors:", result.stderr)
