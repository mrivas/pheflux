# Examples

The use of ```pheflux.py``` is exemplified in the following files:

- tutorial.ipynb: Example run of ```pheflux.py``` using Jupyter notebook.
- example.sh: Example run of ```pheflux.py``` using Linux bash script. It calls the ```pheflux-terminal.py``` file. 

Both examples use the information stored in the ```data``` folder and save the outputs in the ```results `` folder. The inputs and outputs are described below.

## Inputs

Four inputs are needed:

- ```inputFile```       : Name of the csv file with required information.
- ```resultsDir```      : Name of the directory where output files will be stored.
- ```verbosity```       : Verbose mode.
- ```prefix_log_file``` : Prefix log output file.


### Input file

The input file is a comma-separated file containing: the name of the organism, the condition, medium, gene expression file, and genome-scale metabolic network. 

More than one row can be defined in the same file. Pheflux will be run on the data provided for each line. For example, the input file ```examples/data/InputData_example.csv``` contains the following:

| **Organism** | **Condition**          | **GeneExpFile**                                                                                        | **Medium**                           | **Network**           |
|--------------|------------------------|--------------------------------------------------------------------------------------------------------|--------------------------------------|-----------------------|
| Scerevisiae  | Glucose                | data/transcriptomes/Scerevisiae_Expfile_Glucose.csv                                                    | data/mediums/Scerevisiae_Medium.csv  | data/gems/iMM904.xml  |
| Ecoli        | Glucose                | data/transcriptomes/Ecoli_Expfile_Glucose.csv                                                          | data/mediums/Ecoli_Medium.csv        | data/gems/iJO1366.xml |
| Homo_sapiens | PrimaryTumor_StageIIIA | data/transcriptomes/002aa584-7dd3-4c62-952b-b24e5858482b/a5dc521e-bee4-489c-8679-d4b90a327d33.FPKM.txt | data/mediums/Hsapiens_cellMedium.csv | data/gems/Recon3D.xml |

### Results directory

A string describing the folder where the output files will be stored. 

### Verbosity

Indicates whether a verbose output (```verbosity = True```) will be printed on terminal or not (```verbosity = False```).

### Prefix of Log file

Pheflux also produces a log file with summary statistics of all computations. The user-provided prefix is used to name this file. For example, if ```prefix = example``` an ```example_record_XXXX.log.csv``` file will be created, where ```XXXX``` is a random four-character tag.

## Outputs

Two output files are generated:

- ```Organism_Condition_STATUS.fluxes.csv```: Comma-separated file containing the fluxome estimation. 
- ```Organism_Condition_STATUS.Log.csv```: Comma-separated file containing various statistics.

For both files, the words ```Organism```, and ```Condition``` are extracted from the ```examples/data/InputData_example.csv``` file, and the ```STATUS``` word indicates ```IPOPT``` finalization condition.
