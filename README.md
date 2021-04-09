# Omics dashboard

Perform differential expression tests for RNASeq and visualize results with tunable 
control over p-value thresholds.

## Getting started
Start with raw count data from RNASeq and format a data file
with factors and features for all samples.

### Prerequisites

Data files not included in the distribution will be automatically downloaded 
when running GO and KEGG annotation. If the internet is unavailable, download 
ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz and extract to the data/ directory.  

### Installing

pip install -I .  
jupyter nbextension enable --py --sys-prefix widgetsnbextension  
jupyter nbextension enable --py --sys-prefix clustergrammer_widget  

### Execution

To run omics tools at scale, run the following command:

```python omics_tools/run_omics.py --input_counts_file <PATH_TO_COUNTS_FILE> --config_file <PATH_TO_CONFIG_FILE> --output_dir <PATH_TO_OUTPUT_DIR>```


The configs listed are the following:
1. ```input_counts_file```: this is a pointer to the raw counts file on your system
2. ```config_file```: this is a configuration file that describes all the comparisons that need to be made. Keys and values in config MUST line up with columns in your ```input_counts_file```
3. ```output_dir```: Path to your output directory

You can find example config files in ```tests/config```. Here is an example of executing for the Bacillus Inducer 1.0 ER from DARPA's SD2 program.

```python omics_tools/run_omics.py --input_counts_file omics_tools/tests/data/experiment.ginkgo.29422_ReadCountMatrix_preCAD_transposed.csv --config_file omics_tools/tests/config/Bacillus_Inducer_1_0.json --output_dir test_new```


### Deployment

The visualization aspect can be run on jupyter or the dataframe exported to the clustergrammer2
web application. 

If you want to use this tool to perform differential expression testing then there are some options:
1. Run the tests within python using the rpy2 R interface. For 200+ tests, use on a workstation (24+ cpu).
   
2. If you have access to a HPC cluster an option to produce script files is available. This is 
   probably the fastest compute option, but requires knowledge of a job scheduler and a few 
   extra steps. 
   
   Modify ```differential_expression.make_hpc_de_files``` with appropriate system paths.
   
### Versioning

Git

### Authors

* Alexander Cristofaro
* Mohammed Eslami
* George Zheng

### License

This project is licensed under the MIT License - see the LICENSE.md file for details

### Acknowledgements

- Clustergrammer2 for the heatmap visualization.
- GOATOOLS for gene ontology annotation.
- ClusterProfiler for KEGG annotation.   
- edgeR for differential expression.
