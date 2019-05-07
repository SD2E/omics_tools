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

### Deployment

The visualization aspect can be run on jupyter or the dataframe exported to the clustergrammer2
web application. 

If you want to use this tool to perform differential expression testing then there are some options:
1. Run the tests within python using the rpy2 R interface, best when run on a multi-core machine. 
   
2. If you have access to a HPC cluster an option to produce script files is available. This is 
   probably the fastest compute option, but requires knowledge of a job scheduler and external
   file manipulation. 
   
   Modify ```differential_expression.make_hpc_de_files``` with appropriate system paths.
   
### Versioning

Git

### Authors

* Alexander Cristofaro

### License

This project is licensed under the MIT License - see the LICENSE.md file for details

### Acknowledgements

- Clustergrammer2 for the heatmap visualization.
- GOATOOLS for gene ontology annotation.
- ClusterProfiler for KEGG annotation.   
- edgeR for differential expression.
