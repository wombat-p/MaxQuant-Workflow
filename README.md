# Workflow to analyze label-free data with MaxQuant
This workflow is based on Nextflow, running with SDRF implemented. Normalization and statistical comparisons using NormalyzerDE are conducted on the MaxQuant results.

## Content
_Nextflow folder_: The workflow implementation

_data folder:_ additional data to run the UPS data set (NOT including the RAW files)

_Results folder_: Results for the UPS data set.


## Getting started

The Nextflow script can run workflows from different configuration files, such as
1) with SDRF file (raw files can be given as parameter or are download from the location specified in the sdrf file):
a) SDRF file + fasta file
b) SDRF file + fasta file + experimental design file (will overwrite experimental design in sdrf)
c) SDRF file + fasta file + experimental design file + yaml parameter file (will overwrite default and sdrf parameters)

2) without SDRF file:
a) Raw files + fasta file + yaml parameter file
b) Raw file + fasta file + yaml parameter file + experimental design file

The SDRF can be found under annotated projects, and for the PXD001819, the file is added under data.
URL for SDRF files: https://github.com/bigbio/proteomics-metadata-standard/tree/master/annotated-projects

An experimental design file for the NormalyzerDE part can also be found in the data folder. 

For the specification of the yaml parameter file, see https://github.com/bigbio/proteomics-metadata-standard/blob/master/sdrf-proteomics/Data-analysis-metadata.adoc and the example in the Nextflow/data folder

## Run benchmarking data set

Download the raw files from PRIDE: http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001819

Run the workflow, giving the following parameters:

1) The SDRF.tsv file. Relative paths works fine.
2) The fasta file, and this needs a absolute file path, not relative. 
3) (Optional) The absolute path to the experimental design file 
4) (Optional) Parameter file (in yaml format)
5) (Optional) The group comparisons to perform in NormalyzerDE. Without this parameter, the comparison will be versus 
the first condition

Just make sure to update the paths in the configuration file, and then run as

```
nextflow run main.nf --fasta 'https://raw.githubusercontent.com/wombat-p/MaxQuant-Workflow/dev/data/yeast_UPS.fasta' --sdrf 'https://raw.githubusercontent.com/wombat-p/MaxQuant-Workflow/dev/data/sdrf_UPS.tsv' --max_cpus 8 --max_memory 8GB -profile docker  -with-report -with-trace -with-timeline
```

Still missing: parameter for tolerances, enzymes and PTMs, possibility to run without sdrf file

