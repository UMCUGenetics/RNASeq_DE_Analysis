# RNASeq_DE_Analysis
R-script pipeline for DE analysis of RNAseq data as generated by USEQ

These are the scripts used by the UBEC to perform the basic Differential Expression analysis on RNASeq data. 

The main workbook is the Rmd file, containing the various steps in pre- and post-processing and which is used as the basis to generate the final HTML report.

The DE_helper_functions.R file contains most of the functions used by the Rmd file and is in a seperate file so it might also be used for other projects.

## Preparation

First of all make sure all the necessary libraries are installed, these are also listed in the DE_helper_functions.R file.

### Basic package installation

By default the differential expression analysis contains a the heatmap and PCA plot of the samples and the results of the DE analysis using DESeq2 including the heatmaps, volcano plots, MA-plots and tables of the differntially expressed genes.

Install the required packages for the basis DE analysis:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("BiocParallel")
  BiocManager::install("biomaRt")
  BiocManager::install("DESeq2")
  BiocManager::install("clusterProfiler")
    
  install.packages("tidyverse")	#purrr, ggplpot2, dplyr, tidyr, readr, tibble, stringr, forcats. https://www.tidyverse.org/packages/
  install.packages("RColorBrewer")
  install.packages("pheatmap")
  install.packages("ggplot2")
  install.packages("ggrepel")
  install.packages("DT")
  install.packages("ggpubr")
  install.packages('rmarkdown')
  install.packages('plotly')
```

### Annotation packages

While .Rdata files which can be used in the basic DE analysis are supplied, when not used the following annotation packages are required (depending on the species of interest). They are really required when one wants to run the extended GO-Term and pathway analysis as described next.

```
#depends on target species
#human and mouse, most commmon species
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Hs.eg.db")
#dog
#BiocManager::install("org.Cf.eg.db")
```

### Annotation .Rdata files

These annotation data files can be used so that the supporting packages are not needed, nor the need to download the data on the go from the internet...

** add links to library archives **

### Extended package installation

Currently in development is also the inclusion of GO-Term and Pathway analysis modules. For these GOSeq and Gage are being used, while this is in development and not yet used one can also try to run this part and for that the following libraries are needed.

Install the required packages for the extended (pathway and go-term) DE analysis:

```
  BiocManager::install("gage")
  BiocManager::install("gageData")
  BiocManager::install("pathview")
  
  BiocManager::install("goseq")
  #already installed by gage and/or goseq, included for completeness
  #BiocManager::install("GO.db")
```

### Pre-installed library archives

To make things easier we also supply the library folders containing the packages above in a downloadable archive, by passing the location of this library folder to the script the user should not need to install these themselves. Only packages for Windows and CentOS(7) as installed using R 4.2.2 are supplied but hopefully these are usable for other 4.x versions of R and multiple unix based systems as well. This has not been tested and is not guaranteed though!

** add links to library archives **

## Usage

In order to run the DE analysis one needs a readcount file containing the RAW reads, a metadata file containing the experiment setup and sample information and a config file containing the configuration of the DE analysis.

### Readcounts

While not limited to, the current version is based on the readcounts as generated by the Subread FeatureCounts method... If another method or tool is used, please adjust the script as necessary. Start by looking at the 'readCountsFile' and 'removeColumnsFromCountTable' functions (as defined in the helper functions)

### Metadata

The metadata file should be a very simple, tab-delimited, textfile looking like this (the header and columnn names are required):

```
sampleid    condition
C1          case
C2          case
C3          case
WT1         control
WT2         control
WT3         control
...
```

The values in the sampleid column should match the names of the corresponding columns in the readcount file.
The values in the condition column are the conditions, or groups, to be used in the  comparisons in the DE analyis

### Config

Finally one needs to supply the config.txt file, which will tell the script where to actually find the readcounts, metadata, comparisons to perform and also specify other information...

The name should be config.txt and the location should be the same folder as the main .Rmd (and DE_helper_functions.R) file.

```
library_folder="/path/to/R/custom_lib_folder/"
project_name="PID1234"
metadata_file="/path/to/PID1234_metadata.txt"
readcount_file="/path/to/PID1234_processed_exon_featureCounts.raw.txt"
species="grch38"
annotation_data="/path/to//biomart_annodb_grch38.RData"
comparision="
condition,case,control
"
```

All values should be surrounded by quotation marks, the .Rmd script will simply read in these parameters and set them accordingly.

The project_name, metadata_file and readcount_file parameters hopefully speak for themselves and need no further explanations (apart from not using special character in the project_name value).

The library_folder can be used to point to a location where all packages are installed, for example the locatio of the (extracted) pre-downloaded library archive. When the library_folder is empty or left out, the script will simply assume all necessary packages are installed in the default library folder...

The annotatio_data can be used to point to one of the pre-downloaded annotation .Rdata files, if not it will try to download the necessary annotation data from biomaRt or annotationDb. In that case the species parameter also needs to be suppied.

comparision: This is used to define the comparisons to be done for this analysis. Each comparison is put on a seperate, comma-seperated, line. The first element is 'condition' matching the 'condition' column from the metadata file. Next are the two groups to compare, in the example above case vs control.

When one has multiple groups to compare, they can simply be specified on multiple lines. For example:

```
comparision="
condition,group_1,group_2
condition,group_1,group_3
condition,group_2,group_3
"
```

The script will loop through these automatically and generate the report for all of them in one go.

### Parameters

Some other parameters can also be set in the config.txt file, if not specified the script will use it's default settings:

```
#fdr.cutoff is the false discovery (or p-adjusted) cutoff used for differentially expressed gene detection
fdr.cutoff=0.05 

#fc.cutoff is the fold-change cutoff  used for differentially expressed gene detection
fc.cutoff=1.5 

#group.replicates.by can be used to specify by which column technical replicates should be grouped, to be used when an experiment contains 
#both biological as technical replicates. In that case the technical replicates will first be collpased into a single value, after which the biological
#replicates are used as actual replicates
group.replicates.by=NULL

#a boolean value specifying whether GO-term analysis should be performed (in beta/development, use at your own risk!)
perform.go=FALSE 

#a boolean value specifying whether KEGG-pathway analysis should be performed (in beta/development, use at your own risk!)
perform.kegg=FALSE 

#a boolean value specifying whether the top 3 pathways from the KEGG-pathway analysis should be plotted
plot.kegg=FALSE 

#a boolean value specifying whether the various plots as generated in the report should also be exported to pdf format
export.plots.pdf=TRUE 

#a boolean value specifying whether the various plots as generated in the report should also be exported to svg format
#set to FALSE because filesizes for SVG can be quite large
export.plots.pdf=FALSE 

#a boolean value specifying whether only the samples as defined in the metadatafile should be included for all processing steps (including normalization,e tc)
samples.in.comp=TRUE 
```

## Running

In order to finally run/knit the script one can either open the .Rmd file in R-Studio and choose to knit the document from there, or open an R-Session and execute the rmarkdown::render command listed below.

Pandoc is used for the final knitting, while this is (or should be) already included with RStudio you most likely will need to first install it before trying to knit the .Rmd file to a HTML-document directly in R. Please find the proper way of doing this for your OS.

For example: 

```
#ubuntu
sudo apt-get install pandoc

#centos 
yum install pandoc
```

When pandoc is installed, the HTML report can be generated by typing the following command(s)

```
#.libPaths command is only needed when you are using a custom lib folder for this project
.libPaths("/path/to/custom/libs")

#actually generate the HTML file (replace XXX_DE_Analysis.Rmd with the name of your Rmd file)
rmarkdown::render("XXXX_DE_Analysis.Rmd", "html_document")
```

