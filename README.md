MultiDCoX - Manual

Example Data:

Please refer to an "example" folder data for a list of files and its format required for the input to MDCoX program.

Remarks

All data files’ contents must be tab-separatedThe gene expression file’s first column should be a gene symbols (name) or probe Id. We do not encourage empty gene name or symbol neither the gene name with NA label.  If the user gives the gene name or symbol with NA we will exclude these rows for further processing by MDCoX.

Steps To Run MDCoX:

1)  Re-organize the input files (gene expression data and attribute file) such that samples are organized based on the group or stratum.

a) Input:  Attribute file  Gene expression data file  outputFolder/

b) Command:  Rscript  preProcessAtrributeAndExpression example/attributes.txt example/gene_expressionData.txt preprocessed/

c) Output:  processedExpressionMDCX.txt  processedSampleAttributeMDCX.txt  processedGroupsMDCX.txt  GroupsInfoMDCX.txt
    Remarks:  The attribute file must contain header of GeneName (1st column) and related sample IDs. The GeneName column could be the value of gene name or probe ID. The GroupsInfoMDCX.txt file contains all the groups (stratum-wise) information, the group number and related samples’ columns number of the processedExpressionMDC.txt file. You can to refer to this file for group (stratum-wise) information.

2)  Run the MultiDCoX algorithm

a) Input: processedSampleAttributeMDCX.txt processedExpressionMDCX.txt processedGroupsMDCX.txt outputFolder/
Remarks: The input files for MDCoX program is the output files from pre-processing in Step 1.

b) Command: Rscript MDCoX.R processedSampleAttributeMDCX.txt processedExpressionMDCX.txt processedGroupsMDCX.txt results/
Remarks:  --help For help how to use it, type: Rscript MDCoX.R

c) Output:  MultiDCoxGenesetResults_15_29_59_07_10_2016.txt  Threshold_Density_Plot.png

3)  To retrieve a particular geneset statistic (the related covariates’ model statistic) and its factor-wise or group-wise visualization plots.

a) Input:  A_gene_list_file processedSampleAttributeMDCX.txt processedExpressionMDCX.txt processedGroupsMDCX.txt outputFolder/
Remarks:  For a gene list file input, this can be a list of gene name or probe Id or genes’ row number based on the output results from MDCoX program. If you have duplicate gene name, we encourage you to supply a list of row number output (GIDX) from the MDCoX.R results. Please refer to example/genelist folder for the format (i.e genelist.txt file or genelistRowNumber.txt file).

b) Command: Rscript MGeneset.R genelistRowNumber processedSampleAttributeMDCX.txt processedExpressionMDCX.txt processedGroupsMDCX.txt results/
Remarks:  --help For help how to use it, type:   Rscript MGeneset.R --help

c) Output:  GenesetResults_10_52_39_09_10_2016.txt  COVARIATES_Plot_10_52_39_09_10_2016.png  STRATUM_Plot_10_52_39_09_10_2016.png



Please contact the author of this script Herty Liany, email: e0146315@u.nus.edu if you encountered any problems in running this script or to report the errors. Thank you.
#This program is Free for Non-Commercial Use.#Copyright Year(2017)
