MultiDCoX - Manual
------------------

Example Datasets
-------------------
Please refer to an "example" folder  for a list of files and format required for the input to MDCoX program.

Remarks:
All data files’ contents must be tab-separatedThe gene expression file’s first column should be a gene symbols (name) or probe Id. We do not encourage empty gene name or symbol neither the gene name with NA label.  If the user gives the gene name or symbol with NA we will exclude these rows for further processing by MDCoX.

Simulated Data with Co-expressed Genesets
-----------------------------------------
The simulated data (in "example" folder) contains a few co-expressed genesets under different factors, the true positive genesets with its indices as follows:<br />

a) For ER+,P53+ only:<br />
   4,6,12,14,18,20,29,36,40,44,49,61,72,74,81,84,89,97,99,362

b) For ER- only:<br />
   7,76,180,195,217,328,373,393,400,401,411,432,456,471,480,482,489,490,494,497

Notes: You will expect the simulated results (for different factors) contain the genesets (indices) as above.

Steps To Run MultiDCoX
---------------------
1)  Re-organize the input files (gene expression data and attribute file) such that samples are organized based on the group or stratum.

a) Input: <br />
Attribute file<br />
Gene expression data file <br />
outputFolder/ <br />

b) Command:  Rscript  preProcessAttributeAndGeneExpression.R example/attributes.txt example/gene_expressionData.txt preprocessed/

c) Output:<br />
processedExpressionMDCX.txt <br />
processedSampleAttributeMDCX.txt <br />
processedGroupsMDCX.txt<br />
GroupsInfoMDCX.txt<br />
<br />Remarks:  The attribute file must contain header of GeneName (1st column) and related sample IDs. The GeneName column could be the value of gene name or probe ID. The GroupsInfoMDCX.txt file contains all the groups (stratum-wise) information, the group number and related samples’ columns number from the processedExpressionMDC.txt file.

2)  Run the MultiDCoX algorithm

a) Input:<br />
processedSampleAttributeMDCX.txt<br />
processedExpressionMDCX.txt<br />
processedGroupsMDCX.txt<br />
outputFolder/ <br />
<br />Remarks: The input files for MDCoX program is the output files from pre-processing Step 1.

b) Command: Rscript MDCoX.R processedSampleAttributeMDCX.txt processedExpressionMDCX.txt processedGroupsMDCX.txt results/ 
<br />Remarks:  --help For help how to use it, type: Rscript MDCoX.R --help

c) Output:<br />
MultiDCoxGenesetResults_15_29_59_06_10_xxx.txt<br />
Threshold_Density_Plot.png<br />

3)  To retrieve a particular geneset statistic (the related covariates’ model statistic) and its factor-wise or group-wise visualization plots.

a) Input:  <br />
A_gene_list_file<br />
processedSampleAttributeMDCX.txt <br />
processedExpressionMDCX.txt<br />
processedGroupsMDCX.txt <br />
outputFolder/
<br /><br />Remarks:  For a gene list file input, this can be a list of gene name or probe Id or genes’ row number based on the output results from MDCoX program. If you have duplicate gene name, we encourage you to supply a list of row number output (GIDX) from the MDCoX.R results. Please refer to example/genelist folder for the format (i.e genelist.txt file or genelistRowNumber.txt file).

b) Command: Rscript MGeneset.R genelistRowNumber processedSampleAttributeMDCX.txt processedExpressionMDCX.txt processedGroupsMDCX.txt results/
<br />Remarks:  --help For help how to use it, type:   Rscript MGeneset.R --help

c) Output:  <br />
GenesetResults_10_52_39_09_10_xxx.txt  <br />
COVARIATES_Plot_10_52_39_09_10_xxx.png  <br />
STRATUM_Plot_10_52_39_09_10_xxx.png
<br />
<br />
<br />
Please contact the author of this program Herty Liany, email: e0146315@u.nus.edu if you encountered any problems. <br /> Thank you.
<br />#(2017)
