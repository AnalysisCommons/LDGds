# LDGds
## Calculate linkage disequilibrium from a GDS file.
Maintainer: Analysis Commons  
Version: 0.1

## Description:
Generate LD measures from genotypes in [Genomic Data Structure (GDS)](https://www.biostat.washington.edu/sites/default/files/modules/GDS_intro.pdf) format. This workflow will return LD information for a set of defined samples over a set of variants or a defined variant range. A flat file of LD values and a simple visualization are returned

#### What data are required for this workflow to run?
This workflow requires genotype files in GDS format (\*.gds) and a reference variant or genomic region.

#### What does this workflow output?
An M x N matrix of pairwise LD values where N = number of variants within the genomic region. M = 1 if a reference variant is provided. M = N otherwise.

### Workflow Inputs
**Bold** inputs are required; *Italic* inputs are optional. **NOTE**: one of [reference variant, interval, rsid file] is required. All may be provided resulting in LD values for the reference variant within the interval or for rsids.

- **this_gds_file**: [file, \*.gds] GDS file of genotypes per sample.
- *this_sample_ids_file*: [file, default = all samples] File of sample IDs desired for LD calculation. This file should contain one sample ID per line with no header.
- *this_ref_var*: [string, chr:pos] Genetic variant for which LD should be calculated. If provided, output is a row vector with pairwise LD with this variant in each row entry. Variant format should be 'chromosome:position'. Any punctuation seperator may be used. Only the first two values separated by punctuation will be considered.
- *this_rsid_file*: [file] A file with a list of rsids, one per line to calculate LD from.
- *this_interval*: [string, chr:start:end] Genomic interval for whcih LD should be calculated. If provided, LD will be calculated for only those variants falling within this interval. Interval format should be 'chromosome:start:end'. Any punctuation seperators may be used and need not match. Only the first three values separated by punctuation will be considered.
- *this_half_interval*: [int, default = 25kb] 1/2 of desired interval length if no interval is provided. When only a reference variant is provided, this value will be added and subtracted from the reference variant position to define the interval end and start, respectively. 
- *this_min_mac*: [int, default = 0] Minimum minor allele count for variant to be included in LD calculation.
- *this_max_mac*: [int, default = inf] Maximum minor allele count for variant to be included in LD calculation.
- *this_min_maf*: [int, default = 5%] Minimum minor allele frequency for variant to be included in LD calculation.
- *this_max_maf*: [int, default = 1] Maximum minor allele frequency for variant to be included in LD calculation.
- *this_ld_method*: [string, default = 'r'] LD calculation method. This value refers to the output LD values. Refer to documentation for the [SNPRelate package](https://bioconductor.org/packages/release/bioc/html/SNPRelate.html), specifically the function *snpgdsLDMat*. Possible values are: composite, correlation, r, and dprime with reasonable abbreviations accepted.
- **this_out_pref**: [string] Prefix for output file.
- *this_memory*: [int, default = 5GB] Amount of memory to request for computation in GB.
- *this_disk*: [int, default = size(this_gds_file) + 20 GB] Amount of disk space to request for computation in GB.
  
### Workflow Output

- **ld_file**: [file, \*.csv] An M x N diagonal matrix of pairwise LD values where N = number of variants within the genomic region. M = 1 if a reference variant is provided. M = N otherwise.

## Tutorial
Also included in this repository are genotypes from the [1000 Genomes Project](). We can use this data to run the LD calculation workflow.

### Prerequisites
Workflows written in WDL require an execution engine to run. We recommend Cromwell but others may work as well. This workflow also utilizes Docker, which must be installed prior to running. Java version 8 or greater must also be installed for Cromwell to run. Below are links to each software required. (Note: *calculate_LD.R can be run outside of the WDL workflow. It can be used like any other script, provided that the package dependencies are met.*) 

* [WDL](https://software.broadinstitute.org/wdl/documentation/quickstart)
* [Cromwell](http://cromwell.readthedocs.io/en/develop/)
* [Java8](https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)
* [Docker](https://www.docker.com/)

### Execution
Once you have downloaded the Cromwell .jar file, open a terminal window and clone the github repository:

```bash
git clone https://github.com/tmajaria/LDGds.git
```

Navigate to the **LDGds** directory and take a look at the contents:

```bash
cd LDGds/
ls
```

The output should be:
```
README.md
LICENSE	
Dockerfile
calculate_LD.R
LDGds.wdl
test_inputs.json
inputs
ouputs
```

Within the *inputs* directory, you'll find two files: a .gds file with genotypes and a .txt file with sample ids. These files will be used as inputs to the workflow. 


Workflow inputs are defined explicitely by a .json file. Within this file, each of the required inputs are defined by key-value pairs, matched to the variable within the WDL file. Optional inputs are also included. We can take a look at the test_inputs.json file that will be used in this tutorial.

```bash
cat test_inputs.json
```

The output should be:
```
{
	"LD_wf.this_interval": "20.33500000.33600000",
	"LD_wf.this_max_maf": 1,
	"LD_wf.this_min_maf": 0.2,
	"LD_wf.this_sample_ids_file": "inputs/YRI_sample_ids.txt",
	"LD_wf.this_ld_method": "r",
	"LD_wf.this_visualization": "T",
	"LD_wf.this_out_pref": "LDGds_test",
	"LD_wf.this_gds_file": "inputs/1KG_phase3_chr20_subset.gds"
}
```

These inputs show that we will be calculating LD within a 100kb region of chromosome 20, including variants that have a minor allele frequency greater than 20% in a specific set of individuals, those from the YRI subpopulation as defined by 1000 Genomes. The R² LD measure will be used and a visualization will be generated. Both output files will have the prefix *LDGds_test*.

At this point, we're ready to run the workflow. You will need the filepath for the Cromwell .jar file. For this tutorial, let's assume the path is: */my_path/cromwell-43.jar*. We can run the workflow using the syntax:

```bash
java -jar /my_path/cromwell-43.jar run LDGds.wdl -i test_inputs.json
```

Once the workflow is running, lots of messages will output to your terminal screen. Eventually, the workflow will finish and a final output message will be printed, something like:

```
[2019-07-05 14:24:21,11] [info] SingleWorkflowRunnerActor workflow finished with status 'Succeeded'.
{
  "outputs": {
	"LD_wf.workflow_message": "Workflow completed",
    "LD_wf.ld_file": "LDGds/cromwell-executions/LD_wf/0be34c4f-9543-431b-8b9a-71472eb09863/call-calculate_LD/execution/glob-beff34be37d6fa48589cc2bee8d5da74/LDGds_test.csv",
    "LD_wf.ld_plot": "LDGds/cromwell-executions/LD_wf/0be34c4f-9543-431b-8b9a-71472eb09863/call-calculate_LD/execution/glob-209f5e029d2b8b379d112924438467fe/LDGds_test.png"
  },
  "id": "0be34c4f-9543-431b-8b9a-71472eb09863"
}
```

Here, we can see three different outputs. The first is a message that tells us whether the workflow actually ran. There is only one case where the workflow will complete without running any computation: when no reference variant, interval, or rsid file is provided. In this case, the workflow message will be *Reference variant, interval, and rsid file were unspecified, no computation was initiated.* and all other outputs will be *NULL*.

The second output is a file containing the LD values for each variant included in the analysis. This will always be a symmetric matrix. The last output is a visualization of the LD structure with variants in higher LD having more red values and those in lower LD closer to white. The output visualization for this tutorial should look like:

![myimage](https://github.com/tmajaria/LDGds/blob/master/outputs/LDGds_test.png)






