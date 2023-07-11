
# miRNAs_SRP

  <p align="center">
    A pipeline for global analysis of miRNA-mediated stress response in plants
    <br />
    <a href="https://github.com/antoglz/miRNAs_SRP"><strong>Explore the docs »</strong></a>
    <br />
    <br /> 
    <a href="https://github.com/antoglz/miRNAs_SRP/issues">Report Bug</a>
    ·
    <a href="https://github.com/antoglz/miRNAs_SRP/issues">Request Feature</a>
  </p>
</p>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li> <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#languages">Languages</a></li>
        <li><a href="#pipeline-description">Pipeline Description</a></li>
        <ul>
            <li><a href="#phase-i">Phase I: Download, pre-processing and filtering</a></li>
            <li><a href="#phase-ii">Phase II: Analysis of small RNAs</a></li>
            <li><a href="#phase-iii">Phase III: Identification and classification of stress-responsive miRNAs</a></li>
        </ul>
      </ul>
    <li><a href="#requirements">Requirements</a></li>
      <ul>
        <li><a href="#software">Software</a></li>
        <li><a href="#standardized-metadata-structure">Standardized Metadata Structure</a></li>
        <li><a href="#additional-information">Additional Information</a></li>
      </ul>
    </li>
    <li> <a href="#usage">Usage</a>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>


## About The Project
<p style="text-align: justify;">
This repository contains the scripts and part of the additional files needed to perform a global analysis of the miRNA-mediated stress response in plants. These scripts are part of a pipeline for downloading, processing, filtering and analysing miRNA-seq libraries deposited in the <a href=https://www.ncbi.nlm.nih.gov/sra>Sequence Read Archive (SRA)</a> of the National Center for Biotechnology Information (NCBI). This process involves the identification of stress-responsive miRNAs in multiple and diverse plant species, their classification according to their range of stress response (General, Intermediate and Specific) and the generation of the data required to construct the miRNA-mediated stress response network in plants.
</p>


### Languages

- [Bash](https://www.gnu.org/software/bash/manual/bash.html) (v.4.4.20)
- [Python](https://www.python.org) (v.3.8.13)
- [R](https://www.r-project.org) (v.4.0.0)

### Pipeline Description

<p style="text-align: justify;">
The scripts deposited in this repository are intended to be run on a high-performance computing (HPC) cluster using the <a href=https://slurm.schedmd.com/documentation.html>Simple Linux Utility for Resource Management (SLURM)</a> scheduling manager. Two types of scripts are used: main scripts and submission scripts. The main scripts carry out all tasks related to data management, filtering, and analysis. Meanwhile, submission scripts (in bash) are responsible for requesting the necessary resources and executing the main script with which it is associated. Submission scripts contain the suffixes "parallel" or "exe". Those scripts with the "parallel" suffix execute the main script simultaneously for the data of multiple plant species, while those with the "exe" suffix execute the script sequentially. The pipeline can be divided into three phases: downloading, pre-processing and filtering of data (orange), analysis of small RNAs (green) and identification and classification of stress-responsive miRNAs (red).
</p>

<p align="center">
  <img src="https://github.com/antoglz/miRNAs_SRP/blob/main/workflow_img.png" class="center" height="600" width="800" >
</p>

<h4 id="phase-i">Phase I: Download, pre-processing and filtering</h4>

<p style="text-align: justify;">
This pipeline starts with <code>01-Prefetch.sh</code>, a script designed to download the miRNA-seq libraries and perform part of the quality control. These libraries are in compressed FASTQ format and are used as input for <code>02-Trimming.sh</code>. This script is responsible for cleaning and filtering the sequences according to predefined criteria, also removing any remaining adapters from the sequencing process. The resulting libraries are passed through <code>03-Filter_by_depth_rep.py</code>, which converts them to FASTA format and filters them according to their sequencing depth and the number of replicates in their sample group. The <code>04-Counts_and_sRNADatabase.py</code> script filters the generated FASTA files by removing rRNA, tRNA, pnRNA and pnoRNA sequences.
</p>

<h4 id="phase-ii">Phase II: Analysis of small RNAs</h4>

<p style="text-align: justify;">
The <code>04-Counts_and_sRNADatabase.py</code> script also counts the remaining sequences in each library, creating two types of absolute count tables: one that includes the sequences present in at least one of the study libraries and one that contains only those sequences present in all the study libraries. The resulting tables can become extremely large depending on the study, which can be a problem for further analysis. To address this situation, <code>05-Divide_filter_by_experiments.py</code> fragments the tables of those studies that consider multiple factors or different levels of a single factor. Each of the substudies generated includes a control group and all the treatment groups with which it is contrasted, thus reducing the complexity of the tables. The <code>06-Diff_exp_analysis.r</code> script uses the tables including the sequences present in all the substudy libraries to perform an exploratory analysis of the data. It then selects the substudies considered valid in this analysis and uses the tables with the sequences that appear in at least one of the substudy libraries to carry out a differential expression analysis.
</p>

<h4 id="phase-iii">Phase III: Identification and classification of stress-responsive miRNAs</h4>
<p style="text-align: justify;">
The differentially expressed sequences of each substudy are annotated using <code>07-miRNAs_identification.sh</code>. The identified miRNAs are grouped into families by <code>08-Group_miRNAs_by_family.r</code>, also checking whether the expression of miRNAs belonging to the same family shows a common expression pattern or, on the contrary, follows different trends. Finally, <code>09-Build_presence_absence_table.py</code> creates a binary table representing the stressful situations in which at least one member of the annotated miRNA families exhibits differential expression. This table is used by the <code>10-Analysis_results.r</code> script to generate the results that were used to construct the miRNA stress response network in plants.
</p>

## Requirements

### Software

#### External Software

- [Bowtie](https://bowtie-bio.sourceforge.net/index.shtml) (v.1.3.1)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (v.2.3.5.1)
- [fastp](https://github.com/OpenGene/fastp) (v.0.21.0)
- [FastQC](https://github.com/s-andrews/FastQC) (v.0.11.9)
- [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki) (v.2.11.2)

#### Python Libraries

[argparse](https://docs.python.org/es/3/library/argparse.html) (v.1.1), [Biopython](https://biopython.org) (v.1.80), [Numpy](https://numpy.org) (v.1.23.1), [Pandas](https://pandas.pydata.org) (v.1.4.0rc0), [SQLite3](https://docs.python.org/3/library/sqlite3.html) (v.3.39.3)

#### R Libraries

[argparse](https://cran.r-project.org/web/packages/argparse/index.html) (v.2.2.2), [Ckmeans.1d.dp](https://github.com/cran/Ckmeans.1d.dp) (v.4.3.4), [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) (v.2.10.0), [DESeq2](https://github.com/mikelove/DESeq2) (v.1.38.2), [dplyr](https://dplyr.tidyverse.org/index.html) (v.1.1.1), [ff](https://github.com/truecluster/ff) (v.4.0.9), [ggplot2](https://ggplot2.tidyverse.org) (v.3.4.1), [ggthemes](https://cran.r-project.org/web/packages/ggthemes/index.html) (v.4.2.4), [grid](https://github.com/cran/grid/tree/master) (v.4.1.2), [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) (v.2.3), [htmltools](https://github.com/rstudio/htmltools) (v.0.5.5), [Plotly](https://plotly.com/r/) (v.4.10.1), [stringr](https://stringr.tidyverse.org) (v.1.5.0), [tibble](https://tibble.tidyverse.org) (v.3.2.1)

### Standardized Metadata Structure
<p style="text-align: justify;">
Studies stored in SRA have associated metadata that provide descriptive information about their samples. Although there are certain fields of information common to all studies, many of them differ. This fact complicates the automation of data processing and analysis. Therefore, a standardized structure was designed that must be applied to all studies to be used in the analysis. This structure uses the same information fields for all metadata and follows specific guidelines to represent the information for each sample and the relationships between them. The metadata tables are made up of 13 columns that include information on the plant species studied (<code>Species</code>), the study to which each sample belongs (<code>Project</code>) and its own identifier (<code>Run</code>), whether the sample belongs to a treatment or control group (<code>Condition</code>), the replicate number (<code>Replicate</code>), the type (<code>Stress</code>) and level of stress (<code>Level</code>), the exposure time (<code>Time</code>), the plant cultivar (<code>Cultivar</code>), the tissue used (<code>Tissue</code>), the age or stage of development of the plant (<code>Age/Dev-state</code>), the genotype (<code>Genotype</code>: whether resistant or susceptible) and the type of library used (<code>Library</code>: paired-end or single-end).
</p>

In many studies different control groups can be found that will be contrasted with a specific set of treatment groups. Therefore, it is necessary to specify in these cases which groups of samples are related. This is indicated in the <code>Stress</code>, <code>Level</code> and <code>Time</code> columns of the control samples by using the ":" character, which is used to separate the level of the factor associated with the control samples from the level of the factor associated with the treatment samples to be compared (e.g., <code>22°C:4°C</code>). In addition, the "-" character is used to separate words within the same cell, the "&" character is used to separate the names of different stress-inducing agent when they are applied in combination (e.g., <code>cold&salt</code>) and the term "neg" is used when negative values are represented. Finally, when no information is available for any of the fields used, this is indicated by the column identifier (<code>Level</code>: L, <code>Time</code>: T, <code>Cultivar</code>: CV, <code>Age/Dev-stage</code>: D, <code>Genotype</code>: G, <code>Library</code>: Li) followed by a period and a 0 (e.g., <code>L.0</code>). The columns <code>Specie</code>, <code>Project</code>, <code>Run</code>, <code>Condition</code>, <code>Stress</code>, <code>Replicate</code> and <code>Tissue</code> should always contain information about the samples.

Here is a simple example:

|Species |Project  |Run  |Condition  | Stress | Replicate | Level | Time | Cultivar | Tissue | Age/Dev-state | Genotype | Library |
|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
| arth | PRJNAXXXXX | SRR0000001 | treated | drought | 1 | 30%SWC | 6h | Col-0 | leaves | seedling | G.0 | single
arth |PRJNAXXXXX | SRR0000002| treated | drought | 2 | 30%SWC | 6h | Col-0 | leaves | seedling | G.0 | single
arth |PRJNAXXXXX | SRR0000003| treated | drought | 3 | 30%SWC | 6h | Col-0 | leaves | seedling | G.0 | single
arth |PRJNAXXXXX | SRR0000004| control | drought | 1 | 70%:30%SWC | 6h | Col-0 | leaves | seedling | G.0 | single
arth |PRJNAXXXXX | SRR0000005| control | drought | 2 | 70%:30%SWC | 6h | Col-0 | leaves | seedling | G.0 | single
arth |PRJNAXXXXX | SRR0000006| control | drought | 3 | 70%:30%SWC | 6h | Col-0 | leaves | seedling | G.0 | single
arth |PRJNAXXXXX | SRR0000007| treated | drought | 1 | 20%SWC | 6h | Col-0 | leaves | seedling | G.0 | single
arth |PRJNAXXXXX | SRR0000008| treated | drought | 2 | 20%SWC | 6h | Col-0 | leaves | seedling | G.0 | single
arth |PRJNAXXXXX | SRR0000009| treated | drought | 3 | 20%SWC | 6h | Col-0 | leaves | seedling | G.0 | single
arth |PRJNAXXXXX | SRR0000010| control | drought | 1 | 70%:20%SWC | 6h | Col-0 | leaves | seedling | G.0 | single
arth |PRJNAXXXXX | SRR0000011| control | drought | 2 | 70%:20%SWC | 6h | Col-0 | leaves | seedling | G.0 | single
arth |PRJNAXXXXX | SRR0000012| control | drought | 3 | 70%:20%SWC | 6h | Col-0 | leaves | seedling | G.0 | single

### Additional Information

There are three FASTA files that have not been deposited in this repository due to space constraints, but are essential for executing this pipeline. Two of these files are used for annotating differentially expressed sRNA sequences and contain mature plant miRNA sequences deposited in the miRBase and PmiREN databases. The other one includes RNA sequences from the RNACentral database and is used to filter out unwanted RNA sequences such as rRNA, tRNA, snRNA, and snoRNA.

On the other hand, the accession lists of the studies used in the analysis, the standardised metadata for each of them and the <code>species_id.csv</code> file are made available in the <code>Additional_files</code> directory. In this study, each plant species is identified by a unique 4-letter code (e.g., <i>Arabidopsis thaliana</i> = arth). The <code>species_id.csv</code> file establishes correspondence between these codes, the identifiers used in the databases mentioned above and the scientific names of the species. It is worth mentioning that this file only contains information on the species included in the analysis, therefore, it should be updated in case new species are added to the study.

## Usage

Run any of the submission scripts using the SLURM workload manager:

```
sbatch script_example_parallel.sh
```

or

```
sbatch script_example_exe.sh
```

## Contact

Antonio González Sánchez - gonsanan@alumni.uv.es
