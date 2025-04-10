# CNPI

Copy Number Private Investigator, CNPI, is a copy number analysis toolkit developed by the Tychele N. Turner, Ph.D. Lab at Washington University in St. Louis Medical School.

**Lead Developer:** Jack Ustanik (jacku[at]wustl.edu)

<div align="center">
    <img src=https://github.com/TNTurnerLab/CNPI/blob/main/Photos/abstract4.png  alt="Graphical Abstract" width="600" />
</div>

## Description

This program (CNPI) is designed to read through gzipped output files generated using the QuicK-mer2 program (quickmer2.gz files) that contain copy number information at specific windows and matching these windows up with start and stop locations from a bed file containing regions of interest. In the example directory, there is a file called `RefSeq_Curated.bed` that is a bed file containing coordinates for RefSeq genes in the genome. This file was generated from the Table Browser in the UCSC Genome Browser with useful information including the region name, chromosome, and start and stop information. Relevant information from the quickmer2.gz file includes chromosome, start and stop location of windows, and count number for each window. Please note the genome build for both files must be the same. In the example, we use build 38 of the human genome.

## Background

DNA Copy Number, in the human genome, is typically two as one copy of each chromosome is inherited from each parent. Most regions of the genome have a copy number of two except in males where there is a copy of one for chromosomes X and Y. This is not always the case as there are variations within everyone's genome as some regions have copy numbers deviating from the expected value. Therefore chromosomes, genes, or regions will have varied copy numbers. Still typically around one or two however they can vary as duplications and deletions are copy number variants (CNVs). CNVs can be in inherited trait from a parent or they can be de novo in the child. Our goal is to identify abnormal regions of a person’s genome based upon all this information given copy number data and reference regions provided by the user of this program. Genotyping across regions, karyotyping, summary statistics, and images can be produced to do all these tasks. The hope is to figure out which parent the CNVs are coming from and alterations that lead to phenotypic consequences. 

## Goals of CNPI

1. Scanning through both the region of interest file (RefSeq file in the example) and the .gz file (QuicK-mer2 output file) and creating a copy number tool for each region. Windows are generally smaller than the regions so many windows will be included within each region location. 
2. Statistics will be carried out for each region including copy number average and weighted average (weight is the length of windows included in a region), standard deviation, coefficient of variation, size of region, and total amount of windows included within each region.
3. Karyotyping of individual based upon the average of all CN windows for an entire chromosome.
4. Identify Chromosomal Sex based off of X and Y chromosome copy number.
5. It is thought that region locations with an abnormal high or low count number should be investigated more and could be indicative of phenotypic consequences.
6. In the case of count numbers that are abnormally high these regions should be recorded with the possibility of a closer look into why.

## Example File input
     
   [Example Data](https://github.com/TNTurnerLab/CNPI/tree/main/Example_Data/Example_Reference_Files)

<div align="center">

## Reference.bed Example
   
|chromosome | start | stop | Region Name |
|---------|:---:|:---:|:-----------------:|
|chr1	|11873	|14409	|DDX11L1|
|chr1	|14361	|29370	|WASH7P|
|chr1	|17368	|17436	|MIR6859-1|
|chr1	|29773	|35418	|MIR1302-2HG|
|chr1	|30365	|30503	|MIR1302-2|
|chr1	|34610	|36081	|FAM138A|

Any Tab Delimited File in this format will work as long as the first 4 columns are as shown above 

## Quickmer.bed Example

|chromosome | start | stop | Copy Number Estimate |
|---------|:---:|:---:|:-----------------:|
|chr1	|0	|54484	|2.081097|
|chr1	|54484	|60739	|2.930447|
|chr1	|60739	|68902	|2.506133|
|chr1	|68902	|82642	|2.202954|
|chr1	|82642	|88348	|2.219991|
|chr1	|88348	|108310	|2.107214|

</div>


## Overall Notes
   - Genetic Variation (Gen_variation) options include REFERENCE, DUPLICATION, DELETION, TRISOMY  
   - Genetic Variation (Gen_variation) is based on a copy number of two for chromosomes 1-22 in all individuals and for females on the X chromosome, and one for males on the X and Y chromosomes. For a duplicated or deleted chromosome, the program bases the Gen_variation on chromosome 1 or 3.
     
   - The following commands are used to sort and filter .bed and .bed.gz files for either reference or genome files   
     - Sorting to ensure data is in order for maximum speed. Filtering to get rid of unnecessary information within files.
      
   - For Sorting and filtering .bed files
     
    grep -E '^(chr[1-9][0-9]*|chrX|chrY)\b' File.bed | sort -k1,1V -k2,2n > sorted_filtered.bed
     
   - For Sorting and filtering .bed.gz files

    zgrep -E '^(chr[1-9][0-9]*|chrX|chrY)\b' File.bed.gz | sort -k1,1V -k2,2n | gzip > sorted_filtered.bed.gz

   - For Sorting and filtering GFF3/GTF files

    gzcat Unsorted.gff3.gz | grep -v "#" | awk '$1 ~ /^(1?[0-9]|2[0-2]|X|Y)$/' | sort -k1,1V -k2,2n | gzip > GFF3_sorted_filtered.gz

   - For Sorting and Filtering VCF files

    cat file.vcf | grep -v "#" | grep -v "IMPRECISE"| grep -v "LowQual" | grep -Ev "MAPQ=[0-9]($|[^0-9])|MAPQ=1[0-9]($|[^0-9])" | awk '$1 ~ /^chr(1?[0-9]|2[0-2]|X|Y)$/' | sort -k1,1V -k2,2n > VCF_sorted_filtered.txt
         
   - Sorted and Filtered GFF3, GTF, and VCF files still need to be parsed with the Parsing.cpp program

    g++ -std=c++11 Parsing.cpp -o parsing -lz && ./parsing [-c] [-g] [v] [m]

    -c: Indicates that file is a QuicK-mer2 or Bed File
    -g: Indicates that file is a GTF or GFF3 file
    -v: Indicates that file is a VCF file
    -m: If included the following label is added to chromosomes in parsed VCF output: chr
    
   Chromosomal sex is based off of possibly XY combinations. Where if X is only present sex is female and when a Y is present sex is male  


## Usage

This program is designed to read through .gz files containing count number information at specific windows and matching these windows up with region locations of a ref_seq file. RefSeq_Curated.bed file was taken from NCBI with useful information being region name, chromosome, and start and stop information. Relevant information from the .gz file included chromosome, start and stop location of windows, and count number for each window. 

## The following information will be recorded for each reference region:  
|Chromosome|	Start|	Stop|	Region|	CN_Average|	Weighted_Avg|	CN_SD|	CN_SV|	Region_Size|	Total_Windows|	Gen_Variation|
|----------|:-----------:|:----:|:-----:|:---------:|:-----------:|:-----:|:------:|:------------:|:--------------:|:------------:|

## Compile and Run Commands:   
    g++ -std=c++20 CNPI.cpp -o CNPI -lz  
    ./CNPI -d {region.bed} -g {quickmer2.gz} 
### As a single command:
    g++ -std=c++20 CNPI.cpp -o CNPI -lz && time ./CNPI -d {sorted_annotated_region.bed} -g {quickmer2.gz file} -n {number_of_reference_regions} -o {output_file_names}

## Options   
    -d or -bed_gz_path: .bed or .bed.gz file with regions to match copy number windows up against - Required!  
    -g or -gz_path: .bed or .bed.gz file with regions containing cn windows (quickmer2 output) - Required!  
    -n or -bed_path_rows: number of rows in reference file - Optional  
    -c or -cn_rows: number of rows in quickmer2 file - Optional 
    -o or -output_name: preferred name of output file  
    -r or -record: returning a record from reference file to the terminal - Optional  
    -s of -average_check: Will Execute Count Abnormalities function. Input a number and all regions that are greater than the number from the expected Copy number will be saved to a file  
    -b of -sd_abnormalities: Will Execute Standard Deviation Abnormalities function. Input a number and all regions having a greater standard deviation greater than it will be saved to a file  
    -e or -chromosome: For displaying a certain chromosome to the terminal  
    -p or -start_stop; For displaying a certain range withing a chromosome to the terminal. Also need to input -e or -chromosome to run
    -w or -weighted: Using weighted average for Gen_Variation
    -l or -deletion: Customizable threshold for deletion value. Default is 1.3
    -u or -duplication: Custimizable threshold for duplication value. Default is 2.7

<div align="center">

## Example Genotype.txt output
|Chr | Start | Stop | Region | CN_Average | Weighted_Avg | CN_SD | CN_CV | Region_Size | Total_Windows | Gen_Variation |
|--- |:-----:|:----:|:------:|:----------:|:------------:|:-----:|:-----:|:-----------:|:-------------:|:-------------:|
|chr1 | 11873 | 14409 | DDX11L1 | 2.081 | 2.081 | nan | nan | 2536 | 1 | REFERENCE |
|chr1 |14361 |29370	|WASH7P	|2.081	|2.081	|nan	|nan	|15009	|1	|REFERENCE|
|chr1 |17368 |17436	|MIR6859-1	|2.081	|2.081	|nan	|nan	|68	|1	|REFERENCE|
|chr1 |29773 |35418	|MIR1302-2HG	|2.081	|2.081	|nan	|nan	|5645	|1	|REFERENCE|
|chr1 |30365 |30503	|MIR1302-2	|2.081	|2.081	|nan	|nan	|138	|1	|REFERENCE|
|chr1 |34610 |36081	|FAM138A	|2.081	|2.081	|nan	|nan	|1471	|1	|REFERENCE|
|chr1 |65418 |71585	|OR4F5	|2.355	|2.203	|0.214	|9.087	|6167	|2	|REFERENCE|

## Example Karyotype.txt

### 46,XY
|Chr | Start_Pos | End_Pos | CN_Avg | Tot_Windows| SD|
|---|:-------:|:-------:|:-----:|:---------:|:---:|
|chr1	|0	|248944960	|2.007	|179591	|0.255|
|chr2	|0	|242183392	|1.995	|195597	|0.242|
|chr3	|0	|198173712	|1.984	|162500	|0.328|
|chr4	|0	|190191424	|1.948	|155341	|0.311|
|chr5	|0	|181356432	|1.984	|145554	|0.297|
|chr6	|0	|170739136	|1.987	|139223	|0.288|
|chr7	|0	|159333744	|1.981	|122770	|0.281|
|chr8	|0	|145076864	|1.977	|118263	|0.277|
|chr9	|0	|138233904	|1.984	|90671	|0.274|
|chr10	|0	|133787152	|1.997	|106366	|0.272|
|chr11	|0	|135076032	|1.985	|107079	|0.268|
|chr12	|0	|133263864	|1.995	|106197	|0.265|
|chr13	|0	|114351408	|1.946	|80972	|0.266|
|chr14	|0	|106883664	|1.999	|71915	|0.263|
|chr15	|0	|101980848	|2.027	|63245	|0.266|
|chr16	|0	|90215248	|2.031	|59384	|0.415|
|chr17	|0	|83245312	|2.051	|59497	|0.413|
|chr18	|0	|80261320	|1.967	|63126	|0.405|
|chr19	|0	|58607060	|2.012	|39356	|0.402|
|chr20	|0	|64333792	|2.012	|50089	|0.400|
|chr21	|0	|46681576	|2.021	|27792	|0.402|
|chr22	|0	|50799832	|2.043	|26463	|0.403|
|chrX	|0	|156029984	|1.000	|114705	|1.042|
|chrY	|0	|56884848	|1.028	|10049	|1.016|
   
</div>

## Process
[How CNPI Code Runs](https://github.com/TNTurnerLab/CNPI/blob/main/Photos/Algorithm_Process.svg)


# Python Plotting

Paired with C++ codes for custom visualization of copy number and genotype output

## Required
  - Karyotype file input to tell the program if it is 46XY 46XX etc..   
  - quickmer.gz file where windows stats will computed against  

## Goals
1. Using python to start creating figures based on copy number data  
2. Plotting program can output all values outside of a certain threshold  
3. Can also save regions or a certain number that are consecutively outside of the threshold. Up to 3 fails within the chunks  
   writes this information to a file  
4. Creates images based off of specific ranges that are specified

Files need to be the same length when plotting a trio
Use sort and filter commands to get rid of unnecessary lines

## Usage

### Ran As:
    python3 CNPI_plotting.py -f sorted_filtered.bed.gz -r Karyotype.txt
### Possible Commands
    -f File: Required: An Individual's Quickmer Copy Number Data File: Required  
    -p file1 file2: For Plotting Duos and Trios. One or Two separate files: Files need to be the same length when plotting a trios. Use sort and filter commands to get rid of unnecessary lines  
    -r Reference: Background information regarding the sex of the individual: Found in the the Karyotype file: Required  
    -W MINWINDOW: The minimum amount of windows that are consecutively outside of the 1.5 and 2.5 window. Outside of this range indicates an alteration from normal copy number  
    -I WINDOWBUFFER: The Copy Number may oscillate around the 1.5 or 2.5 threshold and if it jumps within the okay region and then back out again this is a buffer of windows that can consecutively fall within the okay region before falling out again  
    -se SELECTCHRM: If you want to visualize a particular chromosome you can indicate which here: Optional: start and stop positions requred 
    -start STARTPOS: Start Position of user specified range. Must also include Chromosome to plot: Optional 
    -stop STOPPOS: Stop Position of user specified range. Must also include chromosome to plot: Optional  
    -gstat GSTAT_TXT: For complementing plots with Genotype_Stats.txt: Optional for seeing region statistics
    -o output: Specifying Output Image Label
    -min minimum: For changing default minimum abnormal value (1.5 default)
    -max --maximum: For changing default maximum abnormal value (2.5 default)

  
    usage: CNPI_plotting.py -f FILE -r REFERENCE [-w MINWINDOW] [-i WINDOWBUFFER] [-se SELECTCHRM] [-start STARTPOS] [-stop STOPPOS] [-gstat GSTAT_TXT] [-p [file1 [file2]]] [-h]

# Individual-level Copy Number Scores (ICNS scoring)

Python program that scores individuals based on normal copy number percentile distributions.
It assigns scores to regions within an individual's Genotype.txt file that show abnormally high or low copy number values. A score of 2 is given to regions significantly below the normal percentile range, while a score of 1 is given to regions significantly above it. The program provides both per-region and per-chromosome scoring breakdowns to help identify abnormal genomic areas and distinguish affected chromosomes.

## Usage

### Ran As:
    python3 ICNS_Functions.py [Create_Scoring] or [CreatingGenePercentiles] or [CN_GeneMatrix]

### Possible Commands:
CN_GeneMatrix: Function for concatonation Copy Number Average Scores from Genotype.txt files that can later be used for making Percentile Values

    -g reference: Genotype.txt file to be fed into Copy Number Average Matrix
    -c condition: For naming CN Matrix

CreatingGenePercentiles: Finding Percentiles of regions within Copy Number Matrix file

    -g reference: Matrix of Copy Number average data for creating CN percentiles
    -c condition: For naming CN percentile file

Create_Scoring: Function that scores individuals based upon region Copy Number Percentile files.

    -d distribution_percentiles: File created from the CreatingGenePercentiles function. Copy Number percentiles of regions
    -t testing_data: Genotype.txt file to create ICNS score for based off of Copy Number Percentiles file
    -o output: Providing a specific name for output files
    -m males: Adjusting Scoring for males for X and Y chromosomes
    -f females: Skipping scoring for Y chromosome in females
    -l deletion: Deletion cutoff value for scoring. Default is 1.5
    -u duplication: Duplication cutoff value for scoring. Default is 2.5

## Dockers

### Linux

    docker run -v "$(pwd):/app:ro" -it jackust/cn_docker:CNPI_Linux_V1.0 

### Mac

    docker run -v "$(pwd):/app:ro" -it jackust/cn_docker:CNPI_Mac_V1.0

## AWS Lambda

CNPI.cpp has the capability of running for less than $0.01 per sample on AWS Lamda
### Steps
1. With AWS account go to: Lambda -> Functions -> Create Function
2. Choose a function name and select Container image
3. Use the following container image URI:
   
            851725497249.dkr.ecr.us-east-2.amazonaws.com/tnt_test@sha256:bcf1b5ba8b6adbde47aa935b8cd27c9386f2fa3fa4a4b80abcfc1384329d9302
4. Keep all other configurations as default and select Create function
5. Once the function is created update Configurations
   - -> Configuration -> General configuration -> Increase Memory, Ephemeral storage, and Timeout as needed
   - -> Permissions -> select Role Name link -> policies -> select AmazonS3FullAccess -> actions -> attach to function
6. Back in your created function update -> Test - Event Json to:

        {
          "s3_bucket": "cnpi-input-bucket",
          "s3_key": ["Input Reference.bed (WGS_sorted_filtered.bed)","Input QuicKmer2.bed.gz (NA12878_sorted_filtered.bed.gz)","Name of output bucket (cnpi-bucket)"],
          "arguments":{
            "d": "Input Reference.bed (WGS_sorted_filtered.bed)",
            "g": "Input QuicKmer2.bed.gz (NA12878_sorted_filtered.bed.gz)",
            "o": "choose output label for files: /tmp/name (/tmp/NA12878)",
            "n": "optional argument",
            "c": "optional argument",
            "r": "optional argument",
            "t": "optional argument",
            "s": "optional argument",
            "l": "optional argument",
            "u": "optional argument",
            "h": "optional argument",
            "w": "optional argument"}
        }


# Example Pictures

## Plotting of a Region Outside of Normal (1.3 - 2.7) Threshold

 - Dash Lines at the bottom represent locations of quickmer window readings
 - Colored boxes within graph represent different genes corresponding to the genes within the table on the right
 - Red dash lines representing within normal copy number range (1.3 - 2.7)
 - The blue line in the middle at 2 representing a copy number value of 2
 - Plotted blue line throughout the graph representing the copy number across chromosome

<div align="center">


### Duplication Event Along Chromosome 16
![Duplication Event Along Chromosome](https://github.com/TNTurnerLab/CNPI/blob/main/Photos/Plotting/Chr16_NA12878_Duplication.png)

</div>

## Plotting Based off an Inputted Range

- Red dash lines representing within normal copy number range (1.3 - 2.7)
- The blue line in the middle at 2 representing a copy number value of 2
- Plotted blue line throughout the graph representing the copy number across inputted range

<div align="center">


### Copy Number Across Chromosome 6
![Normal Chromosome 6](https://github.com/TNTurnerLab/CNPI/blob/main/Photos/Plotting/NA12878_Chr6.png)

### Abnormal Copy Number Across Chromosome 12
![Entire Duplicated Chromosome](https://github.com/TNTurnerLab/CNPI/blob/main/Photos/Plotting/NA12739_chr12.png)

</div>
