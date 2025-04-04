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
     
   Chromosomal sex is based off of possibly XY combinations. Where if X is only present sex is female and when a Y is present sex is male  


## Usage

This program is designed to read through .gz files containing count number information at specific windows and matching these windows up with region locations of a ref_seq file. RefSeq_Curated.bed file was taken from NCBI with useful information being region name, chromosome, and start and stop information. Relevant information from the .gz file included chromosome, start and stop location of windows, and count number for each window. 

## The following information will be recorded for each reference region:  
|Chromosome|	Start|	Stop|	Region|	CN_Average|	Weighted_Avg|	CN_SD|	CN_SV|	Region_Size|	Total_Windows|	Gen_Variation|
|----------|:-----------:|:----:|:-----:|:---------:|:-----------:|:-----:|:------:|:------------:|:--------------:|:------------:|

## Compile and Run Commands:   
    g++ -std=c++11 CNPI.cpp -o CNPI -lz  
    ./CNPI -d {region.bed} -g {quickmer2.gz} 
### As a single command:
    g++ -std=c++11 CNPI.cpp -o CNPI -lz && time ./CNPI -d {sorted_annotated_region.bed} -g {quickmer2.gz file} -n {number_of_reference_regions} -o {output_file_names}

## Options   
    -d or -bed_gz_path: .bed or .bed.gz file with regions to match copy number windows up against - Required!  
    -g or -gz_path: .bed or .bed.gz file wih regions containing cn windows (quickmer2 output) - Required!  
    -n or -bed_path_rows: number of reference file - Optional  
    -c or -cn_rows: number of rows in quickmer2 file - Optional 
    -o or -output_name: preferred name of output file  
    -r or -record: returning a record from reference file to the terminal - Optional  
    -s of -average_check: Will Execute Count Abnormalities function. Input a number and all regions that are greater than the number from the expected Copy number will be saved to a file  
    -b of -sd_abnormalities: Will Execute Standard Deviation Abnormalities function. Input a number and all regions having a greater standard deviation greater than it will be saved to a file  
    -e or -chromosome: For displaying a certain chromosome to the terminal  
    -p or -start_stop; For displaying a certain range withing a chromosome to the terminal. Also need to input -e or -chromosome to run  

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
    -f File: The file of the child or patient that we are trying to plot  
    -p file1 file2: Files containing information for the parents of the child. With 2 parents and a child we can plot trio or a duo  
    -r Reference: background information regarding the chromosomes and sex of the child  
    -W MINWINDOW: The minimum amount of windows that are consecutively outside of the 1.3 and 2.7 window. Outside of this range indicates an alteration from normal copy number  
    -I WINDOWBUFFER: The Copy Number may oscillate around the 1.3 or 2.7 threshold and if it jumps within the okay region and then back out again this is a buffer of windows that can consecutively fall within the okay region before falling out again  
    -se SELECTCHRM: If you want to visualize a particular chromosome you can indicate which here  
    -start STARTPOS: If you want to plot a chromosome and specific start position  
    -stop STOPPOS: If you want to plot a chromosome and specific stop position  
    -gstat GSTAT_TXT: For including transcript regions on the visuals

  
    usage: CNPI_plotting.py -f FILE -r REFERENCE [-w MINWINDOW] [-i WINDOWBUFFER] [-se SELECTCHRM] [-start STARTPOS] [-stop STOPPOS] [-gstat GSTAT_TXT] [-p [file1 [file2]]] [-h]


## Dockers

### Linux

    docker run -v "$(pwd):/app:ro" -it jackust/cn_docker:CNPI_Linux_V1.0 

### Mac

    docker run -v "$(pwd):/app:ro" -it jackust/cn_docker:CNPI_Mac_V1.0

## AWS Lambda

CNPI.cpp has the capability of running at a less that $0.01 cost on AWS Lamda
### Steps
1. With AWS account go to: Lambda -> Functions -> Create Function
2. Choose a function name and select Containter image
3. Use the following container image URI:
   
            851725497249.dkr.ecr.us-east-2.amazonaws.com/tnt_test@sha256:7b856c1e10d19d90a878f4d51de59748127b7046121c91bad02b9a7686189060
4. Keep all other confugurations as default and select Create function
5. Once the function is created update Configurations
   - -> Configuration -> General configuration -> Increase Memory, Ephemeral storage, and Timeout as needed
   - -> Permissions -> select Role Name link -> policies -> select S3 full access -> actions -> attach to function
6. Back in your created function update -> Test - Event Json to:

        {
          "s3_bucket": "cnpi-input-bucket",
          "s3_key": ["Input Reference.bed (WGS_sorted_filtered.bed)","Input QuicKmer2.bed.gz (NA12878_sorted_filtered.bed.gz)","Name of output bucket (cnpi-bucket)"],
          "arguments":{
            "d": "Input Reference.bed (WGS_sorted_filtered.bed)",
            "g": "Input QuicKmer2.bed.gz (NA12878_sorted_filtered.bed.gz)",
            "o": "choose output label for files: /tmp/name (/tmp/NA12878)",
            "n": "40101",
            "c": "2295745",
            "r": "",
            "t": "",
            "s": "",
            "l": "",
            "u": "",
            "h": "",
            "w": ""}
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
