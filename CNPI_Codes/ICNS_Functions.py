import argparse
import pandas as pd
import os
import gzip
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from scipy.stats import wilcoxon
from scipy.stats import mannwhitneyu
import numpy as np


def CN_Gene_Matrix_Function(args):

    genotype_file = args.reference
    condition = args.condition


    print(f'The name of the reference file: {genotype_file}')

    # Read the genotype data
    gene_df = pd.read_csv(genotype_file, delimiter="\t", header=0)

    process_table = f'gene_distribution_{condition}.csv'

    # Create a DataFrame for new data
    gene_rows = pd.DataFrame(gene_df['Region'])
    
    # Create a new DataFrame with Region as index
    df = pd.DataFrame(index=gene_rows['Region'])

    # Check if the output file exists
    file_exists = os.path.isfile(process_table)

    # Create a Series for CN_Average
    gene_avg_series = gene_df['CN_Average']  # This is a Series now

    file_basename = os.path.basename(genotype_file)[:-13]

    # Read existing data if file exists
    if file_exists:
        df = pd.read_csv(process_table, index_col=0, delimiter="\t")
        # Concatenate the new row with existing data
        df[file_basename] = gene_avg_series.values
    else:
        # If the file doesn't exist, just use the new_row
        df[file_basename] = gene_avg_series.values

    # Write DataFrame to CSV
    df.to_csv(process_table, sep="\t")

    print(df)

def Creating_Percentiles_Gene_Data_Function(args):

    genotype_file = args.reference
    condition = args.condition
    qmer = args.qmer_windows

    print(f'The name of the reference file: {genotype_file}')

    # Read the genotype data
    if qmer is not None:
        gene_df = pd.read_csv(genotype_file, delimiter=",", index_col=[0,1,2])
    else:
        gene_df = pd.read_csv(genotype_file, delimiter="\t", index_col=0)

    print(gene_df)

    index_values = gene_df.index.tolist()
    print(index_values)
    

    percentile_data = f'Percentile_Data_{condition}.txt'

    df = pd.DataFrame(columns=['.01 Percentile', '99.9 Percentile', 'Mean', 'Total'])


    for row_index in range(len(gene_df)):
        row = gene_df.iloc[row_index]
        percentile_1 = np.percentile(row, .01).round(5)
        percentile_99 = np.percentile(row, 99.9).round(5)
        mean = row.mean().round(5)
        total = row.count()

        # Create a new row as a DataFrame
        new_row = pd.DataFrame({
            '.01 Percentile': [percentile_1],
            '99.9 Percentile': [percentile_99],
            'Mean': [mean],
            'Total': [total]
        })
        df = pd.concat([df, new_row], ignore_index=True)
        #print(new_row)

    print(df)

    df.index=index_values

    df.to_csv(percentile_data, sep="\t")


def Distribution_Scores_and_Cutoff_Y(args):

    percentiles = args.distribution_percentiles
    testing_data = args.testing_data
    output = args.output
    male_sex = args.males
    female_sex = args.females
    dup_val = args.duplication
    del_val = args.deletion

    chr_name = 'Chr'
    if female_sex is not None:
        sex = 'Female'
    if male_sex is not None:
        sex = 'Male'

    percentile_df = pd.read_csv(percentiles, sep="\t")
    testing_data_df = pd.read_csv(testing_data, sep="\t")
    

    print(percentile_df)
    print(testing_data_df)

    percentile_test = ['region']
    percentile_test_df = percentile_df[percentile_test]

    percentile_test_df['scores'] = 0
    print(percentile_test_df)

    chromosome_df = pd.DataFrame({'chromosomes' :['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                       'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']})

    chromosome_df['counts'] = 0

    for row_index in range(len(percentile_df)):
        if ((percentile_df.iloc[row_index, percentile_df.columns.get_loc('region')] == testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('Region')]) 
            and (percentile_df.iloc[row_index, percentile_df.columns.get_loc('.1 Percentile')] > testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('CN_Average')]) 
            and (testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('CN_Average')] < (del_val-1))
            and (testing_data_df.iloc[row_index, testing_data_df.columns.get_loc(chr_name)] == "chrX" or testing_data_df.iloc[row_index, testing_data_df.columns.get_loc(chr_name)]=="chrY")
            and male_sex is not None):
            print(f"Under .1 Percentile for {percentile_df.iloc[row_index, percentile_df.columns.get_loc('region')]}. The CN was {testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('CN_Average')]}. The .1 Percentile is {percentile_df.iloc[row_index, percentile_df.columns.get_loc('.1 Percentile')]}")
            percentile_test_df.iloc[row_index, percentile_test_df.columns.get_loc('scores')] = 2
            for rows in range(len(chromosome_df)):
                if chromosome_df.iloc[rows, chromosome_df.columns.get_loc('chromosomes')] == testing_data_df.iloc[row_index, testing_data_df.columns.get_loc(chr_name)]:
                    chromosome_df.iloc[rows, chromosome_df.columns.get_loc('counts')] += 2

        elif ((percentile_df.iloc[row_index, percentile_df.columns.get_loc('region')] == testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('Region')]) 
            and (percentile_df.iloc[row_index, percentile_df.columns.get_loc('99.9 Percentile')] < testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('CN_Average')])
            and (testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('CN_Average')] > (dup_val-1))
            and (testing_data_df.iloc[row_index, testing_data_df.columns.get_loc(chr_name)] == "chrX" or testing_data_df.iloc[row_index, testing_data_df.columns.get_loc(chr_name)]=="chrY")
            and male_sex is not None):
            print(f"Over 99.9 Percentile for {percentile_df.iloc[row_index, percentile_df.columns.get_loc('region')]}. The CN was {testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('CN_Average')]}. The 99.9 Percentile is {percentile_df.iloc[row_index, percentile_df.columns.get_loc('99.9 Percentile')]}")
            percentile_test_df.iloc[row_index, percentile_test_df.columns.get_loc('scores')] = 1
            for rows in range(len(chromosome_df)):
                if chromosome_df.iloc[rows, chromosome_df.columns.get_loc('chromosomes')] == testing_data_df.iloc[row_index, testing_data_df.columns.get_loc(chr_name)]:
                    chromosome_df.iloc[rows, chromosome_df.columns.get_loc('counts')] += 1
        
        elif ((percentile_df.iloc[row_index, percentile_df.columns.get_loc('region')] == testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('Region')]) 
            and (percentile_df.iloc[row_index, percentile_df.columns.get_loc('.1 Percentile')] > testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('CN_Average')])
            and (testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('CN_Average')] < del_val)):
            '''and (testing_data_df.iloc[row_index, testing_data_df.columns.get_loc(chr_name)] != "chrX" and testing_data_df.iloc[row_index, testing_data_df.columns.get_loc(chr_name)]!="chrY")):'''
            print(f"Under .1 Percentile for {percentile_df.iloc[row_index, percentile_df.columns.get_loc('region')]}. The CN was {testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('CN_Average')]}. The .1 Percentile is {percentile_df.iloc[row_index, percentile_df.columns.get_loc('.1 Percentile')]}")
            percentile_test_df.iloc[row_index, percentile_test_df.columns.get_loc('scores')] = 2
            for rows in range(len(chromosome_df)):
                if chromosome_df.iloc[rows, chromosome_df.columns.get_loc('chromosomes')] == testing_data_df.iloc[row_index, testing_data_df.columns.get_loc(chr_name)]:
                    chromosome_df.iloc[rows, chromosome_df.columns.get_loc('counts')] += 2

        elif ((percentile_df.iloc[row_index, percentile_df.columns.get_loc('region')] == testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('Region')]) 
            and (percentile_df.iloc[row_index, percentile_df.columns.get_loc('99.9 Percentile')] < testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('CN_Average')])
            and (testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('CN_Average')] > dup_val)):
            '''and (testing_data_df.iloc[row_index, testing_data_df.columns.get_loc(chr_name)] != "chrX" and testing_data_df.iloc[row_index, testing_data_df.columns.get_loc(chr_name)]!="chrY")):'''
            print(f"Over 99.9 Percentile for {percentile_df.iloc[row_index, percentile_df.columns.get_loc('region')]}. The CN was {testing_data_df.iloc[row_index, testing_data_df.columns.get_loc('CN_Average')]}. The 99.9 Percentile is {percentile_df.iloc[row_index, percentile_df.columns.get_loc('99.9 Percentile')]}")
            percentile_test_df.iloc[row_index, percentile_test_df.columns.get_loc('scores')] = 1
            for rows in range(len(chromosome_df)):
                if chromosome_df.iloc[rows, chromosome_df.columns.get_loc('chromosomes')] == testing_data_df.iloc[row_index, testing_data_df.columns.get_loc(chr_name)]:
                    chromosome_df.iloc[rows, chromosome_df.columns.get_loc('counts')] += 1
                
    print(percentile_test_df)
    print(chromosome_df)
    score = percentile_test_df['scores'].sum()
    print(score)

    score_file = output + "_score_breakdown.txt"
    gene_score_file = output + "_Gene_scores.txt"

    chromosome_df.to_csv(score_file, sep='\t',index=False)

    with open(score_file, 'a') as f:
        f.write(f'All\t{score}')

    percentile_test_df.to_csv(gene_score_file, sep='\t', index=None)  



def main():
    ##This is the organized distribution program
    ##Function is to process CNPI output data

    parser = argparse.ArgumentParser(description='Choose Function That you would like to perform')

    subparsers = parser.add_subparsers(dest='function', help='Available subcommands'
    '\n\nhello')

    #creating gene and chromosome breakdowns. Basis of all these functions
    dist_parser = subparsers.add_parser('Create_Scoring', help='Creating Gene and Chromosome based scoring from _Genotype.txt files.')
    dist_parser.add_argument('-d', '--distribution_percentiles', type=str, required=True, help='Distribution Percentiles')
    dist_parser.add_argument('-t', '--testing_data', type=str, required=True, help='testing_data')
    dist_parser.add_argument('-o','--output', type=str, required=False, help='output_name')
    dist_parser.add_argument('-m','--males',type=str, required=False, help='Sex')
    dist_parser.add_argument('-f','--females',type=str, required=False, help='Sex')
    dist_parser.add_argument('-l','--deletion', type=float, default=1.5, help='Deletion Value')
    dist_parser.add_argument('-u','--duplication', type=float, default=2.5, help='Duplication Value')

    #creating percentile.txt Arg Parser
    dist_parser = subparsers.add_parser('CreatingGenePercentiles', help='Finding Percentiles of a population')
    dist_parser.add_argument('-g', '--reference', type=str, required=True, help='Genotype.txt file')
    dist_parser.add_argument('-c', '--condition', type=str, required=True, help='Name of condition of percentile.txt file')
    dist_parser.add_argument('-q', '--qmer_windows', type=str, required=False, help='Parsing Qmer Data Matrix')

    #creating CN Gene Matrix Arg Parser
    dist_parser = subparsers.add_parser('CN_GeneMatrix', help='Creating CN gene matrix for all individuals in a population')
    dist_parser.add_argument('-g', '--reference', type=str, required=True, help='Genotype.txt file')
    dist_parser.add_argument('-c', '--condition', type=str, required=True, help='Name of condition of percentile.txt file')

    
    args = parser.parse_args()


    if args.function == 'Create_Scoring':

        #finding the gene and chromosome scores based upon a genotype.txt file
        '''ran as:
        python3 ICNS_Functions.py Create_Scoring -t 1kg_genomes_and_cnpi_output/NA12878_Genotype.txt -d Percentile_Data/Percentile_data.1_99.9_parents/Percentile_Data_.1_99.9_females.txt -o hi
        '''

        Distribution_Scores_and_Cutoff_Y(args)
        print("Ran the Distribution Scores and Cutoff Function")
    

    if args.function == 'CreatingGenePercentiles':

        #creating percentile files using gene matrix and input

        '''ran as:
        python3 ICNS_Functions.py CreatingGenePercentiles -g Percentile_Data/gene_distribution_filtered_male_parents.csv -c male_parents
        '''

        Creating_Percentiles_Gene_Data_Function(args)


    if args.function == 'CN_GeneMatrix':

        #creating the gene matrix that will be used by the CreatingGenePercentiles function

        '''rans as:
        python3 ICNS_Functions.py CN_GeneMatrix -g 1kg_genomes_and_cnpi_output/NA12878_Genotype.txt -c males_filtered_parents
        '''

        CN_Gene_Matrix_Function(args)


main()