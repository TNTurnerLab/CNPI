import pandas as pd
import gzip
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import numpy as np
import math
import itertools
import sys
from collections import deque
from matplotlib.collections import LineCollection



def go_plot(combined_data, kary_array, parent1_type, parent2_type, gstats, gstat_array, select_chrm, skip, parent_len, start_pos, stop_pos, outputs, output, del_cn, dup_cn):

    print("The combined data is")


    plt.rcParams["font.family"] = "Arial"
    plt.rcParams['font.size'] = 16
    
    # Convert kary_array to Pandas DataFrame
    df = pd.DataFrame(combined_data, columns=['Chr', 'Lowest_Start_Pos', 'Highest_End_Pos', 'CN'])
    last_row = df.iloc[-1]

    if(skip):
        df2 = pd.DataFrame(combined_data, columns=['Parent1_Chr', 'Parent1_Lowest_Start_Pos', 'Parent1_Highest_End_Pos', 'Parent1_CN'])
        last_row2 = df2.iloc[-1]
        if(parent_len != 1):
            df3 = pd.DataFrame(combined_data, columns=['Parent2_Chr', 'Parent2_Lowest_Start_Pos', 'Parent2_Highest_End_Pos', 'Parent2_CN'])
            last_row3 = df3.iloc[-1]


    # Create a new row with the extracted values
    new_row = pd.DataFrame({
        'Chr': [last_row['Chr']],
        'Lowest_Start_Pos': [last_row['Highest_End_Pos']],
        'Highest_End_Pos': [last_row['Highest_End_Pos']],  # Assuming this should be same as Lowest_Start_Pos
        'CN': [last_row['CN']]
    })
    df = pd.concat([df, new_row], ignore_index=True)
    
    if(skip):
    # Create a new row with the extracted values
        new_row2 = pd.DataFrame({
            'Parent1_Chr': [last_row2['Parent1_Chr']],
            'Parent1_Lowest_Start_Pos': [last_row2['Parent1_Highest_End_Pos']],
            'Parent1_Highest_End_Pos': [last_row2['Parent1_Highest_End_Pos']],  # Assuming this should be same as Lowest_Start_Pos
            'Parent1_CN': [last_row2['Parent1_CN']]
        })
        df2 = pd.concat([df2, new_row2], ignore_index=True)

        if(parent_len !=1):
        # Create a new row with the extracted values
            new_row3 = pd.DataFrame({
                'Parent2_Chr': [last_row3['Parent2_Chr']],
                'Parent2_Lowest_Start_Pos': [last_row3['Parent2_Highest_End_Pos']],
                'Parent2_Highest_End_Pos': [last_row3['Parent2_Highest_End_Pos']],  # Assuming this should be same as Lowest_Start_Pos
                'Parent2_CN': [last_row3['Parent2_CN']]
            })
            df3 = pd.concat([df3, new_row3], ignore_index=True)

    # Convert numeric columns to appropriate types
    numeric_cols = ['Lowest_Start_Pos', 'Highest_End_Pos', 'CN']
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce')
    if(skip):
        numeric_cols1 = ['Parent1_Lowest_Start_Pos', 'Parent1_Highest_End_Pos', 'Parent1_CN']
        df2[numeric_cols1] = df2[numeric_cols1].apply(pd.to_numeric, errors='coerce')
        if(parent_len!=1):
            numeric_cols2 = ['Parent2_Lowest_Start_Pos', 'Parent2_Highest_End_Pos', 'Parent2_CN']
            df3[numeric_cols2] = df3[numeric_cols2].apply(pd.to_numeric, errors='coerce')

    highest_y_chr = df['CN'].max()
    if(skip):
        highest_y_parent1 = df2['Parent1_CN'].max()  # Get maximum y value for parent1_range
        if(parent_len!=1):
            highest_y_parent2 = df3['Parent2_CN'].max()  # Get maximum y value for parent2_range}
    
    

    highest_y = highest_y_chr
    
    if gstats:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))  # 1 row, 2 columns

    
    if not gstats:
        fig, ax1 = plt.subplots(figsize=(22, 8))

    # Plotting with Seaborn
    sns.lineplot(x='Lowest_Start_Pos', y='CN', data=df, color='b', linestyle='-', linewidth=.4, label="Individual", ax=ax1)
    if(skip):
        sns.lineplot(x='Parent1_Lowest_Start_Pos', y='Parent1_CN', data=df2, color='green', linestyle='-', linewidth=.5, label=parent1_type, ax=ax1)
        highest_y = max(highest_y_chr, highest_y_parent1)
        if(parent_len!=1):
            sns.lineplot(x='Parent2_Lowest_Start_Pos', y='Parent2_CN', data=df3, color='darkorange', linestyle='-', linewidth=.5, label=parent2_type, ax=ax1)
            highest_y = max(highest_y_chr, highest_y_parent1, highest_y_parent2)


    # Customizing the plot using Matplotlib
    ax1.set_title('Copy Number Estimates '+ kary_array[0][0], fontsize=18)
    ax1.tick_params(axis='both', labelsize=14)
    ax1.set_xlabel('Chr'+ select_chrm[3:] +' Position', fontsize=16)
    ax1.set_ylabel('Copy Number', fontsize=16)
    ax1.grid(True)

    print(f"This highest is {highest_y}")

    if(highest_y<2.5 or math.isnan(highest_y)):
        highest_y=2.5
    if(highest_y>5 or math.isnan(highest_y)):
        highest_y=5    
    print(f"This highest is {highest_y}")
    plt.ylim(bottom=0, top=highest_y+1)
    if gstats:
        ax1.set_ylim(bottom=-0.2, top=highest_y+1)
    if not gstats:
        ax1.set_ylim(bottom=0, top=highest_y+1)
    ax1.set_xlim(start_pos, stop_pos)



    print(kary_array[0][0])
    print(f'The Chromosome is {select_chrm[3:]}')
    if ((kary_array[0][0] == "46,XY" or kary_array[0][0] == "47,XXY" or kary_array[0][0] == "47,XY+9" or kary_array[0][0] == "47,XY+21" or kary_array[0][0] == "47,XY+15" or
    kary_array[0][0] == "47,XY+14" or kary_array[0][0] == "47,XY+12" or kary_array[0][0] == "47,XY+13" or kary_array[0][0] == "47,XY+18" or kary_array[0][0] == "48,XXXY") and select_chrm[3:] == "Y"):
        ax1.axhline(y=1, color='black', linestyle='-', linewidth=2, label='CN=1')
        ax1.axhline(y=dup_cn-1, color='r', linestyle='--', linewidth=2, label=f'CN={dup_cn-1}')
        ax1.axhline(y=del_cn-1, color='r', linestyle='--', linewidth=2, label=f'CN={del_cn-1}')
    elif((kary_array[0][0]=="46,XY" or kary_array[0][0] == "47,XXY" or kary_array[0][0] == "47,XY+9" or kary_array[0][0] == "47,XY+21" or kary_array[0][0] == "47,XY+15" or
    kary_array[0][0] == "47,XY+14" or kary_array[0][0] == "47,XY+12" or kary_array[0][0] == "47,XY+13" or kary_array[0][0] == "47,XY+18") and select_chrm[3:]=="X"):
        ax1.axhline(y=1, color='black', linestyle='-', linewidth=2, label='CN=1')
        ax1.axhline(y=dup_cn-1, color='r', linestyle='--', linewidth=2, label=f'CN={dup_cn-1}')
        ax1.axhline(y=del_cn-1, color='r', linestyle='--', linewidth=2, label=f'CN={del_cn-1}')
    else:
         # Add horizontal lines
        ax1.axhline(y=2, color='black', linestyle='-', linewidth=2, label='CN=2')
        ax1.axhline(y=dup_cn, color='r', linestyle='--', linewidth=2, label=f'CN={dup_cn}')
        ax1.axhline(y=del_cn, color='r', linestyle='--', linewidth=2, label=f'CN={del_cn}')
   
    # Add vertical lines at each data point
    for pos in df['Lowest_Start_Pos']:
        ax1.plot([pos, pos], [-0.3, 0.0], color='gray', linestyle='-', linewidth=1)

    labels = []
    unique_rows = []
    go_label =0
    # Add solid bars underneath the graph

    if gstats:
        colors = ['purple', 'blue', 'green', 'red', 'orange', 'cyan', 'darkviolet', 'lavender', 'dodgerblue', 'salmon', 'silver', 'chocolate',
                  'yellow', 'slategray', 'hotpink', 'lime', 'mediumturquoise', 'moccasin', 'firebrick', 'gold', 'peru', 'darkseagreen', 'slategray']
        color_cycle = itertools.cycle(colors)
        for _, row in df.iterrows():
            for line in gstat_array:
                if ((row[0] == line[0]) and 
                    ((float(line[1]) <= start_pos and float(line[2]) >= stop_pos) or
                     (float(line[1]) >= start_pos and float(line[2]) <= stop_pos) or
                     (float(line[1]) >= start_pos and float(line[1]) <= stop_pos and float(line[2]) >= stop_pos) or
                     (float(line[1]) <= start_pos and float(line[2]) >= start_pos and float(line[2]) <= stop_pos))):
                    
                    labels.append((line))

                    go_label =1

                    # Get the next color from the cycle
                    color = next(color_cycle)  # Default to 'purple' if color_cycle is exhausted
                    
                    # Plot solid bar
                    lines = [[(float(line[1]), -0.1), (float(line[2]), -0.1)]]  # list of line segments
                    lc = LineCollection(lines, linewidths=12, colors=color)  # control linewidth here
                    ax1.add_collection(lc)
    if gstats:
        ax1.set_position([0.05, 0.1, 0.5, 0.8])
        ax2.set_position([0.65, 0, 0.25, 0.8])
    if not gstats:
        ax1.set_position([0.04, 0.1, .95, 0.8])
    

    if(gstats and go_label==1):
        # To store unique rows
        # Create the plot with an extended figure size
        seen_rows = set()

        for row in labels:
            row_tuple = tuple(row)  # Convert list to tuple to make it hashable
            if row_tuple not in seen_rows:
                seen_rows.add(row_tuple)
                unique_rows.append(row)

        unique_df = pd.DataFrame(unique_rows, columns=['Chr','Start','Stop','Name','CN_Avg','Weighted_Avg','CN_SD','CN_CV','Region_Size','Regions','Gen_Variation'])
        
        plot_df = unique_df[['Name','CN_Avg','Regions','Gen_Variation']]
        plot_df = plot_df.rename(columns={'Regions': 'Windows'})


        

        # Round 'CN_Avg' and 'CN_SD' to 4 decimal places
        plot_df['Windows'] = pd.to_numeric(plot_df['Windows'], errors='coerce')
        plot_df['CN_Avg'] = pd.to_numeric(plot_df['CN_Avg'], errors='coerce')
        plot_df['CN_Avg'] = plot_df['CN_Avg'].round(4)

         # Convert DataFrame to list of lists for plt.table
        table_data = plot_df.values.tolist()
        col_labels = plot_df.columns.tolist()
        
        ax2.text(0.5, 1.1, 'Region Information', ha='center', va='center', transform=ax2.transAxes, fontsize=20)

        # Plot the table on the second subplot (ax2)
        ax2.axis('tight')
        ax2.axis('off')  # Hide axes
        table = ax2.table(cellText=table_data, colLabels=col_labels, loc='center', cellLoc='center')#, bbox=[-0.1,-0.1,1.3,1.2])
        table.set_fontsize(20)
        table.auto_set_column_width([i for i in range(len(col_labels))])

        # Force cell height to fit font
        for (row, col), cell in table.get_celld().items():
            cell.set_height(0.05)  # Adjust as needed
            cell.set_text_props(va='center', ha='center', wrap=True)
        
        
        
    current_xlim = ax1.get_xlim()
    if current_xlim[0] < 0:
        ax1.set_xlim(left=0)
    # Show legend
    ax1.legend(loc='upper left', bbox_to_anchor=(.78, 1.01), fontsize=14)

    # Display the plot
    if outputs == False:
        filename = select_chrm+"_"+str(start_pos)+"-"+str(stop_pos)+".png"
    elif outputs == True:
        filename = output+"-"+select_chrm+"_"+str(start_pos)+"-"+str(stop_pos)+".png"
    

    plt.savefig(filename)
    
    plt.show()
    sys.exit(0)
    return 0


def plot_function(consecutive_regions, df_combined,skip, case, gstats, gstat_array, unique_rows, parent1_type, parent2_type, parent_len, pre_array, post_array, del_cn, dup_cn):
    converted_array= []
    last =0

    plt.rcParams["font.family"] = "Arial"
    plt.rcParams['font.size'] = 13

    for row in consecutive_regions:
        columns = row.split('\t')
        converted_array.append(columns)

    for row in converted_array:
        last = int(row[2])

    highest_y_chr =0
    highest_y_parent1 =0
    highest_y_parent2 =0

    converted_pre_array = []
    for row in pre_array:
        columns = row.split('\t')
        converted_pre_array.append(columns)


    converted_array = converted_pre_array + converted_array + post_array
    

    for row in converted_array:
        last = int(row[2])

    # Filter data for chr1 and specific range of positions
    chr_range = df_combined[(df_combined['Chr'] == converted_array[0][0]) & 
                            (df_combined['Lowest_Start_Pos'] >= int(converted_array[0][1])) &
                            (df_combined['Highest_End_Pos'] <= last)]
    highest_y_chr = chr_range['CN'].max()

    if(skip):
        parent1_range = df_combined[(df_combined['Parent1_Chr'] == converted_array[0][0]) & 
                                    (df_combined['Parent1_Lowest_Start_Pos'] >= int(converted_array[0][1])) &
                                    (df_combined['Parent1_Highest_End_Pos'] <= last)]
        highest_y_parent1 = parent1_range['Parent1_CN'].max()  # Get maximum y value for parent1_range

        if parent_len!=1:
            parent2_range = df_combined[(df_combined['Parent2_Chr'] == converted_array[0][0]) & 
                                        (df_combined['Parent2_Lowest_Start_Pos'] >= int(converted_array[0][1])) &
                                        (df_combined['Parent2_Highest_End_Pos'] <= last)]
            highest_y_parent2 = parent2_range['Parent2_CN'].max()  # Get maximum y value for parent2_range


    highest_y = max(highest_y_chr, highest_y_parent1, highest_y_parent2)

    
    # Plotting with Seaborn
    sns.lineplot(x='Lowest_Start_Pos', y='CN', data=chr_range, color='b', linestyle='-', linewidth=.75, label='Individual')
    if(skip):
        sns.lineplot(x='Parent1_Lowest_Start_Pos', y='Parent1_CN', data=parent1_range, color='green', linestyle='-', linewidth=.75, label=parent1_type)
        if parent_len !=1:
            sns.lineplot(x='Parent2_Lowest_Start_Pos', y='Parent2_CN', data=parent2_range, color='darkorange', linestyle='-', linewidth=.75, label=parent2_type)
    

    # Customizing the plot using Matplotlib
    plt.title('Copy Number Estimates for Chromosome ' + converted_array[0][0] + " for " + case)
    plt.xlabel('Chr Start Position')
    plt.ylabel('Copy Number')
    plt.grid(True)

    if(highest_y<2.5):
        highest_y=2.5

    print(f'The highest Y is {highest_y}')
    if np.isnan(highest_y) or np.isinf(highest_y):
        highest_y=10

    plt.ylim(bottom=-0.35, top=highest_y+1)

    label1= 'CN=2'
    label2= f'CN={dup_cn}'
    label3= f'CN={del_cn}'
    y1 =  2
    y2 = dup_cn
    y3 = del_cn
    # Add horizontal lines
    if(case=="46,XY" and (converted_array[0][0] == "chrY" or converted_array[0][0] == "chrX")):
        label1 = 'CN=1'
        label2 = f'CN={dup_cn-1}'
        label3 = f'CN={del_cn-1}'
        y1 = 1
        y2 = dup_cn-1
        y3 = del_cn-1


    plt.axhline(y=y1, color='black', linestyle='-', linewidth=2, label=label1)
    plt.axhline(y=y2, color='r', linestyle='--', linewidth=2, label=label2)
    plt.axhline(y=y3, color='r', linestyle='--', linewidth=2, label=label3)

    # Add vertical lines at each data point
    for pos in chr_range['Lowest_Start_Pos']:
        plt.plot([pos, pos], [-0.3, 0], color='gray', linestyle='-', linewidth=1)


    labels = []
    go_label =0
# Add solid bars underneath the graph
    if gstats:
        colors = ['purple', 'blue', 'green', 'red', 'orange', 'cyan', 'darkviolet', 'lavender', 'dodgerblue', 'salmon', 'silver', 'chocolate',
                  'yellow', 'slategray', 'hotpink', 'lime', 'mediumturquoise', 'moccasin', 'firebrick', 'gold', 'peru', 'darkseagreen', 'slategray']
        color_cycle = itertools.cycle(colors)
        for row in converted_array:
            line_break = 0
            for line in gstat_array:
                if line_break ==1:
                    break
                if (row[0] == line[0] and 
                    ((float(line[1]) <= float(row[1]) and float(line[2]) >= float(row[2])) or
                     (float(line[1]) >= float(row[1]) and float(line[2]) <= float(row[2])) or
                     (float(line[1]) >= float(row[1]) and float(line[1]) <= float(row[2]) and float(line[2]) >= float(row[2])) or
                     (float(line[1]) <= float(row[1]) and float(line[2]) >= float(row[1]) and float(line[2]) <= float(row[2])))):
                    
                    labels.append((line))
                    go_label =1

                    # Get the next color from the cycle
                    color = next(color_cycle)  # Default to 'purple' if color_cycle is exhausted
                    
                    # Plot solid bar
                    #plt.plot([float(line[1]), float(line[2])], [-0.1, -0.1], color=color, linewidth=1)
                    lines = [[(float(line[1]), -0.1), (float(line[2]), -0.1)]]  # list of line segments
                    lc = LineCollection(lines, linewidths=10, colors=color)  # control linewidth here
                    plt.gca().add_collection(lc)
                    line_break =1
                    

    if(gstats and go_label==1):
        # To store unique rows

        seen_rows = set()

        for row in labels:
            row_tuple = tuple(row)  # Convert list to tuple to make it hashable
            if row_tuple not in seen_rows:
                seen_rows.add(row_tuple)
                unique_rows.append(row)

                    
    plt.xlim(float(converted_array[0][1]), float(last))
    plt.show()


def main():
    ##plotting time

    parser = argparse.ArgumentParser(description=('''\
    CNPI Plotting:
    Creating figures based on quickmer copy number data.
    
    Plotting can output abnormal copy number ranges outside of the 1.5-2.5 threshold. A default amount of 60 regions is used.
    The user can specify amount of quickmer windows to include within each ababnormal Copy Number range.
    A default buffer of 3 is used.
    Amount of windows in a row falling withing normal range before being abnormal again. Helpful for copy number oscilations.
    
    Plotting is also for plotting user specified copy number ranges. For plotting entire chomosomes or regions withing a chromosome.
    
    Can include Genotype Statistics produced from CNPI.cpp to complement copy number images.
    '''))
    parser.add_argument('-f', '--file', type=str, required=True, help="Required: An Individual's Quickmer Copy Number Data File: Required")
    parser.add_argument('-p', '--pair', nargs="*", type=str, required=False, metavar=('file1', 'file2'), help='Optional: For Plotting Duos and Trios. One or Two separate files: Files need to be the same length when plotting a trios. Use sort and filter commands to get rid of unnecessary lines')
    parser.add_argument('-r', '--reference', type=str, required=True, help='Path to the reference file: Required')
    parser.add_argument('-w', '--minWindow', type=str, required=False, default=60,help='<Min Window size of of abnormal copy number regions: Optional: Default 60')
    parser.add_argument('-i', '--windowBuffer', type=str, required=False, default=3, help='Amount of windows that can fall within normal range when plotting abnormal copy number ranges 1.5-2.5 rule: Optional: Default 3')
    parser.add_argument('-se', '--selectChrm', type=str, required=False, help='Chromosome to plot. For user specificed copy number images: Optional')
    parser.add_argument('-start', '--startPos', type=str, required=False, help='Start Position of user specified range. Must also include Chromosome to plot: Optional')
    parser.add_argument('-stop', '--stopPos', type=str, required=False, help='Stop Position of user specified range. Must also include chromosome to plot: Optional')
    parser.add_argument('-gstat', '--gstat_txt', type=str, required=False, help='For complementing plots with Genotype_Stats.txt: Optional for seeing region statistics')
    parser.add_argument('-o','--output', type=str, required=False, help='Output Image Label')
    parser.add_argument('-min','--minimum', type=float, required=False, default=1.5, help='Output Image Label')
    parser.add_argument('-max','--maximum', type=float, required=False, default=2.5, help='Output Image Label')
    #parser.add_argument('-b','--basename',type=str,required=False,help='basname')



    args = parser.parse_args()

    skip = True
    pickChrm = False
    tryStart = False
    tryStop = False
    gstats = False
    parent_len =0
    outputs = False
    output = ""
    

    if args.pair is not None and (len(args.pair) == 2 or len(args.pair) ==1):
        if (len(args.pair)==2):
            gz_parent1 = args.pair[0]
            gz_parent2 = args.pair[1]
            parent_len =2
        elif (len(args.pair)==1):
            gz_parent1 = args.pair[0]
            parent_len = 1

    else:
        skip = False

    if args.output is not None:
        outputs = True
        output=args.output


    if args.selectChrm is not None:
        pickChrm=True
        select_chrm = args.selectChrm

    if args.startPos is not None:
        tryStart=True
        start_pos = int(args.startPos)

    if args.stopPos is not None:
        tryStop=True
        stop_pos = int(args.stopPos)

    if args.gstat_txt is not None:
        gstats=True
        gstat_path = args.gstat_txt

    del_cn = args.minimum

    dup_cn = args.maximum
    print(f'The duplication cn is {dup_cn} and the deletion cn is {del_cn}')

        
    # Access the file paths provided
    gz_path = args.file
    if(skip and len(args.pair)==2):
        gz_parent1 = args.pair[0]
        gz_parent2 = args.pair[1]
    elif(skip and len(args.pair)==1):
        gz_parent1 = args.pair[0]
        
    file_path = args.reference
    pass_limit = int(args.minWindow)
    fail_limit = int(args.windowBuffer)

    
    # Process the files
    print(f'Person of Interest: {gz_path}')
    if(skip):
        if(parent_len==1):
            print(f'Parent 1: {gz_parent1}')
        elif(parent_len==2):
            print(f'Parent 1: {gz_parent1}')
            print(f'Parent 2: {gz_parent2}')
    print(f'Reference file: {file_path}')
    print(f'Pass Limit: {pass_limit}')
    print(f'Fail Limit: {fail_limit}')



    gz_array = []
    kary_array = []
    gz_parent1_array = []
    gz_parent2_array = []
    gstat_array = []
 

    with open(file_path, 'r') as file:
        # Read the entire contents of the file
        for line in file:
            # Process each line here
            kary_fields = line.strip().split('\t')
            kary_array.append(kary_fields)


    #gz file for plotting
    with gzip.open(gz_path, 'rt') as gz_file:
        for line in gz_file:
            # Process each line here
            fields = line.strip().split('\t')
            gz_array.append(fields)

    if(skip and (parent_len==1 or parent_len==2)):
        #gz file for plotting
        with gzip.open(gz_parent1, 'rt') as gz_file:
            for line in gz_file:
                # Process each line here
                fields = line.strip().split('\t')
                gz_parent1_array.append(fields)
    if(skip and (parent_len==2)):
        #gz file for plotting
        with gzip.open(gz_parent2, 'rt') as gz_file:
            for line in gz_file:
                # Process each line here
                fields = line.strip().split('\t')
                gz_parent2_array.append(fields)

    if(gstats):
        with open(gstat_path, 'r') as stats:
        # Read the entire contents of the file
            for line in stats:
                # Process each line here
                fields = line.strip().split('\t')
                gstat_array.append(fields)

    
    combined_data_sex = {
        'Chr': [row[0] for row in gz_array],
        'Lowest_Start_Pos': [int(row[1]) for row in gz_array],
        'Highest_End_Pos': [int(row[2]) for row in gz_array],
        'CN': [float(row[3]) for row in gz_array],
    }

    if(skip and (parent_len==2)):
        combined_data_sex2 = {
            'Parent1_Chr': [row[0] for row in gz_parent1_array],
            'Parent1_Lowest_Start_Pos': [int(row[1]) for row in gz_parent1_array],
            'Parent1_Highest_End_Pos': [int(row[2]) for row in gz_parent1_array],
            'Parent1_CN': [float(row[3]) for row in gz_parent1_array],
            'Parent2_Chr': [row[0] for row in gz_parent2_array],
            'Parent2_Lowest_Start_Pos': [int(row[1]) for row in gz_parent2_array],
            'Parent2_Highest_End_Pos': [int(row[2]) for row in gz_parent2_array],
            'Parent2_CN': [float(row[3]) for row in gz_parent2_array],
        }
        combined_data_sex.update(combined_data_sex2)

    if(skip and (parent_len==1)):
        combined_data_sex2 = {
            'Parent1_Chr': [row[0] for row in gz_parent1_array],
            'Parent1_Lowest_Start_Pos': [int(row[1]) for row in gz_parent1_array],
            'Parent1_Highest_End_Pos': [int(row[2]) for row in gz_parent1_array],
            'Parent1_CN': [float(row[3]) for row in gz_parent1_array],
        }
        combined_data_sex.update(combined_data_sex2)

    if(pickChrm):
        gz_array = [row for row in gz_array if row[0] == select_chrm]
        gz_parent1_array = [row for row in gz_parent1_array if row[0] == select_chrm]
        gz_parent2_array = [row for row in gz_parent2_array if row[0] == select_chrm]

    extend = 50000
    if tryStart and pickChrm and tryStop:
        extend_start_pos = start_pos
        if(start_pos - extend>=0):
            extend_start_pos= start_pos - extend
        else:
            extend_start_pos=0
        extend_stop_pos = stop_pos + 50000

    if pickChrm and tryStart and tryStop:
        gz_array = [row for row in gz_array if row[0] == select_chrm and int(row[2]) > extend_start_pos and int(row[1]) < extend_stop_pos]
        gz_parent1_array = [row for row in gz_parent1_array if row[0] == select_chrm and int(row[2]) > extend_start_pos and int(row[1]) < extend_stop_pos]
        gz_parent2_array = [row for row in gz_parent2_array if row[0] == select_chrm and int(row[2]) > extend_start_pos and int(row[1]) < extend_stop_pos]

    # Combine all datasets into one for plotting
    combined_data = {
        'Chr': [row[0] for row in gz_array],
        'Lowest_Start_Pos': [int(row[1]) for row in gz_array],
        'Highest_End_Pos': [int(row[2]) for row in gz_array],
        'CN': [float(row[3]) for row in gz_array],
    }


    if(skip and (parent_len==2)):
        combined_data2 = {
            'Parent1_Chr': [row[0] for row in gz_parent1_array],
            'Parent1_Lowest_Start_Pos': [int(row[1]) for row in gz_parent1_array],
            'Parent1_Highest_End_Pos': [int(row[2]) for row in gz_parent1_array],
            'Parent1_CN': [float(row[3]) for row in gz_parent1_array],
            'Parent2_Chr': [row[0] for row in gz_parent2_array],
            'Parent2_Lowest_Start_Pos': [int(row[1]) for row in gz_parent2_array],
            'Parent2_Highest_End_Pos': [int(row[2]) for row in gz_parent2_array],
            'Parent2_CN': [float(row[3]) for row in gz_parent2_array],
        }
        combined_data.update(combined_data2)

    if(skip and (parent_len==1)):
        combined_data2 = {
            'Parent1_Chr': [row[0] for row in gz_parent1_array],
            'Parent1_Lowest_Start_Pos': [int(row[1]) for row in gz_parent1_array],
            'Parent1_Highest_End_Pos': [int(row[2]) for row in gz_parent1_array],
            'Parent1_CN': [float(row[3]) for row in gz_parent1_array],
        }
        combined_data.update(combined_data2)

    # Create a DataFrame for combined data
    df_combined = pd.DataFrame(combined_data)
    df_combined_sex = pd.DataFrame(combined_data_sex)
    
    
    parent1_type = "Parent 1"
    parent2_type = "Parent 2"
    if(skip and pickChrm is None and tryStart is None and tryStop is None):
        average_cn1 = df_combined[df_combined['Parent1_Chr'] == 'chrY']['Parent1_CN'].mean()
        if average_cn1<.7:
            parent1_type = "Mother"
            print(f'{gz_parent1} is the Mother')
        else:
            parent1_type = "Father"
            print(f'{gz_parent1} is the Father')
        if parent_len !=1:
            average_cn2 = df_combined[df_combined['Parent2_Chr'] == 'chrY']['Parent2_CN'].mean()
            if average_cn1>.7:
                parent2_type = "Mother"
                print(f'{gz_parent2} is the Mother')
            else:
                parent2_type = "Father"
                print(f'{gz_parent2} is the Father')
        

    if(skip is not None and pickChrm is not None and tryStart is not None and tryStop is not None):
        average_cn1 = df_combined_sex[df_combined_sex['Parent1_Chr'] == 'chrY']['Parent1_CN'].mean()
        if average_cn1<.7:
            parent1_type = "Mother"
            print(f'{gz_parent1} is the Mother')
        else:
            parent1_type = "Father"
            print(f'{gz_parent1} is the Father')
        if parent_len !=1:
            average_cn2 = df_combined_sex[df_combined_sex['Parent2_Chr'] == 'chrY']['Parent2_CN'].mean()
            if average_cn1>.7:
                parent2_type = "Mother"
                print(f'{gz_parent2} is the Mother')
            else:
                parent2_type = "Father"
                print(f'{gz_parent2} is the Father')
        

    if (pickChrm):
        go_plot(df_combined, kary_array, parent1_type, parent2_type, gstats, gstat_array, select_chrm, skip, parent_len, start_pos, stop_pos, outputs, output, del_cn, dup_cn)

    consecutive_count = 0
    consecutive_regions = []
    unique_rows = []
    consecutive_failures = 0
    previous_row = None
    previous_chrm = ""
    pre_array = []
    post_array = []
    stop =0

    pre_array = deque(maxlen=30)

    file= 'consecutive_fails.txt'
    with open(file, "w") as f:

        f.write(gz_path +"\t"+ kary_array[0][0])

        for i, row in enumerate(gz_array):

            if(previous_row is not None and row[0]!=previous_row[0] and i>0):
                pre_array.clear()
                consecutive_count =0
                consecutive_failures =0
            #elif(previous_row is not None and len(pre_array)>1):
            #    pre_array.clear()

            if (float(row[3]) > dup_cn or float(row[3]) < del_cn) and (row[0] != "chrX" and row[0] != "chrY"):
                if (previous_row is not None and consecutive_failures >0 and consecutive_count !=0 and row[0]==previous_chrm):
                    consecutive_regions.append('\t'.join(previous_row))
                if(row[0]==previous_chrm):
                    consecutive_count += 1
                    consecutive_failures =0
                    consecutive_regions.append('\t'.join(row))

            elif ((kary_array[0][0] == "46,XY") and (float(row[3]) > dup_cn-1 or float(row[3]) < del_cn-1)) and (row[0] == "chrX" or row[0] == "chrY"):
                if (previous_row is not None and consecutive_failures >0 and consecutive_count !=0 and row[0]==previous_chrm):
                    consecutive_regions.append('\t'.join(previous_row))
                if(row[0]==previous_chrm):
                    consecutive_count += 1
                    consecutive_failures=0
                    consecutive_regions.append('\t'.join(row))

            elif ((kary_array[0][0] == "46,XX") and (float(row[3]) > dup_cn or float(row[3]) < del_cn)) and (row[0] == "chrX"):
                            if(previous_row is not None and consecutive_failures >0 and consecutive_count !=0 and row[0]==previous_chrm):
                                consecutive_regions.append('\t'.join(previous_row))
                            if(row[0]==previous_chrm):
                                consecutive_count += 1
                                consecutive_failures=0
                                consecutive_regions.append('\t'.join(row))

            else:
                consecutive_failures +=1
                pre_array.append('\t'.join(row))
                if(row[0]==previous_chrm and stop==0):
                    #pre_array.append('\t'.join(row))
                    #pre_array.append(row)

                    consecutive_regions.append('\t'.join(row))
                    stop=1
                
                # Reset consecutive count and save/print if 20 or more consecutive
                if consecutive_failures > fail_limit:
                    if consecutive_count >= pass_limit:
                        f.write("\nConsecutive regions:")
                        
                        b =i
                        while len(post_array) < 20 and b + 1 < len(gz_array):
                            post_array.append(gz_array[b + 1])
                            b+=1

                        if(row[0]==previous_chrm):
                            print(f'The consecutive count is {consecutive_count} the pass limit is {pass_limit} the consecuitive failures is {consecutive_failures} the fail limit is {fail_limit}')
                            plot_function(consecutive_regions, df_combined, skip, kary_array[0][0], gstats, gstat_array, unique_rows, parent1_type, parent2_type, parent_len, pre_array, post_array, del_cn, dup_cn)

                        pre_array.clear()
                        post_array.clear()

                        for region in consecutive_regions:
                            f.write("\n"+region)

                        if(gstats):
                            f.write("\n\nTranscripts within this region:")
                            for row in unique_rows:
                                f.write("\n"+ "\t".join(row))
                        f.write("\n------------------------------")

                    consecutive_count = 0
                    consecutive_regions = []
                    unique_rows =[]
                    consecutive_failures =0
            previous_row = row
            previous_chrm = row[0]

    f.close() 

    print(kary_array[0][0])


main()
