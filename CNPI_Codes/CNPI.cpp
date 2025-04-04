#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <vector>
#include <numeric>
#include <iomanip>
#include <getopt.h>
#include <cmath>
#include <algorithm>
#include <thread>
#include <sstream>
#include <cstring>
#include <filesystem>



// Function to round a number to a specified number of decimal places
double roundToDecimals(double value, int decimals) {
    double factor = std::pow(10.0, decimals);
    return std::round(value * factor) / factor;
}

// Define long options
//this is for listing arguments with a longer name rather than just having a -letter to do so
static struct option long_options[] = {
    {"gz_path", required_argument, NULL, 'g'},
    {"bed_path_rows", required_argument, NULL, 'n'},
    {"output_name", required_argument, NULL, 'o'},
    {"record", required_argument, NULL, 'r'},
    {"bed_gz_path", required_argument, NULL, 'd'},
    {"help", no_argument, NULL, 'h'},
    {"average_check", required_argument, NULL, 's'},
    {"sd_check", required_argument, NULL, 't'},
    {"chromosome", required_argument, NULL, 'e'},
    {"start_stop", required_argument, NULL, 'p'},
    {"cn_rows", required_argument, NULL, 'c'},
    {"duplication", required_argument, NULL, 'u'},
    {"deletion", required_argument, NULL, 'l'},
    {"weighted", required_argument, NULL, 'w'},
    {0, 0, 0, 0}
};

//finds the karyotype across regions within the bed file
//used the genotype file as input
void regions(std::string **genotype_table_output, std::string **stats_tabs_table, int bed_path_rows, int unique_headers, std::string output_name, bool chr, 
bool start_stop, int start, int end, std::string return_chromosome, std::string* headers){
    
    //for finding the lowest and highest value included within each karyotype
    int high_table[unique_headers];
    int low_table[unique_headers];
    int low = 999999999;
    int chrm =0;

    for(int i=0; i<bed_path_rows-1; i++){
        if(stats_tabs_table[i][1]==stats_tabs_table[i+1][1]){
            if(std::stoi(stats_tabs_table[i][3])>= std::stoi(stats_tabs_table[i+1][3])){
                high_table[chrm]=i;
            }else if(std::stoi(stats_tabs_table[i][3])< std::stoi(stats_tabs_table[i+1][3])){
                high_table[chrm]=i+1;
            }if(std::stoi(stats_tabs_table[i][2]) <=low){
                low = std::stoi(stats_tabs_table[i][2]);
                low_table[chrm] =i;
            }
            
        //iterating to next chromosome    
        }else if(stats_tabs_table[i][1] != stats_tabs_table[i+1][1]){
            headers[chrm] = stats_tabs_table[i][1];
            chrm +=1;
            low = 999999999;
            
            if(chrm>=unique_headers){
                break;
            }
        }
    }

    //outputting certan chromsome to screen need -regions_karyotype -chromsome chr10
    if(chr==true && start_stop==false){
        bool valid_chr = false;
        std::cout<<"\n\nChr\tStart\tStop\tRegion\tCN_Average\tWeighted_Avg\tCN_SD\tCN_CV\tRegions_Size\tTotal_Windows\tGen_Variation";
        for(int i=0; i<bed_path_rows; i++){
            if(return_chromosome ==genotype_table_output[i][1]){
                std::cout<< "\n" << genotype_table_output[i][1] << "\t" << genotype_table_output[i][2] << "\t"<< genotype_table_output[i][3] << "\t"
                << stats_tabs_table[i][0] << "\t"<< genotype_table_output[i][4] << "\t"<< genotype_table_output[i][8] << "\t"<< genotype_table_output[i][5]
                << "\t"<< genotype_table_output[i][9] << "\t" << genotype_table_output[i][6] << "\t" << genotype_table_output[i][7] << "\t" << genotype_table_output[i][10];
                valid_chr = true;
            }
        }
        if(valid_chr==false){
            std::cout<< "\n\nYou Did Not Enter a Valid Chromosome! The Following are valid Chromosomes:\n";
            for(int i=0; i<unique_headers; i++){
                std::cout << headers[i] << " ";
            }
        }    
    }

    //ouputting certain range of a chromsome need -regions_karyotype -chromsome chr10 -start_stop 1000 1000000
    if(chr && start_stop){
        bool valid_range = false;
        std::cout<<"\n\nChr\tStart\tStop\tRegion\tCN_Average\tWeighted_Avg\tCN_SD\tCN_CV\tRegions_Size\tTotal_Windows\tGen_Variation";
        for(int i=0; i<bed_path_rows; i++){
            if(return_chromosome ==stats_tabs_table[i][1] && (std::stoi(stats_tabs_table[i][2]) >= start) && (std::stoi(stats_tabs_table[i][3])<= end)){
                std::cout<< "\n" << genotype_table_output[i][1] << "\t" << genotype_table_output[i][2] << "\t"<< genotype_table_output[i][3] << "\t"
                << stats_tabs_table[i][0] << "\t"<< genotype_table_output[i][4] << "\t"<< genotype_table_output[i][8] << "\t"<< genotype_table_output[i][5]
                << "\t"<< genotype_table_output[i][9] << "\t" << genotype_table_output[i][6] << "\t" << genotype_table_output[i][7] << "\t" << genotype_table_output[i][10];  
                valid_range = true;
            }
        }
        if(valid_range==false){
            std::cout << "\n\nYour Chromosome or Chromosme range was incorrect! The following are valid Chromosomes and regions ranges:";
            for(int i=0; i<unique_headers; i++){
                std::cout << "\n" << headers[i] << "\t" << stats_tabs_table[low_table[i]][2] << "\t" << stats_tabs_table[high_table[i]][3];
            }
        }
    }

    std::cout<< "\n\nRegion Output complete!";
}


//find regions of the bed file that have standard deviations above a certain value
void count_sd(float score_check, std::string output_name, int bed_path_rows, std::string **genotype_table_output , std::string **stats_tabs_table, bool output_provide, std::string filename){

    std::cout<< "\n\nSD Count Abnormalities:"; 
    std::string file;

    if(output_provide){
        file = output_name + "_SD_Abnormality_Output.txt";
    }else{
        file = filename + "_SD_Abnormality_Output.txt";
    }

    
    std::ofstream f(file);

    f<< "Chr\tStart\tStop\tRegion\tCN_Average\tWeighted_Avg\tCN_SD\tCN_CV\tRegions_Size\tTotal_Windows\tGen_Variation";

    //write to file if sd from genotype file is above threshold and is not equal to nan
    for (int i = 0; i < bed_path_rows - 1; i++) {
        if (genotype_table_output[i][0] == "chrY") {
            if (std::stof(genotype_table_output[i][5]) >= (score_check) && (genotype_table_output[i][5]!= "nan")) {
                f << "\n" << genotype_table_output[i][1] << "\t"<< genotype_table_output[i][2] << "\t" << genotype_table_output[i][3]
                 << "\t"<< stats_tabs_table[i][0] << "\t"<< genotype_table_output[i][4] << "\t"<< genotype_table_output[i][8] << "\t"<< genotype_table_output[i][5]
                 << "\t"<< genotype_table_output[i][9] << "\t" << genotype_table_output[i][6] << "\t" << genotype_table_output[i][7] << "\t" << genotype_table_output[i][10];
            }
        } else {
            if (std::stof(genotype_table_output[i][5]) >= (score_check) && (genotype_table_output[i][5]!= "nan")) {
                f << "\n" << genotype_table_output[i][1] << "\t"<< genotype_table_output[i][2] << "\t" << genotype_table_output[i][3]
                 <<  "\t"<< stats_tabs_table[i][0]<< "\t"<< genotype_table_output[i][4] << "\t"<< genotype_table_output[i][8] << "\t"<< genotype_table_output[i][5]
                 << "\t"<< genotype_table_output[i][9] << "\t" << genotype_table_output[i][6] << "\t" << genotype_table_output[i][7] << "\t" << genotype_table_output[i][10];
            }
        }
    }

    f.close();
    std::cout << "\nAbnormalities based off of copy number standard deviation complete! The Score Check is " << score_check;

}


//find count abnormalities if the average of the region is above or below 2 or 1 by a selected value
void count_abnormalities(std::string **genotype_table_output, float score_check, std::string output_name, int bed_path_rows, std::string **stats_tabs_table, bool output_provide, std::string filename){
    
    std::cout<< "\n\nAverage Count Abnormalities:";

    //write it to a file
    std::string file;

    if(output_provide){
        file = output_name + "_Abnormality_Output.txt";
    }else{
        file = filename + "_Abnormality_Output.txt";
    }
    
    
    std::ofstream f(file);

    f<< "Chr\tStart\tStop\tRegions\tCN_Average\tWeighted_Avg\tCN_SD\tCN_CV\tRegions_Size\tTotal_Windows\tGen_Variation";

    int female=1;
    int counts =0;
    float total =0;

    for(int i=0; i<bed_path_rows; i++){
        if(genotype_table_output[i][0]=="chrY"){
            total += std::stof(genotype_table_output[i][4]);
            counts +=1;
        }

    }
    if(total/counts < 1.3){female=0;}

    //based upon if if it is a male or female do the following
    for (int i = 0; i < bed_path_rows - 1; i++) {
        if (genotype_table_output[i][0] == "chrY") {
            if(std::stof(genotype_table_output[i][4]) <.05){break;}
            female=0;
            if (std::stof(genotype_table_output[i][4]) >= (1 + score_check) || std::stof(genotype_table_output[i][4]) <= (1 - score_check)) {
                f << "\n" << genotype_table_output[i][1] << "\t"<< genotype_table_output[i][2] << "\t" << genotype_table_output[i][3]
                 <<  "\t"<< stats_tabs_table[i][0] << "\t"<< genotype_table_output[i][4] << "\t"<< genotype_table_output[i][8] << "\t"<< genotype_table_output[i][5]
                 << "\t"<< genotype_table_output[i][9] << "\t" << genotype_table_output[i][6] << "\t" << genotype_table_output[i][7] << "\t" << genotype_table_output[i][10];
            }
        }else if(genotype_table_output[i][0]== "chrX" && female==1){
            if (std::stof(genotype_table_output[i][4]) >= (2 + score_check) || std::stof(genotype_table_output[i][4]) <= (2 - score_check)) {
                f << "\n" << genotype_table_output[i][1] << "\t"<< genotype_table_output[i][2] << "\t" << genotype_table_output[i][3]
                 <<  "\t"<< stats_tabs_table[i][0] << "\t"<< genotype_table_output[i][4] << "\t"<< genotype_table_output[i][8] << "\t"<< genotype_table_output[i][5]
                 << "\t"<< genotype_table_output[i][9] << "\t" << genotype_table_output[i][6] << "\t" << genotype_table_output[i][7] << "\t" << genotype_table_output[i][10];
            }
        }else if(genotype_table_output[i][0]== "chrX" && female==0){
            if (std::stof(genotype_table_output[i][4]) >= (1 + score_check) || std::stof(genotype_table_output[i][4]) <= (1 - score_check)) {
                f << "\n" << genotype_table_output[i][1] << "\t"<< genotype_table_output[i][2] << "\t" << genotype_table_output[i][3]
                 << "\t"<< stats_tabs_table[i][0] << "\t"<< genotype_table_output[i][4] << "\t"<< genotype_table_output[i][8] << "\t"<< genotype_table_output[i][5]
                 << "\t"<< genotype_table_output[i][9] << "\t" << genotype_table_output[i][6] << "\t" << genotype_table_output[i][7] << "\t" << genotype_table_output[i][10];
            }
        }else{
            if (std::stof(genotype_table_output[i][4]) >= (2 + score_check) || std::stof(genotype_table_output[i][4]) <= (2 - score_check)) {
                f << "\n" << genotype_table_output[i][1] << "\t"<< genotype_table_output[i][2] << "\t" << genotype_table_output[i][3]
                 << "\t"<< stats_tabs_table[i][0] << "\t"<< genotype_table_output[i][4] << "\t"<< genotype_table_output[i][8] << "\t"<< genotype_table_output[i][5]
                 << "\t"<< genotype_table_output[i][9] << "\t" << genotype_table_output[i][6] << "\t" << genotype_table_output[i][7] << "\t" << genotype_table_output[i][10];
            }
        }
    }

    f.close();
    std::cout << "\nAbnormalities based off copy number average complete! The Score Check is " << score_check;

}


//standard deviation of numbers in a vector
double get_standard_dev(std::vector<double> &stdnums, float average){

    double accum = 0.0;
    std::for_each(std::begin(stdnums), std::end(stdnums), [&](const double d) {
    accum += (d - average) * (d - average);
    });
    float stdev=0;
    return sqrt(accum / (stdnums.size()-1));
    
}

void gz_stats(int &duplicated_autosome, std::string &type, int cn_rows, int unique_headers, std::string* cn_header_table, float **cn_float_tabs_table, std::string gz_path,
std::string output_name, bool chr, bool start_stop, int start, int end, std::string return_chromosome, std::string* headers, int* high_table, int* low_table, bool output_provide, std::string filename){

    
    //std::string headers[unique_headers];
    //std::cout<< "\nUnique Headers: " << unique_headers;

    float **karyotype_table = new float *[unique_headers];
    for (int i = 0; i < unique_headers; i++) {
        karyotype_table[i] = new float[7];
    }

    for(int i=0; i<unique_headers; i++){
        for(int j=0; j<7; j++){
            karyotype_table[i][j] =0.0;
        }
    }

    int karyotype_totals[unique_headers];
    for(int i=0; i<unique_headers; i++){
        karyotype_totals[i] =0;
    }

    int low = 999999999;
    int chrm =0;

    int length=0;
    std::vector<double> stdnums(length,0);  //finding standard deviation of chromosomes

    int b =0;
    int A = start;
    int B = end;
    float region_total =0;
    int regions =0;
    //for finding highest,lowest window, size of windows per chromsome, average, total windows
    for(int i=0; i<cn_rows; i++){
        if(cn_header_table[i]==cn_header_table[i+1]){
            if((cn_float_tabs_table[i][1])>= (cn_float_tabs_table[i+1][1])){
                karyotype_table[chrm][1] = (cn_float_tabs_table[i][1]);
                high_table[chrm]=i;
            }else if((cn_float_tabs_table[i][1])< (cn_float_tabs_table[i+1][1])){
                karyotype_table[chrm][1] = (cn_float_tabs_table[i+1][1]);
                high_table[chrm]=i+1;
            }if((cn_float_tabs_table[i][0]) <=low){
                karyotype_table[chrm][0] = (cn_float_tabs_table[i][0]);
                low = (cn_float_tabs_table[i][0]);
                low_table[chrm] =i;
            }

            
            if(chr & start_stop & cn_header_table[i]== return_chromosome){
                //Check if A is between C and D
                int C= cn_float_tabs_table[i][0];
                int D= cn_float_tabs_table[i][1];
                if ((A >= C && A <= D)||(B >= C && B <= D)||(A<=C && B>=D)||(A>=C && B<=D)){
                    region_total += cn_float_tabs_table[i][2];
                    regions += 1;
                }
            }

            stdnums.push_back(cn_float_tabs_table[i][2]);

            karyotype_table[chrm][2] += (cn_float_tabs_table[i][2]);
            karyotype_table[chrm][3] += 1;
            karyotype_table[chrm][5] += cn_float_tabs_table[i][2] * (cn_float_tabs_table[i][1] - cn_float_tabs_table[i][0]);
            karyotype_table[chrm][6] += (cn_float_tabs_table[i][1] - cn_float_tabs_table[i][0]);
            karyotype_totals[chrm] = karyotype_totals[chrm]+1;

        }else if(cn_header_table[i] != cn_header_table[i+1]){
            karyotype_table[chrm][2] += (cn_float_tabs_table[i][2]);
            karyotype_table[chrm][3] += 1;
            karyotype_totals[chrm] = karyotype_totals[chrm]+1;


            stdnums.push_back(cn_float_tabs_table[i][2]);
            float average = karyotype_table[chrm][2]/karyotype_totals[chrm];
            //float average = karyotype_table[chrm][2]/karyotype_table[chrm][3];
            karyotype_table[chrm][4] = get_standard_dev(stdnums, average);
            length=0;
            std::vector<double> stdnums(length,0);

            headers[chrm] = cn_header_table[i];
            chrm +=1;
            low = 999999999;
            
            if(chrm>=unique_headers){
                break;
            }
            b=0;
        }
        b++;
    }
    


    if(start_stop & chr){
        std::cout << "\n\nThe Average Copy number across windows of the specified region is:\t" << region_total/regions;
    }    
    

    int chrm_num =0;
    std::string label;
    bool trisomy = false;
    bool xxy_possibility = false; //CHECK this with actual data
    float xx_cn =0.0, x_cn =0.0;

    //is 46xy of 46xx etc
    for (int i=0; i<unique_headers; i++){
        if(headers[i] != "chrX" && headers[i] != "chrY"){
            if(karyotype_table[i][2]/karyotype_table[i][3] > 1.3 && karyotype_table[i][2]/karyotype_table[i][3] < 2.7){
                chrm_num += 2;
            }else if(karyotype_table[i][2]/karyotype_table[i][3] > 2.7){
                chrm_num += 3;
                duplicated_autosome = i+1;
                trisomy =true;
            }    
        }else if(headers[i] == "chrX"){
            if(karyotype_table[i][2]/karyotype_table[i][3] > .3 && karyotype_table[i][2]/karyotype_table[i][3] < 1.7){
                label+= "X";
                chrm_num +=1;
                x_cn = karyotype_table[i][2]/karyotype_table[i][3];
            }else if(karyotype_table[i][2]/karyotype_table[i][3] >1.7 && karyotype_table[i][2]/karyotype_table[i][3]<2.5){
                label += "XX";
                chrm_num +=2;
                xxy_possibility = true;
                xx_cn = karyotype_table[i][2]/karyotype_table[i][3];
            }else if(karyotype_table[i][2]/karyotype_table[i][3] >2.5 && karyotype_table[i][2]/karyotype_table[i][3]<3.5){
                label += "XXX";
                chrm_num +=3;
                std::cout << "\n\nXXX Karyotype the copy number is:\t" << karyotype_table[i][2]/karyotype_table[i][3];
            }else if(karyotype_table[i][2]/karyotype_table[i][3] >3.5 && karyotype_table[i][2]/karyotype_table[i][3]<4.5){
                label += "XXXX";
                chrm_num +=4;
                std::cout << "\n\nXXXX Karyotype the copy number is:\t" << karyotype_table[i][2]/karyotype_table[i][3];
            }
        }else if(headers[i] == "chrY"){
            if(karyotype_table[i][2]/karyotype_table[i][3] > .3 && karyotype_table[i][2]/karyotype_table[i][3] < 1.7){
                label+= "Y";
                chrm_num +=1;
                if(xxy_possibility==true){std::cout << "\n\nThe XX Copy Number is:\t" << xx_cn <<"\tY Karyotype the copy number is:\t" << karyotype_table[i][2]/karyotype_table[i][3];}
            }else if(karyotype_table[i][2]/karyotype_table[i][3] >1.7){
                label += "YY";
                chrm_num +=2;
                std::cout << "\n\nX Karyotype the copy number is:\t" << x_cn << "\tThe Y copy number is:\t" <<karyotype_table[i][2]/karyotype_table[i][3];
                if(xxy_possibility==true){std::cout << "\n\nThe XX Copy Number is:\t" << xx_cn <<"\tY Karyotype the copy number is:\t" << karyotype_table[i][2]/karyotype_table[i][3];} 
            }else if(karyotype_table[i][2]/karyotype_table[i][3]<.3 && (x_cn>.3 && x_cn<1.7)){
                std::cout << "\n\nThe X Copy Number is:\t" << x_cn <<"\tY Karyotype the copy number is:\t" << karyotype_table[i][2]/karyotype_table[i][3];
            }
        }
    }    
    std::string file;

    if(output_provide){
        file = output_name + "_Karyotype.txt";
    }else{
        file = filename + "_Karyotype.txt";
    }
    
    type = label;

    std::ofstream f(file);

    f << chrm_num << "," <<label;

    std::cout<< "\n\nThe Karyotype is: " << chrm_num << ","<< label;

    f<< "\nChr  Start_Pos  End_Pos  CN_Avg  Tot_Windows SD";

    for(int i=0; i<unique_headers; i++){
        f << std::fixed << std::setprecision(3) <<"\n" << headers[i] << "\t" << int(karyotype_table[i][0])<< "\t"<< int(karyotype_table[i][1]) 
        << "\t"<<(karyotype_table[i][2]/karyotype_table[i][3]) /*<< "\t" << (karyotype_table[i][5]/karyotype_table[i][6]) */<< "\t"<< int(karyotype_table[i][3]) << "\t" << karyotype_table[i][4];
    }

    if(trisomy){
        std::cout<< "+"<<duplicated_autosome <<"\nThere is trisomy of chromosome " << duplicated_autosome << "\nCopy number for this chromosome is:\t" << karyotype_table[duplicated_autosome-1][2]/karyotype_table[duplicated_autosome-1][3];
        f<< "\nThere is trisomy of chromosome " << duplicated_autosome;

        if(duplicated_autosome==21){
            std::cout<< "\nIndividual Likely Has Down Syndrome";
            f << "\nIdividual Likely Has Down Syndrome";
        }else if(duplicated_autosome==18){
            std::cout<< "\nIndividual Likely Has Edwards Syndrome";
            f << "\nIdividual Likely Has Edwards Syndrome";
        }else if(duplicated_autosome==13){
            std::cout<< "\nIndividual Likely Has Patau Syndrome";
            f << "\nIdividual Likely Has Patau Syndrome";
        }else if(duplicated_autosome==8){
            std::cout<< "\nIndividual Likely Has Warkany Syndrome 2";
            f << "\nIdividual Likely Has Warkany Syndrome 2";
        }
    }

    std::cout<< "\n\nKaryotyping Across Copy Numbers Complete!";

}


int main(int argc, char *argv[]){

    std::cout << "\n\t\t\t--------------------------------------------------------------------"
    << "\n\t\t\t| Welcome to the Copy Number Private Investigator (CNPI) Program!! |"
    <<"\n\t\t\t--------------------------------------------------------------------\n\n";

    //initilazation of arguments
    std::string gz_path;    // Declare string for gz_path outside the if block
    int bed_path_rows = 0;  // Initialize bed_path_rows outside the if block
    std::string output_name;    // Specifies ouput prefix if provided
    std::string return_record;  // Returns a specific record to terminal if provided
    std::string return_chromosome; // Returns a specific chromosome to screen if provided
    std::string bed_gz_path;    // for reference bed file. Function can open .bed and .bed.gz
    bool record =false, abnormalities =false, sd =false, start_stop=false, chr =false, run_0 =false, bed_provide =false, cn_fill=false, output_provide=false, weighted=false;
    float score_check =0, sd_check =0;
    int start =0, end =999999999; //for start and stop locations of a chromosome if provided
    int cn_rows =0;

    float del = 1.3;
    float dup = 2.7;


    //arguments short names and argument manager
   int opt;
    while ((opt = getopt_long_only(argc, argv, "g:n:o:d:s:t:e:p:r:c:u:l:wxh", long_options, NULL)) != -1) {
        switch (opt) {
            case 'g':
                gz_path = optarg;
                break;
            case 'n':
                bed_path_rows = std::stoi(optarg);
                bed_provide =true;
                break;
            case 'o':
                output_name = optarg;
                output_provide = true;
                break;
            case 'r':
                record =true;
                return_record = optarg;
                break;
            case 'd':
                bed_gz_path = optarg;
                break;
            case 's':
                abnormalities =true;
                score_check = std::stof(optarg);
                break;
            case 't':
                sd=true;
                sd_check = std::stof(optarg);
                break;
            case 'u':
                dup = std::stof(optarg);
                std::cout << "\nThe Duplication Threshold is: "<< dup;
                break;
            case 'l':
                del = std::stof(optarg);
                std::cout << "\nThe Deletion Threshold is: "<< del;
                break;
            case 'e':
                chr = true;
                return_chromosome = optarg;
                break;
            case 'c':
                cn_fill = true;
                cn_rows = std::stoi(optarg);
            case 'w':
                weighted = true;
                break;
            case 'p':
                // First argument after -s is the start location
                start = atoi(optarg);
                start_stop =true;
                // Second argument after -s is the end location (if available)
                if (optind < argc && argv[optind][0] != '-') {
                    end = atoi(argv[optind]); // Use optind for end value
                    optind++; // Move to the next argument (end value)
                }
                break;
            case 'h':
                std::cout<< "\nThis program is designed to produce copy number estimates at user specified regions. Matching copy number estimate windows with region locations of a Reference file and providing summary information for Reference File regions." 
                << " Reference files containing chromosome, start and stop, and region information. Relevant information from the Quickmer CN Estimate file included chromosome, start and stop location of windows, and copy number for each window."
                << " Information that will be recorded into an ouptut file will be the following:\n\nChromosome|\tStart|\tStop|\tRegion|\tCN_Average|\tWeighted_Avg|\tCN_SD|\tCN_CV|\tRegion_Size|\tTotal_Windows|\tGen_Variation"
                << "\n\nChromosome:\t\tReferring to the chromosome a specific reference region comes from"
                << "\nStart and Stop:\t\tReferring to the start and stop regions withing a chromosome that a regions comes from"
                << "\nRegions:\t\tThe reference region that is being characterized. For example a Gene, Exon, or a Promoter"
                << "\nCN_Average:\t\tThe Average of all the QuicKmer-2 Windows that overlap with a reference region"
                << "\nWeighted_Average:\tThe weighted average of the QuicKmer-2 windows that overlap with a reference region. The size of each window is taken into account"
                << "\nCN_SD:\t\t\tThe standard deviation of all the windows that overlap with a reference region"
                << "\nCN_CV:\t\t\tThe coefficient of variation of all the widows that overlap with a reference region"
                << "\nRegion_Size:\t\tIn base pairs the size of a reference region"
                << "\nTotal_Windows:\t\tThe amount of QuicKmer-2 windows that overlap with a reference region"
                << "\nGen_Variation:\t\tGenetic Variation of a reference region. Possible options being Reference, Deletion, Duplication, Trisomy, NA (Y chromosome reading in females)"
                << "\n\nUsage:"
                << "\n-d or -bed_gz_path : .bed or .bed.gz file with regions to match cn windows up against (also know as the reference file). Without file header: <RefSeq_Curated.bed> required!\t-d RefSeq_Curated_6_18_24.bed"
                << "\n\t\tShould include the following: Chromosome, Start, Stop, and a region label"
                << "\n-g or -gz_path : .bed.gz or .bed file containing CN windows. Without file header: <identifier.qm2.CN.1k.bed.gz> required!\t\t\t\t\t\t\t-g HG00700.qm2.CN.1k.bed.gz"
                << "\n\t\tShould include the following: Chromosome, Start, Stop, and Copy Number Estimate"
                << "\n-n or -bed_path_rows <number of lines in bed_path_rows> Optional!\t\t\t\t\t\t\t\t\t\t\t\t\t\t-n 101172"
                << "\nIf -n provided: Parsing of the reference file is faster"
                << "\n-c or -cn_rows: For inputting the amount of QuicK-mer2 rows into the program. Optional!\t\t\t\t\t\t\t\t\t\t\t\t-c 2300000"
                << "\nIf -c provided: Parsing of the QuicKmer-2 file is faster"
                << "\n-o or -output_name <preferred name of outputfile> optional!"
                << "\nIf -o provided: output file will be the following: <preferred name of output-file>identifier.txt.\t\t\t\t\t\t\t\t\t\t-o pizza.txt"
                << "\nIf -o not provided: output file will be the following: Scanning_1_read_Output<identifier.qm2.CN.1k.bed.gz>.txt\tScanning_1_read_Outputpizza.qm2.CN.1k.bed.gz.txt"
                << "\n-r or -record returning a specific region from the refseq file to the terminal. Not Required\t\t\t\t\t\t\t\t\t\t\t-r NR_024321.1"
                << "\n-s or -average_check: Will Execute Count Abnormalities function. Input a count abnormality threshold value to check each record with\t\t\t\t\t\t-s .1"
                << "\n-t or -sd_check: Will Execute the Standard Deviation Abnormalities function. Input a standard deviation abnormality threshold to check each record with\t\t\t\t-t .3"
                << "\n-e or -chromosome: For displaying a certain chromosome to the screen.\t\t\t\t\t\t\t\t\t\t\t\t\t\t-e chr7"
                << "\n-p or -start_stop: For displaying a certain range within a chromosome. Need to input -e or -chromosome to run\t\t\t\t\t\t\t\t\t-p 10000 1000000 (input -e chr)"
                << "\n-w or -weighted: Using weighted average for Gen_Variation"
                << "\n-l or -deletion: Customizable threshold for deletion value. Default is 1.3"
                << "\n-u or -duplication: Custimizable threshold for duplication value. Default is 2.7"
                << "\n\n";
                return 0; 
        }
    }

    
    std::cout << "\n____________________________________________________________________________________________________"
    "\n____________________________________________________________________________________________________";
    std::cout << "\n\nThe files are " << bed_gz_path << " and " << gz_path;
    
    
    //next number of files are for finding the amount of rows and reading the gz counts file
    gzFile cn_table = gzopen(gz_path.c_str(), "rb");
    if (!cn_table) {
        std::cerr << "\nError opening input file!\nCheck if you are entering a valid file or a valid file format!";
        return 0;
    }

    //unzipping and finding lines of .gz file
    char buffer[16384]; // Adjust buffer size as needed
    int bytes_read;
    std::string lines;

    if (cn_fill == false){
        while ((bytes_read = gzread(cn_table, buffer, sizeof(buffer))) > 0) {
            for (int i = 0; i < bytes_read; ++i) {
                if (buffer[i] != '\n') {
                    lines += buffer[i];
                }else{
                    lines.clear(); // Clear the line for the next iteration
                    cn_rows +=1;
                }
            }
        }
        //set gz file back to 0
        (gzrewind(cn_table));
    }

    //setting up cn tables
    std::string* cn_header_table = new std::string[cn_rows]();

    //save information from cn file
    //first read in as a string and then columns 2-4 saved as floats 
    
    float**cn_float_tabs_table = new float*[cn_rows];
   	std::string **cn_tabs_table = new std::string*[cn_rows];
	for(int i=0; i<cn_rows; i++){
		cn_tabs_table[i]=new std::string[4];
        cn_float_tabs_table[i]=new float[3];
        for(int j=0; j<4; j++){
            if(j<3){
                cn_float_tabs_table[i][j]=0.0;
            }
            cn_tabs_table[i][j] = "";    
        }
	}


    //populating the gz table
    int cn_unique_headers =1; //At least 1 unique header. Checking to see how many total unique headers there are
    int cn_tab =0;          //necessary because we dont want to read tabs and only want 4 per row
    int current_row =0;     //for iterating the rows
    while ((bytes_read = gzread(cn_table, buffer, sizeof(buffer))) > 0) {
        for (int i = 0; i < bytes_read; ++i) {
            if (buffer[i] != '\n') {
                if(buffer[i]== '\t'){
                    cn_tab+=1;
                    if(cn_tab==1 && current_row>1 &&(cn_tabs_table[current_row][0] != cn_tabs_table[current_row-1][0])){
                        cn_unique_headers +=1; //we found a unique header if conditions are met
                    }
                }
                cn_tabs_table[current_row][cn_tab] +=buffer[i];
            } else {
                current_row +=1;    //iterate to next row
                cn_tab=0;           //set column to 0
            }
        }
    }


    //opening the reference bed file
    gzFile gz_bed_table = gzopen(bed_gz_path.c_str(), "rb");
    if (!cn_table) {
        std::cerr << "\nError opening input file! Check if you are entering a valid file or valid file format!";
        return 0;
    }

    //unzipping or opening refseq file and finding the size
    char buff[1024]; // Adjust buffer size as needed
    int byte_read;

    
    //finding the amount of rows in the bed ncbi refseq file
    //only ran if there has not been a provided amount of rows
    if(bed_provide==false){
        while ((byte_read = gzread(gz_bed_table, buff, sizeof(buff))) > 0) {
            for (int i = 0; i < byte_read; ++i) {
                if (buff[i] != '\n') {
                    lines += buff[i];
                }else{
                    lines.clear(); // Clear the line for the next iteration
                    bed_path_rows +=1; //rows in reference file
                }
            }
        }
        gzrewind(gz_bed_table); //set back to zero
    }

    //table containing final genotype info as strings
    std::string **stats_tabs_table = new std::string*[bed_path_rows];
    //ncbi headers table. Header of each region
    std::string* ncbi_header_table = new std::string[bed_path_rows]();

    //ncbi information calculations as floats
    float **ncbi_float_table = new float*[bed_path_rows];
    //setting up bed table
   	std::string **bed_tabs_table = new std::string*[bed_path_rows];

	for (int i=0; i<bed_path_rows; i++){
		bed_tabs_table[i]=new std::string[4];
        stats_tabs_table[i]=new std::string[8];
        ncbi_float_table[i]=new float[5];
	}

    for(int i=0; i<bed_path_rows; i++){
        for(int j=0; j<8; j++){
            if(j<4){
                bed_tabs_table[i][j]="";
                ncbi_float_table[i][j]=0.0;
            }else if(j<5){
                ncbi_float_table[i][j]=0.0;
            }    
            stats_tabs_table[i][j]="";
        }
    }
    //end setting up bed ncbi table rows and initilizations
    
    //populating the refseq table
    int tabs =0, gz_bed_rows=0, gz_gff3_rows; //tabs and rows when reading file
    int unique_headers =1;  //for finding amount of headers. At least 1 unique header
    
    while ((byte_read = gzread(gz_bed_table, buff, sizeof(buff))) > 0) {
        for (int i = 0; i < byte_read; ++i) {
            if(buff[i] == '\t') {
                tabs +=1;       //tab increased so put information at next tab
                
            }else if(buff[i]=='\n'){
                if(gz_bed_rows>1 &&(bed_tabs_table[gz_bed_rows][0] != bed_tabs_table[gz_bed_rows-1][0])){
                        unique_headers +=1;     //found additional unique header from reference bed file
                    }
                gz_bed_rows+=1; //iterate to next row
                tabs=0;         //set tabs back to zero
            }else if(tabs<4){
                bed_tabs_table[gz_bed_rows][tabs] +=buff[i];
            }
        }
    }
    
    //populating cn header and cn float table. Useful for math
    for(int i =0; i<cn_rows; i++){
        cn_header_table[i]=cn_tabs_table[i][0];
        cn_float_tabs_table[i][0] = std::stof(cn_tabs_table[i][1]);
        cn_float_tabs_table[i][1] = std::stof(cn_tabs_table[i][2]);
        cn_float_tabs_table[i][2] = std::stof(cn_tabs_table[i][3]);
    }
    for(int i=0; i<bed_path_rows; i++){
           //populating ncbi header and floats tables
        ncbi_header_table[i]=bed_tabs_table[i][0];
        ncbi_float_table[i][0] = std::stof(bed_tabs_table[i][1]);
        ncbi_float_table[i][1] = std::stof(bed_tabs_table[i][2]);
    }

    delete[] cn_tabs_table;

    //record info into the stats tabs table
    for(int i=0; i<bed_path_rows; i++){
        stats_tabs_table[i][0]= bed_tabs_table[i][3];
        stats_tabs_table[i][1]= bed_tabs_table[i][0];
        stats_tabs_table[i][2]= bed_tabs_table[i][1];
        stats_tabs_table[i][3]= bed_tabs_table[i][2];
        stats_tabs_table[i][6] = std::to_string(std::stoi(bed_tabs_table[i][2])-std::stoi(bed_tabs_table[i][1]));
    }

    delete[] bed_tabs_table;
    
    int start_stops = cn_unique_headers+1; // indexing to find start and stop locations: was previously
    //find index locations of each new chromosome
    int* cn_index_table = new int[start_stops]();
    std::string *cn_header_names_table = new std::string[cn_unique_headers]();

    //populating a table of all the cn headers and locations of headers
    //populating the unique headers into a table
    int indexing =0;
    int cn_header=0;
    for(int i =0; i< cn_rows; i++){
        if(i==0){
            cn_index_table[0] =indexing;
            cn_header_names_table[cn_header] = cn_header_table[i];
            indexing +=1;
            cn_header +=1;
        }else if(cn_header_table[i] == cn_header_table[i-1]){
            continue;
        }else{
            if(cn_header_table[i] != cn_header_table[i-1]){
                cn_index_table[indexing] = i; //new found header position
                cn_header_names_table[cn_header] = cn_header_table[i]; //new found header
                indexing+=1;
                cn_header+=1;
            }
        }
        if(indexing==start_stops-1){
            cn_index_table[indexing]= cn_rows;
            break;
        }
    }

    //displaying header amount and headers
    std::cout<< "\n\nThere are "<< cn_unique_headers << " CN Chromosomes: ";
    for(int i=0; i<cn_unique_headers; i++){
        std::cout<< " "<< cn_header_names_table[i];
    }


    //find index locations of each new chromosome
    //saving the names of all the unique headers 
    int* ncbi_bed_index_table = new int[start_stops]();
    std::string *ncbi_header_names_table = new std::string[unique_headers]();

    int ncbi_header=0;
    int index_it =0;
    for(int i =0; i< bed_path_rows; i++){
        if(i==0){
            ncbi_bed_index_table[i] =index_it;
            ncbi_header_names_table[ncbi_header] = ncbi_header_table[i];
            index_it +=1;
            ncbi_header +=1;
        }else if(ncbi_header_table[i] == ncbi_header_table[i-1]){
            continue;
        }else{
            if(ncbi_header_table[i] != ncbi_header_table[i-1]){
                ncbi_bed_index_table[index_it] = i; //new found header position
                ncbi_header_names_table[ncbi_header] = ncbi_header_table[i]; //new found header
                ncbi_header +=1;
                index_it+=1;
            }
        }
        if(index_it==unique_headers){
            ncbi_bed_index_table[index_it] = bed_path_rows;
            break;
        }
    }

    //dipslaying header amount and headers
    std::cout<< "\nThere are " << unique_headers << " Ref Chromosomes:";
    for(int i=0; i<unique_headers; i++){
        std::cout<< " "<< ncbi_header_names_table[i];
    }

    bed_path_rows = ncbi_bed_index_table[unique_headers];

    //rest of the cn data stored as floats
    int**reorder_cn_index_table = new int*[cn_unique_headers];
	for (int i=0; i<cn_unique_headers; i++){
		reorder_cn_index_table[i]=new int[2];
	}
    for(int i=0; i<cn_unique_headers; i++){
        for(int j=0; j<2; j++){
            reorder_cn_index_table[i][j]=0;
        }
    }
    
    //finding the starts and stops of the headers
    for(int i=0; i<cn_unique_headers; i++){
        reorder_cn_index_table[i][0]=cn_index_table[i];
        reorder_cn_index_table[i][1]=(cn_index_table[i+1]-cn_index_table[i]);
    }


    //if there is a varied amount of headers. If the refseq has less than the gz file match up the unique header location indexes
    bool partial = false;
    if(unique_headers!=cn_unique_headers){
        int match=0;
        std::cout<< "\n\nThe reference file has " << unique_headers << " chromosomes while the genome file has "<< cn_unique_headers;
        for(int i=0; i<unique_headers; i++){
            int reorder=0;

            for(int j=0; j<cn_unique_headers; j++){
                if(ncbi_header_names_table[i] == cn_header_names_table[j]){
                    std::cout<< "\nChromosome match " << ncbi_header_names_table[i]; //<< " and " << cn_header_names_table[j];
                    std::cout<< ":\tGenome file start position: " << reorder_cn_index_table[j][0] << "\tIt is this many windows\t" << reorder_cn_index_table[j][1];
                    reorder_cn_index_table[i][0]=reorder_cn_index_table[j][0];
                    reorder_cn_index_table[i][1]=reorder_cn_index_table[j][1];
                }
            }match++;
        }
        partial=true;
        std::cout<< "\nThere are "<< match << " matching chromosomes between the reference and genome files!";
    }


    delete[] cn_header_names_table;
    delete[] ncbi_header_names_table;

    //rest of the cn data stored as floats
    float**weighted_bed = new float*[bed_path_rows];
	for (int i=0; i<bed_path_rows; i++){
		weighted_bed[i]=new float[2];
	}

    int gz_run_index =0;    //iterating through gz index table
    int gz_run= reorder_cn_index_table[0][0];   //gz position
    int gz_end = reorder_cn_index_table[0][1]+reorder_cn_index_table[0][0]; //end of gz header
    int i=0;
    int index =0;

    float A =0, B =0, C =0, D =0, mean =0;
    int stop =0, stop2 =0, iteration=0;


    for(; i<(bed_path_rows); i++){
        int go =0;

        //if find the correct start and stop index of cn_table with respect to where it should match in the reference table
        if(i==ncbi_bed_index_table[index+1]){index+=1, gz_run=reorder_cn_index_table[index][0], gz_end = gz_run+reorder_cn_index_table[index][1], gz_run_index=gz_run;}
        

        int length=0;
        std::vector<double> stdnums(length,0);  //for holding standard deviation information
        //running thought the cn file
        for(; gz_run < (gz_end);){
           
            //checking to see if headers of the cn file and the ncbi file are the same
            if(ncbi_header_table[i] == cn_header_table[gz_run]){

                //save values
                A = cn_float_tabs_table[gz_run][0];
                B = cn_float_tabs_table[gz_run][1];
                C = ncbi_float_table[i][0];
                D = ncbi_float_table[i][1];


                //Check if A is between C and D
                if ((A >= C && A <= D)||(B >= C && B <= D)||(A<=C && B>=D)||(A>=C && B<=D) && (gz_run!=gz_end-1)){
                    if(go!=1){
                        iteration = gz_run;
                        go=1;
                    }

                    weighted_bed[i][0] = cn_float_tabs_table[gz_run][2]*(B-A);
                    weighted_bed[i][1] = (B-A);
                    ncbi_float_table[i][2] += cn_float_tabs_table[gz_run][2];
                    ncbi_float_table[i][3] += 1;
                    stdnums.push_back(cn_float_tabs_table[gz_run][2]);

                    if(gz_run == (gz_end-1)){ //for finding avearge and sd
                        ncbi_float_table[i][2] = ncbi_float_table[i][2]/ncbi_float_table[i][3];
                        float average = ncbi_float_table[i][2];   
                        ncbi_float_table[i][4]= get_standard_dev(stdnums, average);
                    }
                //if outside calculate Standard Deviation and update indexing of where to scan  
                }else if(A > C && A > D && B > C && B > D){
                    ncbi_float_table[i][2] = ncbi_float_table[i][2]/ncbi_float_table[i][3];
                    float average = ncbi_float_table[i][2];
                    ncbi_float_table[i][4]=get_standard_dev(stdnums, average);
                    if(i<bed_path_rows-1){
                        if(ncbi_float_table[i+1][0]>=ncbi_float_table[i][0]){
                            gz_run=iteration;
                        }else if(ncbi_float_table[i+1][0]<ncbi_float_table[i][0]){
                            gz_run = gz_run_index;
                        }
                    }    
                    break;
               }
            }

            //increment gz run and then test to see if there is an occaision where the ncbi file is out of order. If out of order set index back      
            gz_run++;
            if(i==bed_path_rows-1 && partial==false){
                gz_run=gz_run_index;
                break;
            }
            if(gz_run>=(gz_end)){
                //int b=1; might be bogus but maybe not
                if(ncbi_header_table[i+1]==ncbi_header_table[i]){
                    gz_run=gz_run_index;
                    iteration=gz_run_index;
                }
                break;
            }       
        }

    }

    delete[] ncbi_header_table;


    //reoder the new float data to be in the stats table
    for(int i =0; i<bed_path_rows; i++){
        double rounded_Avg = roundToDecimals(ncbi_float_table[i][2],3);
        std::ostringstream oss1;
        oss1 << std::fixed << std::setprecision(3) << rounded_Avg;
        stats_tabs_table[i][4]= oss1.str();
        double rounded_SD = roundToDecimals(ncbi_float_table[i][4],3);
        std::ostringstream oss2;
        oss2 << std::fixed << std::setprecision(3) << rounded_SD;
        stats_tabs_table[i][5]= oss2.str();
        stats_tabs_table[i][7]= std::to_string(int(ncbi_float_table[i][3])); 
    }

    delete[] ncbi_float_table;

    //karyotyping
    int duplicated_autosome =0; //if an autosome has a 3rd chromosome
    std::string type ="";       //sex chromosome value
    //cn karyotype function
    //std::cout<< "\nUnique Headers: " << cn_unique_headers;

    std::string headers[cn_unique_headers];
    int high_table[cn_unique_headers]; //highest cn value
    int low_table[cn_unique_headers];  //lowest cn value


    std::string file;
    std::string filename;

    if(output_provide==false){
        std::filesystem::path s(gz_path);
        filename = s.filename().string();
    }

    gz_stats(duplicated_autosome ,type, cn_rows, cn_unique_headers, cn_header_table, cn_float_tabs_table, gz_path, output_name,
    chr, start_stop, start, end, return_chromosome, headers, high_table, low_table, output_provide, filename);

    delete[] cn_header_table;
    delete[] cn_float_tabs_table;


    if(output_provide){
        file = output_name + "_Chromosomal_Sex.txt";
    }else{
        file = filename + "_Chromosomal_Sex.txt";
    }
    //file = output_name + "_Chromosomal_Sex.txt";
    std::ofstream s(file);
    std::cout<< "\n\nThe Chromosomal Sex is: "<< type;
    if(s.is_open()){
        if(type=="XX" or type=="X" or type=="XXX"){
            s << "Female";
            std::cout<< " - Female";
        }else if(type=="XY" or type=="XXY" or type=="XYY" or type == "XXXY" or type=="XXXXY" or type=="XXXXXY"){
            s << "Male";
            std::cout<< " - Male";
        }
    }s.close();
    
    //open file and record the stats information
    //std::string file;
    if(output_provide){
        file = output_name + "_Genotype.txt";
    }else{
        file = filename + "_Genotype.txt";
    }

    
    std::ofstream f(file);
    f << std::fixed << std::setprecision(3);

    std::string **genotype_table_output = new std::string*[bed_path_rows];
    for(i=0; i<bed_path_rows; i++){
        genotype_table_output[i]=new std::string[11];
    }


    if(f.is_open()){
        f<< "Chr\tStart\tStop\tRegion\tCN_Average\tWeighted_Avg\tCN_SD\tCN_CV\tRegion_Size\tTotal_Windows\tGen_Variation";
        
        int E = 56887902;
        int F = 57217415;
        int G = 155701382;
        int H = 156030895;
        
        for(int i=0; i<bed_path_rows; i++){
            f<< "\n";
            for(int j =1; j<8; j++){
                f<< stats_tabs_table[i][j];
                genotype_table_output[i][j] = stats_tabs_table[i][j];

                if(j==3){
                    f << "\t" << stats_tabs_table[i][0];
                }

                //labeling pseudoautosomal regions of chrX and chrY with special character *
                if(j==1 and(stats_tabs_table[i][1]=="chrX" or stats_tabs_table[i][1]=="chrY")){

                    A = std::stoi(stats_tabs_table[i][2]);
                    B = std::stoi(stats_tabs_table[i][3]);
                    C = 10000;
                    D = 2781479;

                    //Check if A is between C and D for chrX and chrY 1st pseudoautosomal region
                    if ((A >= C && A <= D)||(B >= C && B <= D)||(A<=C && B>=D)||(A>=C && B<=D)){
                        f <<"*";
                        genotype_table_output[i][j] += "*";
                    }//Check if A is between G and H for chrX 2nd pseudoautosomal region
                    else if ((stats_tabs_table[i][1]=="chrX") && ((A >= G && A <= H)||(B >= G && B <= H)||(A<=G && B>=H)||(A>=G && B<=H))){
                        f <<"*";
                        genotype_table_output[i][j] += "*";
                    }//Check if A is between E and F for chrY 2nd pseudoautosomal region
                    else if ((stats_tabs_table[i][1]=="chrY") && ((A >= E && A <= F)||(B >= E && B <= F)||(A<=E && B>=F)||(A>=E && B<=F))){
                        f <<"*";
                        genotype_table_output[i][j] += "*";
                    }
                }
                else if(j==4){
                    f<< "\t" << weighted_bed[i][0]/weighted_bed[i][1];
                    genotype_table_output[i][8] = std::to_string(weighted_bed[i][0]/weighted_bed[i][1]);
                }    
                else if(j==5){
                    f<< "\t" << 100*(std::stof(stats_tabs_table[i][5])/std::stof(stats_tabs_table[i][4])); //Copy Number Coefficient of variation
                    genotype_table_output[i][9] = std::to_string(100*(std::stof(stats_tabs_table[i][5])/std::stof(stats_tabs_table[i][4])));
                }f<< "\t";
            }

            //double del = 1.3;
            //double dup = 2.7;
            
            //following if loop is for including genetic variation and includes more loops for special cases
            float cn_avg = std::stof(stats_tabs_table[i][4]);
            float weighted_avg = std::stof(genotype_table_output[i][8]);
            
            float threshold = cn_avg;
            if(weighted==true){
                threshold = weighted_avg;
            }

            if(stats_tabs_table[i][1]=="chrX"){
                if(type=="XX" or type=="XXY" or type=="XXYY"){
                    if(threshold<del){
                        f <<"DELETION";
                        genotype_table_output[i][10] = "DELETION";
                    }else if(threshold>del && threshold<dup){
                        f <<"REFERENCE";
                        genotype_table_output[i][10] = "REFERENCE";
                    }else if(std::stof(stats_tabs_table[i][4])>dup){
                        f <<"DUPLICATION";
                        genotype_table_output[i][10] = "DUPLICATION";
                    }
                }else if(type=="XY" or type=="X"){
                    if(std::stof(stats_tabs_table[i][4])<(del-1)){
                        f <<"DELETION";
                        genotype_table_output[i][10] = "DELETION";
                    }else if(std::stof(stats_tabs_table[i][4])>(del-1) && std::stof(stats_tabs_table[i][4])<(dup-1)){
                        f <<"REFERENCE";
                        genotype_table_output[i][10] = "REFERENCE";
                    }else if(std::stof(stats_tabs_table[i][4])>(dup-1)){
                        f <<"DUPLICATION";
                        genotype_table_output[i][10] = "DUPLICATION";
                    }
                }else if(type=="X"){
                    if(std::stof(stats_tabs_table[i][4])<(del-1)){
                        f <<"DELETION";
                        genotype_table_output[i][10] = "DELETION";
                    }else if(std::stof(stats_tabs_table[i][4])>(del-1) && std::stof(stats_tabs_table[i][4])<(dup-1)){
                        f <<"REFERENCE";
                        genotype_table_output[i][10] = "REFERENCE";
                    }else if(std::stof(stats_tabs_table[i][4])>(dup-1)){
                        f <<"DUPLICATION";
                        genotype_table_output[i][10] = "DUPLICATION";
                    }
                }else if(type=="XXX" or type=="XXXY"){
                    if(std::stof(stats_tabs_table[i][4])<(del+1)){
                        f <<"DELETION";
                        genotype_table_output[i][10] = "DELETION";
                    }else if(std::stof(stats_tabs_table[i][4])>(del+1) && std::stof(stats_tabs_table[i][4])<(dup+1)){
                        f <<"REFERENCE";
                        genotype_table_output[i][10] = "REFERENCE";
                    }else if(std::stof(stats_tabs_table[i][4])>(dup+1)){
                        f <<"DUPLICATION";
                        genotype_table_output[i][10] = "DUPLICATION";
                     }
                }else if(type=="XXXX" or type=="XXXXY"){
                    if(std::stof(stats_tabs_table[i][4])<(del+2)){
                        f <<"DELETION";
                        genotype_table_output[i][10] = "DELETION";
                    }else if(std::stof(stats_tabs_table[i][4])>(del+2) && std::stof(stats_tabs_table[i][4])<(dup+2)){
                        f <<"REFERENCE";
                        genotype_table_output[i][10] = "REFERENCE";
                    }else if(std::stof(stats_tabs_table[i][4])>(dup+2)){
                        f <<"DUPLICATION";
                        genotype_table_output[i][10] = "DUPLICATION";
                     }
                }
            }else if(stats_tabs_table[i][1]=="chrY"){
                if(type=="XY" or type=="XXY" or type=="XXXY" or type=="XXXXY" or type=="XXXXXY"){
                    if(std::stof(stats_tabs_table[i][4])<(del-1)){
                        f <<"DELETION";
                        genotype_table_output[i][10] = "DELETION";
                    }else if(std::stof(stats_tabs_table[i][4])>(del-1) && std::stof(stats_tabs_table[i][4])<(dup-1)){
                        f <<"REFERENCE";
                        genotype_table_output[i][10] = "REFERENCE";
                    }else if(std::stof(stats_tabs_table[i][4])>(dup+1)){
                        f <<"DUPLICATION";
                        genotype_table_output[i][10] = "DUPLICATION";
                    }
                }else if(type=="XYY"){
                    if(std::stof(stats_tabs_table[i][4])<del){
                        f <<"DELETION";
                        genotype_table_output[i][10] = "DELETION";
                    }else if(std::stof(stats_tabs_table[i][4])>del && std::stof(stats_tabs_table[i][4])<dup){
                        f <<"REFERENCE";
                        genotype_table_output[i][10] = "REFERENCE";
                    }else if(std::stof(stats_tabs_table[i][4])>dup){
                        f <<"DUPLICATION";
                        genotype_table_output[i][10] = "DUPLICATION";
                    }
                }else if(type=="XX" or type=="X" or type=="XXX"){
                    f<< "NA";
                }
            }else if(stats_tabs_table[i][1]== "chr"+std::to_string(duplicated_autosome)){
                if(std::stof(stats_tabs_table[i][4])<(del+1)){
                    f <<"DELETION";
                    genotype_table_output[i][10] = "DELETION";
                }else if(std::stof(stats_tabs_table[i][4])>(del+1) && std::stof(stats_tabs_table[i][4])<(dup+1)){
                    f <<"TRISOMY";
                    genotype_table_output[i][10] = "TRISOMY";
                }else{
                    f <<"DUPLICATION";
                    genotype_table_output[i][10] = "DUPLICATION";
                }
            }else{
                if(threshold<del){
                    f <<"DELETION";
                    genotype_table_output[i][10] = "DELETION";
                }else if(threshold>del && threshold<dup){
                    f <<"REFERENCE";
                    genotype_table_output[i][10] = "REFERENCE";
                }else{
                    f <<"DUPLICATION";
                    genotype_table_output[i][10] = "DUPLICATION";
                }
            }  
        }
    }
    
    f.close();

    //display a specific record
    if(record){
        std::cout<<"\n\nChr\tStart\tStop\tRegion\tCN_Average\tWeighted_Avg\tCN_SD\tCN_CV\tRegion_Size\tTotal_Windows\tGen_Variation";
        for(int i=0; i<bed_path_rows; i++){
            //std::cout<< "\n" << stats_tabs_table[i][0];
            if(return_record==stats_tabs_table[i][0]){
                std::cout<< "\n" << genotype_table_output[i][1] << "\t"<< genotype_table_output[i][2] << "\t" << genotype_table_output[i][3]
                 << "\t"<< stats_tabs_table[i][0] << "\t"<< genotype_table_output[i][4] << "\t"<< genotype_table_output[i][8] << "\t"<< genotype_table_output[i][5]
                 << "\t"<< genotype_table_output[i][9] << "\t" << genotype_table_output[i][6] << "\t" << genotype_table_output[i][7] << "\t" << genotype_table_output[i][10];
            }
        }
    }

    std::cout<< "\n\nGenotyping Complete!";

    //average abnormalities function
    if(abnormalities){
        count_abnormalities(genotype_table_output , score_check, output_name, bed_path_rows, stats_tabs_table, output_provide, filename);
    }

    //standard deviation abnormalities function
    if(sd){
        count_sd(sd_check, output_name, bed_path_rows, genotype_table_output, stats_tabs_table, output_provide, filename);
    }
    
    //regions function
    if(chr){
        regions(genotype_table_output, stats_tabs_table, bed_path_rows, unique_headers, output_name, chr, start_stop, start, end, return_chromosome, headers);
    }

    delete[] stats_tabs_table;
    delete[] genotype_table_output;

    std::cout << "\n\nProgram Successfully Ran!!\n";
    std::cout << "____________________________________________________________________________________________________"
    "\n____________________________________________________________________________________________________\n";
}
