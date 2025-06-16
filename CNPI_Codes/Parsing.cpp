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


//gzcat File.gz | grep -v "#" | awk '$1 ~ /^(1?[0-9]|2[0-2]|X|Y)$/' | sort -k1,1V -k2,2n | gzip > File.sorted_filtered.gz

static struct option long_options[] = {
    {"cn_path", required_argument, NULL, 'c'},
    {"g_path", required_argument, NULL, 'g'},
    {"vcf_path", required_argument, NULL, 'v'},
    {"CNVnator_path", required_argument, NULL, 'n'},
    {"help", no_argument, NULL, 'h'},
    {0, 0, 0, 0}
};

int main(int argc, char *argv[]){



    std::cout << ("\n\nTime to Parse through a File!!\n\n");


    std::string region_file_path;    // for reference bed file. Function can open .bed and .bed.gz
    int region_file_rows =0;
    bool cn_parse=false, vcf_parse=false, gfile_parse=false, include_chrm=false, CNVnator_parse=false;


    //arguments short names and argument manager
   int opt;
    while ((opt = getopt_long_only(argc, argv, "c:g:v:n:mh", long_options, NULL)) != -1) {
        switch (opt) {
            case 'c':
                region_file_path = optarg;
                cn_parse = true;
                break;
            case 'g':
                region_file_path = optarg;
                gfile_parse =true;
                break;
            case 'v':
                region_file_path = optarg;
                vcf_parse = true;
                break;
            case 'm':
                include_chrm=true;
                break;
            case 'n':
                region_file_path = optarg;
                CNVnator_parse = true;
                break;
            case 'h':
                std::cout<< "\nTime to parse a file"
                << "\n-g to parse a gff3 or gtf file"
                << "\n-v to parse a VCF file"
                << "\n-m Include chromosome labels if they are not within VCF file"
                << "\n-c to parse a Bed File or a QuicK-mer2 file"
                << "\n-n to parse CNVnator output";
                return 0; 
        }
    }


    std::cout << "\n____________________________________________________________________________________________________"
    "\n____________________________________________________________________________________________________";
    std::cout << "\n\nThe files are " << region_file_path;
    
    
    //next number of files are for finding the amount of rows and reading the gz counts file
    gzFile cn_table = gzopen(region_file_path.c_str(), "rb");
    if (!cn_table) {
        std::cerr << "\nError opening input file!\nCheck if you are entering a valid file or a valid file format!";
        return 0;
    }

    //unzipping and finding lines of .gz file
    char buffer[16384]; // Adjust buffer size as needed
    int bytes_read;
    std::string lines;


    //opening the reference bed file
    gzFile gz_bed_table = gzopen(region_file_path.c_str(), "rb");
    if (!cn_table) {
        std::cerr << "\nError opening input file! Check if you are entering a valid file or valid file format!";
        return 0;
    }

    //unzipping or opening refseq file and finding the size
    char buff[1024]; // Adjust buffer size as needed
    int byte_read;

    
    //finding the amount of rows in the bed ncbi refseq file
    //only ran if there has not been a provided amount of rows
    while ((byte_read = gzread(gz_bed_table, buff, sizeof(buff))) > 0) {
        for (int i = 0; i < byte_read; ++i) {
            if (buff[i] != '\n') {
                lines += buff[i];
            }else{
                lines.clear(); // Clear the line for the next iteration
                region_file_rows +=1; //rows in reference file
            }
        }
    }
    gzrewind(gz_bed_table); //set back to zero
    std::cout << "\n\nThe File rows are " << region_file_rows;
    
    
    //setting up bed table
   	std::string **bed_tabs_table = new std::string*[region_file_rows];

	for (int i=0; i<region_file_rows; i++){
		bed_tabs_table[i]=new std::string[4];
	}

    for(int i=0; i<region_file_rows; i++){
        for(int j=0; j<4; j++){
            bed_tabs_table[i][j]=""; 
        }
    }
    //end setting up bed ncbi table rows and initilizations
    if(cn_parse==true){
        //populating the refseq table
        int tabs =0, cn_file_rows=0; //tabs and rows when reading file
        int unique_headers =1;  //for finding amount of headers. At least 1 unique header
        std::vector<std::string> chrms;
        std::vector<int> indexes;

        
        while ((byte_read = gzread(gz_bed_table, buff, sizeof(buff))) > 0) {
            for (int i = 0; i < byte_read; ++i) {
                if(buff[i] == '\t') {
                    tabs +=1;       //tab increased so put information at next tab   
                }else if(buff[i]=='\n'){
                    if(cn_file_rows>1 &&(bed_tabs_table[cn_file_rows][0] != bed_tabs_table[cn_file_rows-1][0])){
                            unique_headers +=1;     //found additional unique header from reference bed file
                            chrms.push_back(bed_tabs_table[cn_file_rows][0]);
                            indexes.push_back(cn_file_rows);
                    }else if (cn_file_rows==0){
                        chrms.push_back(bed_tabs_table[cn_file_rows][0]);
                        indexes.push_back(cn_file_rows);
                    }
                    cn_file_rows+=1; //iterate to next row
                    tabs=0;         //set tabs back to zero

                }else if(tabs<=3){
                   if(tabs==0){
                    bed_tabs_table[cn_file_rows][0] += buff[i];
                   }if(tabs==1){
                    bed_tabs_table[cn_file_rows][1] += buff[i];
                   }if(tabs==2){
                    bed_tabs_table[cn_file_rows][2] += buff[i];
                   }else{
                    bed_tabs_table[cn_file_rows][3] +=buff[i];
                   }
                }
            }
        }
        indexes.push_back(region_file_rows);
        std::cout << "\nThe Unique Headers are: " << unique_headers;

        std::string CN_file;
        CN_file = "Location_reference_file.txt";

        std::ofstream c(CN_file);

        for(int i=0; i<region_file_rows; i++){
            //std::cout<< "\n";
            c << bed_tabs_table[i][0] << "\t" << bed_tabs_table[i][1] << "\t" << bed_tabs_table[i][2] << "\tRegion:" << bed_tabs_table[i][1] << "-" << bed_tabs_table[i][2] << "\n";

        }
        c.close();

           // Printing the vector
        std::string index_file;
        index_file = "CN_index_file.txt";
        std::ofstream i(index_file);
        i << unique_headers << "\n";
        std::cout<< "\n";
        for (int num: indexes) {
            std::cout << num << " ";
            i << num << " ";
        }i << "\n";
        std:: cout<< "\n";
        for (const std::string& chrms : chrms) {
            std::cout << chrms << " ";
            i << chrms << " ";
        }std::cout<< "\n";
        i.close();

    }else if(vcf_parse==true){
        int tabs =0, vcf_file_rows=0; //tabs and rows when reading file
        int unique_headers =1;  //for finding amount of headers. At least 1 unique header
        bool prefix = true, skip=false, stop=false, type_skip=false, type_stop=false;
        std::string val ="";
        std::string svtype="";
        std::vector<std::string> chrms;
        std::vector<int> indexes;

        
        while ((byte_read = gzread(gz_bed_table, buff, sizeof(buff))) > 0) {
            for (int i = 0; i < byte_read; ++i) {
                if(buff[i] == '\t') {
                    tabs +=1;       //tab increased so put information at next tab   
                }else if(buff[i]=='\n'){
                    if(vcf_file_rows>1 &&(bed_tabs_table[vcf_file_rows][0] != bed_tabs_table[vcf_file_rows-1][0])){
                            unique_headers +=1;     //found additional unique header from reference bed file
                            chrms.push_back(bed_tabs_table[vcf_file_rows][0]);
                            indexes.push_back(vcf_file_rows);
                    }else if (vcf_file_rows==0){
                        chrms.push_back(bed_tabs_table[vcf_file_rows][0]);
                        indexes.push_back(vcf_file_rows);
                    }
                    if(bed_tabs_table[vcf_file_rows][2].empty()){
                        bed_tabs_table[vcf_file_rows][2]="NA";
                    }
                    bed_tabs_table[vcf_file_rows][3] = ("Region:" + bed_tabs_table[vcf_file_rows][1] + "-" + bed_tabs_table[vcf_file_rows][2]) + "_" + svtype;
                    vcf_file_rows+=1; //iterate to next row
                    tabs=0;         //set tabs back to zero
                    prefix = true;
                    val="";
                    skip=false;
                    stop=false;
                    svtype="";
                    type_skip=false;
                    type_stop=false;

                }else{
                    if(tabs==1){
                        bed_tabs_table[vcf_file_rows][1] +=buff[i];
                    }else if(tabs==7){
                        val += buff[i];
                        if (val.length() >= 4 && skip==false) {
                            std::string last4 = val.substr(val.length() - 4);
                            if(last4=="END="){
                                if(val.length()>=6 && skip==false){
                                    std::string last6 = val.substr(val.length() - 6);
                                    if(last6!="CIEND="){
                                        skip=true;}
                                }
                            }
                        }else if(skip==true && stop==false){
                            if (buff[i]!=';'){bed_tabs_table[vcf_file_rows][2] +=buff[i];}
                            else if(buff[i]==';'){stop=true;}
                        }

                        if(val.length()>=7 && type_skip ==false){
                            std::string last7 = val.substr(val.length() - 7);
                            if(last7=="SVTYPE="){
                                type_skip=true;
                            }
                        }else if(type_stop==false && type_skip==true){
                            if (buff[i]!=';'){svtype +=buff[i];}
                            else if(buff[i]==';'){type_stop=true;}
                        }
                    }else if(tabs==0){
                        if(prefix==true){
                            if(include_chrm==true){
                                bed_tabs_table[vcf_file_rows][0] = "chr";}
                            prefix=false;}
                        bed_tabs_table[vcf_file_rows][0] +=buff[i];
                    }
                }
            }
        }
        indexes.push_back(region_file_rows);
        std::cout << "\nThe Unique Headers are: " << unique_headers;

        std::string file;
        std::string NA_file;
        file = "VCF_parsed_file.txt";
        NA_file = "VCF_NA_paresed_file.txt";

        std::ofstream f(file);
        std::ofstream n(NA_file);

        for(int i=0; i<region_file_rows; i++){
            //std::cout<< "\n";
            std::string last3 = bed_tabs_table[i][3].substr(bed_tabs_table[i][3].length() - 3);
    
            if(bed_tabs_table[i][2]=="NA" || last3=="INV" || last3=="INS" || last3=="BND"){
                for(int j=0; j<4; j++){
                    n<<bed_tabs_table[i][j] <<"\t";
                }n<<"\n";
            }else{
                for(int j=0; j<4; j++){
                    f<<bed_tabs_table[i][j] <<"\t";
                }f<<"\n";
            }
        }

        f.close();
            // Printing the vector
        std::string index_file;
        index_file = "VCF_index_file.txt";
        std::ofstream i(index_file);
        i << unique_headers << "\n";
        std::cout<< "\n";
        for (int num: indexes) {
            std::cout << num << " ";
            i << num << " ";
        }i << "\n";
        std:: cout<< "\n";
        for (const std::string& chrms : chrms) {
            std::cout << chrms << " ";
            i << chrms << " ";
        }std::cout<< "\n";
        i.close();

    }else if(CNVnator_parse==true){
        int spaces =0, gz_gff3_rows=0; //tabs and rows when reading file
        int unique_headers =1;  //for finding amount of headers. At least 1 unique header
        bool prefix = true, colon=false, hyphen=false;
        std::vector<std::string> chrms;
        std::vector<int> indexes;
    
        while ((byte_read = gzread(gz_bed_table, buff, sizeof(buff))) > 0) {
            for (int i = 0; i < byte_read; ++i) {
                if(buff[i] == ' ') {
                    spaces +=1;       //tab increased so put information at next tab
                }else if(buff[i]=='\n'){
                    if(gz_gff3_rows>=1 &&(bed_tabs_table[gz_gff3_rows][0] != bed_tabs_table[gz_gff3_rows-1][0])){
                            unique_headers +=1;     //found additional unique header from reference bed file
                            chrms.push_back(bed_tabs_table[gz_gff3_rows][0]);
                            indexes.push_back(gz_gff3_rows);
                    }else if (gz_gff3_rows==0){
                        chrms.push_back(bed_tabs_table[gz_gff3_rows][0]);
                        indexes.push_back(gz_gff3_rows);
                    }
                    gz_gff3_rows+=1; //iterate to next row
                    spaces=0;         //set tabs back to zero
                    prefix = true;
                    colon = false;
                    hyphen = false;
                }else if(spaces==0 || spaces==2 || spaces==4){
                    continue;   
                }else if(spaces==1 || spaces==3){
                    if(spaces==1 && buff[i]==':'){
                        colon=true;
                        continue;
                    }else if(spaces==1 && colon==true && buff[i]=='-'){
                            hyphen=true;
                            continue;
                    }else if(spaces==1 && colon==false){
                        bed_tabs_table[gz_gff3_rows][0] +=buff[i];
                    }else if(spaces==1 && colon==true && hyphen==false){
                        bed_tabs_table[gz_gff3_rows][1] +=buff[i];
                    }else if(spaces==1 && colon==true && hyphen==true){
                        bed_tabs_table[gz_gff3_rows][2] +=buff[i];
                    
                    }else if (spaces==3){
                        bed_tabs_table[gz_gff3_rows][3]+=buff[i];
                    }
                    //bed_tabs_table[gz_gff3_rows][tabs] +=buff[i];
                }
            }
        }

        for (int i=0; i<region_file_rows; i++){
            std::cout<< "\n" << i << "\t" << bed_tabs_table[i][0] << "\t" << bed_tabs_table[i][1] << "\t" << bed_tabs_table[i][2] << "\t" << bed_tabs_table[i][3];
        }
        indexes.push_back(region_file_rows);
        std::cout << "\nThe Unique Headers are: " << unique_headers;

        std::string file;
        file ="CNVnator_parsed_file.txt";
        std::ofstream f(file);

        for(int i=0; i<region_file_rows; i++){
            for(int j=0; j<4; j++){
                if(j==0){
                    f<<"chr";
                }
                f <<bed_tabs_table[i][j] <<"\t";
            }f << "\n";
        }
        f.close();

            // Printing the vector
        std::string index_file;
        index_file = "CNVnator_index_file.txt";
        std::ofstream i(index_file);
        i << unique_headers << "\n";
        std::cout<< "\n";
        for (int num: indexes) {
            std::cout << num << " ";
            i << num << " ";
        }i << "\n";
        std:: cout<< "\n";
        for (const std::string& chrms : chrms) {
            std::cout << "chr" << chrms << " ";
            i << "chr" << chrms << " ";
        }std::cout<< "\n";
        i.close();
    }if(cn_parse==true){
        //populating the refseq table
        int tabs =0, cn_file_rows=0; //tabs and rows when reading file
        int unique_headers =1;  //for finding amount of headers. At least 1 unique header
        std::vector<std::string> chrms;
        std::vector<int> indexes;

        
        while ((byte_read = gzread(gz_bed_table, buff, sizeof(buff))) > 0) {
            for (int i = 0; i < byte_read; ++i) {
                if(buff[i] == '\t') {
                    tabs +=1;       //tab increased so put information at next tab   
                }else if(buff[i]=='\n'){
                    if(cn_file_rows>1 &&(bed_tabs_table[cn_file_rows][0] != bed_tabs_table[cn_file_rows-1][0])){
                            unique_headers +=1;     //found additional unique header from reference bed file
                            chrms.push_back(bed_tabs_table[cn_file_rows][0]);
                            indexes.push_back(cn_file_rows);
                    }else if (cn_file_rows==0){
                        chrms.push_back(bed_tabs_table[cn_file_rows][0]);
                        indexes.push_back(cn_file_rows);
                    }
                    cn_file_rows+=1; //iterate to next row
                    tabs=0;         //set tabs back to zero

                }else if(tabs<=3){
                   if(tabs==0){
                    bed_tabs_table[cn_file_rows][0] += buff[i];
                   }if(tabs==1){
                    bed_tabs_table[cn_file_rows][1] += buff[i];
                   }if(tabs==2){
                    bed_tabs_table[cn_file_rows][2] += buff[i];
                   }else{
                    bed_tabs_table[cn_file_rows][3] +=buff[i];
                   }
                }
            }
        }
        
        indexes.push_back(region_file_rows);
        std::cout << "\nThe Unique Headers are: " << unique_headers;

        std::string CN_file;
        CN_file = "Location_reference_file.txt";

        std::ofstream c(CN_file);

        for(int i=0; i<region_file_rows; i++){
            //std::cout<< "\n";
            c << bed_tabs_table[i][0] << "\t" << bed_tabs_table[i][1] << "\t" << bed_tabs_table[i][2] << "\tRegion:" << bed_tabs_table[i][1] << "-" << bed_tabs_table[i][2] << "\n";

        }
        c.close();

           // Printing the vector
        std::string index_file;
        index_file = "CN_index_file.txt";
        std::ofstream i(index_file);
        i << unique_headers << "\n";
        std::cout<< "\n";
        for (int num: indexes) {
            std::cout << num << " ";
            i << num << " ";
        }i << "\n";
        std:: cout<< "\n";
        for (const std::string& chrms : chrms) {
            std::cout << chrms << " ";
            i << chrms << " ";
        }std::cout<< "\n";
        i.close();

    }
}