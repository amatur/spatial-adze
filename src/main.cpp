#include <iostream>
// #include "CLI11.hpp"
#include <string>
#include "allele.hpp"
#include "vcfreader.hpp"
#include "color.hpp"


using namespace std;
void writeOutput(string geofile, GridChar& gridChar, vector<Div>& divs)
{
    std::ifstream inFile(geofile); // Open input file
    std::ofstream outFile("output.tsv"); // Open output file

    if (!inFile) {
        std::cerr << "Error: Unable to open input file." << std::endl;
        return;
    }

    if (!outFile) {
        std::cerr << "Error: Unable to open output file." << std::endl;
        return;
    }

    std::string line;
    std::getline(inFile, line); //  header line
    outFile <<"id\ty\tx\tmean_alpha" << endl; //copy header line to output file

    while (std::getline(inFile, line)) {
        // Read each line from input file
        std::istringstream iss(line);
        std::string col1, col2, col3;
        if (!(iss >> col1 >> col2 >> col3)) {
            std::cerr << "Error: Invalid line format." << std::endl;
            cout<<line<<endl;
            continue;
        }

        double y = std::stod(col2);
        double x = std::stod(col3);
        int gridno = gridChar.getGridNo(x, y);
        if(gridno !=-1){
            // Append "10" as the fourth column
            std::string newLine = col1 + "\t" + col2 + "\t" + col3 + "\t" + std::to_string(divs[gridno].mean_alpha) +"\n";

            // Write modified line to output file
            outFile << newLine;
        }
        
    }
    // Close files
    inFile.close();
    outFile.close();
}

void writeOutputPoplars(string geofile, GridChar& gridChar, vector<Div>& divs)
{
    std::ifstream inFile(geofile); // Open input file
    std::ofstream outFile("output_poplars_20.tsv"); // Open output file

    if (!inFile) {
        std::cerr << "Error: Unable to open input file." << std::endl;
        return;
    }

    if (!outFile) {
        std::cerr << "Error: Unable to open output file." << std::endl;
        return;
    }

    std::string line;
    std::getline(inFile, line); //  header line
    outFile <<"id\ty\tx\tmean_alpha" << endl; //copy header line to output file

    while (std::getline(inFile, line)) {
        // Read each line from input file
        std::istringstream iss(line);
        std::string col1, col2, col3;
        if (!(iss >> col1 >> col2 >> col3)) {
            std::cerr << "Error: Invalid line format." << std::endl;
            cout<<line<<endl;
            continue;
        }

        double x = std::stod(col2);
        double y = std::stod(col3);
        int gridno = gridChar.getGridNo(x, y);
        if(gridno !=-1 && divs[gridno].mean_alpha!=0){
            // Append "10" as the fourth column
            std::string newLine = col1 + "\t" + col2 + "\t" + col3 + "\t" + std::to_string(divs[gridno].mean_alpha) +"\n";

            // Write modified line to output file
            outFile << newLine;
        }
        
    }
    // Close files
    inFile.close();
    outFile.close();


    //make the vectors

    vector<double> v_count(gridChar.numGrids);
    vector<double> v_mean_alpha(gridChar.numGrids);
    for (int i = 0; i < gridChar.numGrids; i++){
        v_count[i] = divs[i].N;
        v_mean_alpha[i] = divs[i].mean_alpha;
        //cout<<v_mean_alpha[i] <<endl;
    }
    // exit(1);
    //draw the dot
    Graphviz gv;
    gv.draw_grid_value(gridChar.a, v_count, v_mean_alpha);
    
}



void adze_main(string vcffile, string geofile, int g){


    Div d;
    int numSnps = 0;
    int numHaps = 0;

    bool haplotype = false;
    int haplen = 1;

    //readAllelesFromVCF(vcffile);
    //readVCF(vcffile, d.Nis);


    //POPLARS
    vcffile= "/Users/amatur/code/spatial-data/poplars/poplars_llhap.txt";
    geofile = "/Users/amatur/code/spatial-data/poplars/poplars_llhap.txt";
    GridChar gridChar(-150, -110, 40, 70,8);
    vector<Div> divs(gridChar.numGrids);    

    if(haplotype){
        read_poplars_haplotype(vcffile, gridChar, divs, numSnps, numHaps);
    }else{
        read_poplars(vcffile, gridChar, divs, numSnps, numHaps);
    }

    cout << "ADZE!!" << endl;
    cout << "vcf: " << vcffile << endl;
    cout << "geo: " << geofile << endl;
    cout << "g: " << g << endl;
    cout << "grid size: " << gridChar.a << endl;



    /*
    //ARABIDOPSIS MAIN
    //GridChar gridChar(-10, 15, 35, 60, 8); //ARABIDOPSIS MAIN
    //GridChar gridChar(-10, 100, -30, 60,9);
    //GridChar gridChar(-10, 55, 35, 100);

    //GridChar gridChar(-10, 15, 35, 60);
    vector<Div> divs(gridChar.numGrids);//0-24
    read012(geofile, vcffile, gridChar, divs, numSnps);
    */


    /* 
    vector<double> sum_log_qs; //size = I
    int numDivs;
    int I = 3;
    vector<vector<int> > qs; // for current g // inner v size = numDivs, outer v size = I
    sum_log_qs.resize(I);
    for(int i = 0; i< I; i++){
        sum_log_qs[i] = 0;
        for(int j = 0; j<numDivs; j++){
            sum_log_qs[i] += log(qs[i][j]);
        }   
    }
    */

    int max_maxg = numHaps;
    // for (int i = 0; i < gridChar.numGrids; i++){
    //      if(divs[i].N >= 2){
    //         if(haplotype){
    //             int s = 0;
    //             for(int x = 0; x < numSnps; x=x+haplen){
    //                 if(divs[i].get_maxg_per_loci(s) < max_maxg){
    //                         max_maxg = divs[i].get_maxg_per_loci(s); 
    //                 }
    //                 s++;
    //             }

    //         }else{
    //             for(int s = 0; s < numSnps; s++){
    //                 if(divs[i].get_maxg_per_loci(s) < max_maxg){
    //                         max_maxg = divs[i].get_maxg_per_loci(s); 
    //                 }
    //             }
    //         }
            
    //      }
    //      cout<<i<<" maxg: "<<max_maxg<<endl;
    // }
    
    max_maxg =2 ;
    //print all element of d.Nis
    for (int i = 0; i < gridChar.numGrids; i++){
        divs[i].gridNo = i;
        if(divs[i].N >= 2){
            double alphasum = 0;

            
            if(haplotype){
                int s = 0;
                for(int x=0; x<numSnps; x=x+haplen){
                    if(i==28 && s<10)
                        cout<<divs[i].compute_alpha(s, max_maxg)<<" ";
                    alphasum+= divs[i].compute_alpha(s, max_maxg);
                    s++;
                }
                
                alphasum/=s*1.0;
            }else{  
                for(int s=0; s<numSnps; s++){
                    if(i==28 && s<10)
                        cout<<divs[i].compute_alpha(s, max_maxg)<<" ";
                    alphasum+= divs[i].compute_alpha(s, max_maxg);
                }
                alphasum/=numSnps*1.0;
            }

            divs[i].mean_alpha = alphasum;
            

            cout<<i<<":"<<divs[i].N << " "<<divs[i].compute_alpha(10)<< " " <<alphasum<<endl;
        }
    }


    writeOutputPoplars(vcffile, gridChar, divs);
    //writeOutput(geofile, gridChar, divs);

    //vector of points
    //vector of ids

}




int main(int argc, char** argv) {
    string vcffile, geofile;
    int g = 5;
    // CLI::App app{"Spatial ADZE"};
    // app.add_option("vcf", vcffile, "VCF")->required();
    // app.add_option("geo", geofile, "GEO")->required();
    // app.add_option("g", g, "Group, g")->required();

    // CLI11_PARSE(app, argc, argv);

    vector<string> args(argv + 1, argv + argc);
    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: adze -v <vcf> -m <metadata> -g <numgroup>" << endl;
            return 0;
        } else if (*i == "-i") {
            vcffile = *++i;
        } else if (*i == "-m") {
            geofile = *++i;
        } else if (*i == "-g") {
            g = std::stoi(*++i);
        } 
    }
    adze_main(vcffile, geofile, g);
    std::cout << "Result: " << g << std::endl;

    return 0;

}
