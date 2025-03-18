#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "gzstream.h"
using namespace std;


void readGeo(string filename, vector<int>& ids){
    // cerr << "Opening " << filename << "...\n";
    // ifstream file(filename);
    // if (!file.is_open()) {
    //     cerr << "Error: Unable to open file " << filename << endl;
    //     return 1;
    // }

    // vector<int> ids;
    // vector<double> x_values;
    // vector<double> y_values;

    // string line;
    // getline(file, line); // Skip header
    // while (getline(file, line)) {
    //     stringstream ss(line);
    //     string id;
    //     double x, y;
    //     if (ss >> id >> x >> y) {
    //         ids.push_back(id);
    //         x_values.push_back(x);
    //         y_values.push_back(y);
    //     }
    // }

    // file.close();
}

int countFields(string& str)
{
    string::const_iterator it;
    int result;
    int numFields = 0;
    int seenChar = 0;
    for (it = str.begin() ; it < str.end(); it++)
    {
        result = isspace(*it);
        if (result == 0 && seenChar == 0)
        {
            numFields++;
            seenChar = 1;
        }
        else if (result != 0)
        {
            seenChar = 0;
        }
    }
    return numFields;
}



// poplars_llhap.txt
/***
 * biallelic read
 */
void read_poplars(string variantfilename, GridChar &gridChar, vector<Div> &divs, int &numSnps, int &numHaps)
{
    // for(map<string,int>::iterator it = Nis_map.begin(); it != Nis_map.end(); it++){
    //     Nis.push_back(it->second);
    // }

    ifstream vcfFile;
    vcfFile.open(variantfilename.c_str());
    if (vcfFile.fail())
    {
        cerr << "Error: Unable to open file " << variantfilename << endl;
        return;
    }

    string line;
    getline(vcfFile, line);
    // cout<<line<<endl;

    numSnps = countFields(line) - 3;
    numHaps = 0;
    // Nis.resize(numSnps);
    //  for(int i = 0; i<numSnps; i++){
    //     Nis[i].resize(4);
    //     Nis[i][0] = 0;
    //     Nis[i][1] = 0;
    //     Nis[i][2] = 0;
    //     Nis[i][3] = 0;

    //  }

    do
    {
        stringstream ss(line);

        string allele;
        int id; //drainage
        double x = -100;
        double y = -100;
        ss >> id >> x >> y;
        // cout<< id << " " << x << " " << y << endl;
        if (x >= gridChar.min_x && x <= gridChar.max_x && y >= gridChar.min_y && y <= gridChar.max_y)
        //if(true)
        {

            int gridno = gridChar.getGridNo(x, y);
            if (false)
                cout << id << " " << x << " " << y << " " << " " << gridno << " ";

            divs[gridno].gridNo = gridno;
            divs[gridno].N++;

            numHaps++;

            if (divs[gridno].Nis_per_loci.size() != numSnps)
            {
                divs[gridno].Nis_per_loci.resize(numSnps);
            }

            //ss >> allele; // skip : read the entry

            //if (false)
            //    cout << allele << endl; // the sample id that is in grid

            for (int i = 0; i < numSnps; i++)
            {
                if (divs[gridno].Nis_per_loci[i].size() != 4)
                { // either 4 or 0
                    divs[gridno].Nis_per_loci[i].resize(4);
                    divs[gridno].Nis_per_loci[i][0] = 0;
                    divs[gridno].Nis_per_loci[i][1] = 0;
                    divs[gridno].Nis_per_loci[i][2] = 0;
                    divs[gridno].Nis_per_loci[i][3] = 0;
                }

                if (ss >> allele)
                {
                    if (allele == "0")
                    {
                        divs[gridno].Nis_per_loci[i][0] += 1;
                    }
                    else if (allele == "1")
                    {
                        divs[gridno].Nis_per_loci[i][1] += 1;
                    }
                    else if (allele == "2")
                    {
                        divs[gridno].Nis_per_loci[i][2] += 1;
                    }
                    else
                    {                                         //-1
                        divs[gridno].Nis_per_loci[i][3] += 1; // missing
                        // cerr << "Error: missing allele " << allele << endl;
                        // return;
                    }
                }
            }
        }
        else
        {
            // for(int i = 0; i<numSnps; i++){
            //     ss >> allele; //skip
            // }
        }
    } while (getline(vcfFile, line));

    vcfFile.close();
    // do the Nis again; filter max missing

    for (int i = 0; i < gridChar.numGrids; i++)
    {
        if (divs[i].N < 5)
        {
            divs[i].turnoff = true;
        }
        int sampleloci = 100;
        if (divs[i].N != 0)
        {
            cout << "grid " << i << ":" << divs[i].N << " " << divs[i].Nis_per_loci[10][0] << " " << divs[i].Nis_per_loci[sampleloci][1] << " " << divs[i].Nis_per_loci[sampleloci][0] << " " << divs[i].Nis_per_loci[sampleloci][3] << endl;
        }
    }
}


// poplars_llhap.txt
/***
 * biallelic read
 */
void read_poplars_haplotype(string variantfilename, GridChar &gridChar, vector<Div> &divs, int &numSnps, int &numHaps)
{
     
    ifstream vcfFile;
    vcfFile.open(variantfilename.c_str());
    if (vcfFile.fail())
    {
        cerr << "Error: Unable to open file " << variantfilename << endl;
        return;
    }

    string line;
    getline(vcfFile, line);
    // cout<<line<<endl;

    numSnps = countFields(line) - 3;
    numHaps = 0;
    // Nis.resize(numSnps);
    //  for(int i = 0; i<numSnps; i++){
    //     Nis[i].resize(4);
    //     Nis[i][0] = 0;
    //     Nis[i][1] = 0;
    //     Nis[i][2] = 0;
    //     Nis[i][3] = 0;

    //  }

    do
    {
        stringstream ss(line);

        string allele;
        int id; //drainage
        double x = -100;
        double y = -100;
        ss >> id >> x >> y;
        // cout<< id << " " << x << " " << y << endl;
        if (x >= gridChar.min_x && x <= gridChar.max_x && y >= gridChar.min_y && y <= gridChar.max_y)
        //if(true)
        {

            int gridno = gridChar.getGridNo(x, y);
            if (false)
                cout << id << " " << x << " " << y << " " << " " << gridno << " ";

            divs[gridno].gridNo = gridno;
            divs[gridno].N++;
            vector<map<string, int> >& haplotype_map = divs[gridno].haplotype_map;

            numHaps++;

            // if (divs[gridno].Nis_per_loci.size() != numSnps)
            // {
            //     divs[gridno].Nis_per_loci.resize(numSnps);
            // }

            //ss >> allele; // skip : read the entry

            //if (false)
            //    cout << allele << endl; // the sample id that is in grid

            int haplen = 1;
            int vectoridx = 0;
            for (int i = 0; i < numSnps; i=i+haplen)
            {
                if(haplotype_map.size() < vectoridx+1)
                    haplotype_map.push_back(map<string, int>());


                allele = "";
                for (int a = 0; a < haplen; a++){
                    string achar;
                    ss >> achar;
                    allele+=achar;
                }

                // allele.resize( haplen );
                // ss.read( &allele[0], haplen );
                //cout<<"sl "<<allele<<endl;
                if (haplotype_map[vectoridx].find(allele) == haplotype_map[vectoridx].end())
                {
                    haplotype_map[vectoridx][allele] = 1;
                }
                else
                {
                    haplotype_map[vectoridx][allele]++;
                }

                vectoridx++;
            }
        }
        else
        {
            // for(int i = 0; i<numSnps; i++){
            //     ss >> allele; //skip
            // }
        }
    } while (getline(vcfFile, line));

    vcfFile.close();
    // do the Nis again; filter max missing

    for (int i = 0; i < gridChar.numGrids; i++)
    {   
        if (divs[i].N < 8)
            continue;
        //cout<< "grid "<< i << " "<<divs[i].haplotype_map.size() <<endl;;

        divs[i].Nis_per_loci.resize(divs[i].haplotype_map.size());
        //iterate through the map haplotype_map
        for (int j = 0; j < divs[i].haplotype_map.size(); j++)
        {
                
            for (map<string, int>::iterator it = divs[i].haplotype_map[j].begin(); it != divs[i].haplotype_map[j].end(); it++)
            {
                //cout << "grid "<< i << " "<< j << ":" << it->first << " " << it->second << endl;
                divs[i].Nis_per_loci[j].push_back(it->second);

            }
            // cout<<divs[i].Nis_per_loci[j].size()<<endl;
        }
       
        
        // // if (divs[gridno].Nis_per_loci.size() != numSnps)
        // //     {
        // //         divs[gridno].Nis_per_loci.resize(numSnps);
        // //     }

        // if (divs[i].N < 5)
        // {
        //     divs[i].turnoff = true;
        // }
        // int sampleloci = 100;
        // if (divs[i].N != 0)
        // {
        //     cout << "grid " << i << ":" << divs[i].N << " " << divs[i].Nis_per_loci[10][0] << " " << divs[i].Nis_per_loci[sampleloci][1] << " " << divs[i].Nis_per_loci[sampleloci][0] << " " << divs[i].Nis_per_loci[sampleloci][3] << endl;
        // }
    }

    for (int i = 0; i < gridChar.numGrids; i++)
    {   
        divs[i].haplotype_map.clear();
    }
}

/***
 * biallelic read
 */
void read012(string geofilename, string variantfilename, GridChar &gridChar, vector<Div> &divs, int &numSnps)
{
    // for(map<string,int>::iterator it = Nis_map.begin(); it != Nis_map.end(); it++){
    //     Nis.push_back(it->second);
    // }
    ifstream geoFile;
    geoFile.open(geofilename.c_str());
    if (geoFile.fail())
    {
        cerr << "Error: Unable to open file " << geofilename << endl;
        return;
    }

    ifstream vcfFile;
    vcfFile.open(variantfilename.c_str());
    if (vcfFile.fail())
    {
        cerr << "Error: Unable to open file " << variantfilename << endl;
        return;
    }

    string geoline;
    string line;
    getline(geoFile, geoline); // skip header

    getline(geoFile, geoline);

    getline(vcfFile, line);
    // cout<<line<<endl;

    numSnps = countFields(line) - 1;
    int numHaps = 0;
    // Nis.resize(numSnps);
    //  for(int i = 0; i<numSnps; i++){
    //     Nis[i].resize(4);
    //     Nis[i][0] = 0;
    //     Nis[i][1] = 0;
    //     Nis[i][2] = 0;
    //     Nis[i][3] = 0;

    //  }

    do
    {
        stringstream ss_geo(geoline);

        stringstream ss(line);

        string allele;
        int id;
        double x = -100;
        double y = -100;
        ss_geo >> id >> y >> x;
        // cout<< id << " " << x << " " << y << endl;
        if (x >= gridChar.min_x && x <= gridChar.max_x && y >= gridChar.min_y && y <= gridChar.max_y)
        {

            int gridno = gridChar.getGridNo(x, y);
            if (false)
                cout << id << " " << x << " " << y << " " << " " << gridno << " ";

            divs[gridno].gridNo = gridno;
            divs[gridno].N++;

            numHaps++;

            if (divs[gridno].Nis_per_loci.size() != numSnps)
            {
                divs[gridno].Nis_per_loci.resize(numSnps);
            }

            ss >> allele; // skip : read the entry

            if (false)
                cout << allele << endl; // the sample id that is in grid

            for (int i = 0; i < numSnps; i++)
            {
                if (divs[gridno].Nis_per_loci[i].size() != 4)
                { // either 4 or 0
                    divs[gridno].Nis_per_loci[i].resize(4);
                    divs[gridno].Nis_per_loci[i][0] = 0;
                    divs[gridno].Nis_per_loci[i][1] = 0;
                    divs[gridno].Nis_per_loci[i][2] = 0;
                    divs[gridno].Nis_per_loci[i][3] = 0;
                }

                if (ss >> allele)
                {
                    if (allele == "0")
                    {
                        divs[gridno].Nis_per_loci[i][0] += 1;
                    }
                    else if (allele == "1")
                    {
                        divs[gridno].Nis_per_loci[i][1] += 1;
                    }
                    else if (allele == "2")
                    {
                        divs[gridno].Nis_per_loci[i][2] += 1;
                    }
                    else
                    {                                         //-1
                        divs[gridno].Nis_per_loci[i][3] += 1; // missing
                        // cerr << "Error: missing allele " << allele << endl;
                        // return;
                    }
                }
            }
        }
        else
        {
            // for(int i = 0; i<numSnps; i++){
            //     ss >> allele; //skip
            // }
        }
    } while (getline(vcfFile, line) && getline(geoFile, geoline));

    vcfFile.close();
    geoFile.close();
    // do the Nis again; filter max missing

    for (int i = 0; i < gridChar.numGrids; i++)
    {
        if (divs[i].N < 5)
        {
            divs[i].turnoff = true;
        }
        int sampleloci = 100;
        if (divs[i].N != 0)
        {
            cout << "grid " << i << ":" << divs[i].N << " " << divs[i].Nis_per_loci[10][0] << " " << divs[i].Nis_per_loci[sampleloci][1] << " " << divs[i].Nis_per_loci[sampleloci][0] << " " << divs[i].Nis_per_loci[sampleloci][3] << endl;
        }
    }
}

// Function to read alleles from a VCF file
void readAllelesFromVCF(const string& filename) {
    //vcftools --gzvcf chr1_1k.vcf.gz --max-missing 0 --remove-indels --remove-filtered-all --recode --out filtered_output

    igzstream vcfFile;
    vcfFile.open(filename.c_str());
    if (vcfFile.fail()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return;
    }
    // if (!vcfFile.is_open()) {
    //     cerr << "Error: Unable to open file " << filename << endl;
    //     return;
    // }

    string line;
    while (getline(vcfFile, line)) {
        cout<<line<<endl;
        // Skip header lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Split the line into fields
        vector<string> fields;
        size_t pos = 0;
        while ((pos = line.find('\t')) != string::npos) {
            fields.push_back(line.substr(0, pos));
            line.erase(0, pos + 1);
        }
        fields.push_back(line); // Last field

        // Assuming alleles are in the 4th and 5th columns (1-based indexing)
        string allele1 = fields[3];
        string allele2 = fields[4];

        // Print alleles
        cout << "Alleles: " << allele1 << ", " << allele2 <<":" <<  fields[5] << " " << fields[6] << endl;
    }

    vcfFile.close();
}

void readVCF(string filename,vector<int>& Nis)
{
    map<string, int> Nis_map;
    bool unphased = true;
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int numMapCols = 9;
    //int fileStart = fin.tellg();
    string line;
    int nloci = 0;
    int previous_nhaps = -1;
    int current_nhaps = 0;

    int skipcount = 0;
    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        
        if (line[0] == '#') {
            continue;
        }
        cout<<line<<endl;
        //getline(fin,line);
        //if(fin.eof()) break;
        nloci++;
        current_nhaps = countFields(line);
        //cout << "nloci: " << current_nhaps << endl;
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps)
        {
            cout << "ERROR: line " << nloci << " of " << filename << " has " << current_nhaps
                 << " fields, but the previous line has " << previous_nhaps << " fields.\n";
                 cout<<line<<endl;
            //throw 0;
        }
        previous_nhaps = current_nhaps;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int nhaps = unphased ? (current_nhaps - numMapCols) : (current_nhaps - numMapCols) * 2;
    int nfields = (current_nhaps - numMapCols);

    cerr << "Loading " << nhaps << " haplotypes and " << nloci << " loci...\n";
    int nloci_before_filtering = nloci;
    nloci = 0;


    string junk;
    char allele1, allele2, separator;
    bool skipLine = false;
    for (int locus = 0; locus < nloci_before_filtering; locus++)
    {
        for (int i = 0; i < numMapCols; i++) {
            fin >> junk;
            if (i == 0 && junk[0] == '#') {
                skipLine = true;
                break;
            }
        }
        if (skipLine) {
            getline(fin, junk);
            skipLine = false;
            locus--;
            continue;
        }
        for (int field = 0; field < nfields; field++)
        {
            fin >> junk;
            allele1 = junk[0];
            separator = junk[1];
            allele2 = junk[2];
            if(allele1!=allele2){
                cout<<allele1<<separator<<allele2<<endl;
            }
            
            if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
            {
                //cout<<"SKIPPING: Alleles must be coded 0/1 only. Found "<<allele1 << " " << allele2 << "\n"<<endl;
                // cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                // cerr << allele1 << " " << allele2 << endl;
                // throw 0;
            }

            //if(separator != '|'){
            //    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
            //    throw 0;
            //}
            if(unphased){
                char allele = '0';
                if (allele1 == '1' && allele2 == '1'){
                    allele = '2';
                }
                else if (allele1 == '1' || allele2 == '1'){
                    allele = '1';

                }
                //else allele = '0' implied;
                
                if(Nis_map.find(""+allele) == Nis_map.end()){
                    Nis_map[""+allele] = 1;
                }else{
                    Nis_map[""+allele]++;
                }
            }  
        }
    }
    //copy content of Nis_map to Nis
    for(map<string,int>::iterator it = Nis_map.begin(); it != Nis_map.end(); it++){
        Nis.push_back(it->second);
    }
 

    fin.close();
}
