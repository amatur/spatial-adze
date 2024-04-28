#include<vector>
#include <bits/stdc++.h>
using namespace std;

class GridChar{
    public:
        int a;
        int numGrids ;
        double min_x;
        double max_x;
        double min_y;
        double max_y;
        GridChar(double min_x, double max_x, double min_y, double max_y, int a = 5){
            this->min_x = min_x;
            this->max_x = max_x;
            this->min_y = min_y;
            this->max_y = max_y;
            this->a = a;
            numGrids = a*a; //ax*ay if not square
        }
        int getGridNo(double x, double y){
            if(x<min_x || x>max_x || y<min_y || y>max_y){
                return -1;
            }
            int gx = (x - min_x)/(max_x - min_x)*a;
            int gy = (y - min_y)/(max_y - min_y)*a;
            return gx + a*gy;
        }
        
};
class Point{
    public:
        //GridChar gridChar;
        int gx;
        int gy;

        vector<Point> getNeighbors(){
            vector<Point> neighbors;
            Point p1 = {gx+1, gy};
            Point p2 = {gx-1, gy};
            Point p3 = {gx, gy+1};
            Point p4 = {gx, gy-1};
            neighbors.push_back(p1);
            neighbors.push_back(p2);
            neighbors.push_back(p3);
            neighbors.push_back(p4);
            //also add diagonal neighbors
            Point p5 = {gx+1, gy+1};
            Point p6 = {gx-1, gy-1};
            Point p7 = {gx+1, gy-1};
            Point p8 = {gx-1, gy+1};
            neighbors.push_back(p5);
            neighbors.push_back(p6);
            neighbors.push_back(p7);
            neighbors.push_back(p8);
            return neighbors;
        }
        int gridNo(){
            //return gx + gridChar.a*gy;
        }
};  


class Sample{
    GridChar gridChar;
    int x; // actual x = from long lat
    int y;
    int div_gridNo; // assign  gridno to each sample
    void setGridNo(){
        this->div_gridNo = x + gridChar.a*y;
    }
};

class Div{
    public:
    int gridNo;
    //int gridX;
    //int gridY;
    bool turnoff = false;
    double mean_alpha = 0;
    int N = 0; // number of individuals in this div, i.e. grid
    vector<vector<int> > Nis_per_loci; // let N=18, type 0: 10 type 1: 5, type 2: 3, 3: missing sum should be N=18; //allele counts
    // void getNeighbors(){
    //     //(x+1, y); // if not empty, add to neighbor list
    // }
    
    unsigned nChoosek( unsigned n, unsigned k )
    {
        if (k > n) return 0;
        if (k * 2 > n) k = n-k;
        if (k == 0) return 1;

        int result = n;
        for( int i = 2; i <= k; ++i ) {
            result *= (n-i+1);
            result /= i;
        }
        return result;
    }

    //assume qs, Nis are already computed
    double compute_pi(int g, int loci, const vector<vector<int> >& qs, vector<double>& sum_log_qs){
        //sum over all p_ig's
        int I = Nis_per_loci[loci].size() - 1;

        double sum = 0;
        for(int i = 0; i<I; i++){
            //compute p_ig
            sum += (1 - qs[i][gridNo])*(exp(sum_log_qs[i]-log(qs[i][gridNo])));
        }
        return sum;
    }

    //allelic richness, populate qs
    double compute_alpha(int g, vector<vector<int> >& qs, int loci){
        //sum over all p_ig's
        int I = Nis_per_loci[loci].size() - 1;
        vector<int> Nis = Nis_per_loci[loci];
        double sum = 0;
        for(int i = 0; i<I; i++){
            //compute p_ig
            qs[i][gridNo] = (nChoosek(N-Nis[i],g)*1.0/nChoosek(N,g));
            sum += 1 - qs[i][gridNo];
        }
        return sum;
    }
    //allelic richness
    double compute_alpha(int loci, int g=5){
        //sum over all p_ig's
        int I = Nis_per_loci[loci].size() - 1;
        vector<int> Nis = Nis_per_loci[loci];
        double sum = 0;
        for(int i = 0; i<I; i++){
            //compute p_ig
            sum += 1 - (nChoosek(N-Nis[i],g)*1.0/nChoosek(N,g));
        }
        return sum;
    }
};

class ADZE_params{
    //int I;
    vector<double> sum_log_qs; //size = I
    //int g;
    //int nloci;
    //int nhaps;
    int numDivs;
    //vector<Div> divs;
    vector<vector<int> > qs; // for current g // inner v size = numDivs, outer v size = I

    //do this after alpha is computed, each time g is updated
    void init_sum_log_qs(int I){
        sum_log_qs.resize(I);
        for(int i = 0; i< I; i++){
            sum_log_qs[i] = 0;
            for(int j = 0; j<numDivs; j++){
                sum_log_qs[i] += log(qs[i][j]);
            }   
        }
    }
};



/*
Output format:


x y gridNo long lat pi alpha

*/