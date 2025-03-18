
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>


using namespace std;

class Color {
private:
    int r;
    int g;
    int b;

public:
    Color(int rr=0, int gg=0, int bb=0) : r(rr), g(gg), b(bb) {}

    std::string rgbToHex() {
        std::stringstream ss;
        ss << "#" << std::hex << std::setw(2) << std::setfill('0') << r
           << std::hex << std::setw(2) << std::setfill('0') << g
           << std::hex << std::setw(2) << std::setfill('0') << b;
        return ss.str();
    }

    // alpha in range 0,1
    Color static get_color_from_value(double alpha, const Color& c0, const Color& c1)
    {
        Color outc;
        outc.r = (1-alpha) * c0.r +  alpha * c1.r;
        outc.b = (1-alpha) * c0.b+  alpha * c1.b;
        outc.g = (1-alpha) * c0.g +  alpha * c1.g;
        return outc;
    }

    /*
        * Get color from value in bwr palette
        * @param v value
        * @param vmin minimum value
        * @param vmax maximum value
        * @return color
    */
    Color static get_color_from_value(double v, double vmin, double vmax)
    {
        // if(v <  vmin){
        //     return Color(255, 255, 255);
        // }

        // clamp value within range
        v = v <= vmin ? vmin
            : v >= vmax ? vmax
            : v;

        const double alpha = (vmax <= vmin)
            ? 0.5                   // avoid divide by zero
            : (v - vmin) / (vmax - vmin);

        const Color blue(0, 0, 255);
        const Color white(255, 255, 255);
        const Color green(0, 255, 0);
        const Color yellow(255, 255, 0);
        const Color purple(160, 32, 240);
    
        const Color red(255, 0, 0); //maxcolor

        //palette br
        //return get_color_from_value(alpha, blue, red);

        //palette yb
        //return get_color_from_value(alpha, yellow, blue);
        

        //palette ypb
        // if (alpha < 0.5) {
        //     return get_color_from_value(alpha * 2, yellow, purple);
        // } else {
        //     return get_color_from_value(alpha * 2 - 1, purple, blue);
        // }
        
        //palette ywb
        if (alpha < 0.5) {
            return get_color_from_value(alpha * 2, yellow, white);
        } else {
            return get_color_from_value(alpha * 2 - 1, white, blue);
        }

        //palette bwr
        // if (alpha < 0.5) {
        //     return get_color_from_value(alpha * 2, blue, white);
        // } else {
        //     return get_color_from_value(alpha * 2 - 1, white, red);
        // }


        //palette bwygr
        // if (alpha < 0.25) {
        //     return get_color_from_value(alpha * 4, blue, white);
        // } else if (alpha < 0.5) {
        //     return get_color_from_value(alpha * 4 - 1, white, yellow);
        // } else if (alpha < 0.75) {
        //     return get_color_from_value(alpha * 4 - 2, yellow, green);
        // }
        // else {
        //     return get_color_from_value(alpha * 4 - 3, green, red);
        // }
    }

    // Overloading << operator for std::cout
    friend std::ostream& operator<<(std::ostream& os, const Color& color) {
        os <<"#" << std::hex << std::setw(2) << std::setfill('0') << color.r
           << std::hex << std::setw(2) << std::setfill('0') << color.g
           << std::hex << std::setw(2) << std::setfill('0') << color.b;
        return os;
    }
};

class Graphviz{
    public:
        string dotfilename="circle_graph.dot";
        string GRAPH_NAME = "Circle";
        ofstream* fout;
        ofstream file;

        Graphviz(){
            init();
        }
        ~Graphviz(){
            close();
        }
        
        void init(){
            file.open(dotfilename);
            if (!file.is_open()) {
                cerr << "Error opening file." << endl;
                return;
            }
            fout= &file;
        }

        void close(){
            file.close();
        }

        int getGridNo(int gx, int gy, int a){
            return gx + a*gy;
        }

        double getGridValue(int gx, int gy, int a, vector<double>& v){
            return v[gx + a*gy];
        }

        void draw_legend(double min, double max){
            ofstream file("legend.dot");

            if (!file.is_open()) {
                cerr << "Error opening file." << endl;
                return;
            }

            ofstream* fout= &file;

            (*fout)<< "digraph "<<GRAPH_NAME<< "{" << endl; 
            (*fout)<< "rankdir=LR;"<<endl;
            (*fout)<< "" <<endl;
            (*fout)<< "node [shape=square, width=0.2, height=0.2, fontsize=12];" <<endl;
            //(*fout)<< "edge [dir=none, arrowhead=empty, arrowtail=empty, arrowtail=normal, arrowhead=normal, labeldistance=2, labelangle=90, fontsize=8];" <<endl;
            (*fout)<< "edge [style=\"invis\"];" <<endl;
            
            
            (*fout)<< "" <<endl;
            
            //do 10 partition

            int num_partition = 6;

            //print [label="A", style="filled", fillcolor="blue", color="black"]; in for loop
            for (int i = 0; i <= num_partition; i++) {
                 const double alpha = min + (max - min)*i/num_partition*1.0;


                string colorhex = Color::get_color_from_value(alpha, min, max).rgbToHex();
                (*fout)<< "a" << i << " [label=\""<<  std::fixed << std::setprecision(1) << alpha <<"\", style=\"filled\", fillcolor=\""<< colorhex <<"\", color=\"black\" ];"<<endl;
            }

            for (int i = 0; i < num_partition; i++) {
                (*fout)<< "a" << i <<"-> "<< "a"<<i+1<<" [label=\"\"];"<<endl;
            }

            (*fout) << "}" << endl;
                
            file.close();
        }

        void draw_grid(int grid_size){

            (*fout)<< "graph "<<GRAPH_NAME<< "{" << endl; // graph Circle {
            (*fout)<< "node [shape=rectangle, width=1, height=1, style=filled];" << endl ; //node [shape=rectangle, width=1, height=1, style=filled];
            
            //print { rank=same; A1; A2; A3; A4; A5; } using for loop
            for (int j = grid_size-1; j >=0; j--) {
                (*fout) << "{ rank=same; ";
                for (int i = 0; i < grid_size; i++) {
                    (*fout) << "X" << i <<"Y"<< j << "; ";
                }
                (*fout) << "}" << endl;
            }



            for (int j = grid_size-1; j >=0; j--) {
                for (int i = 0; i < grid_size; i++) {
                    //string colorhex = "#FF0000";
                    string colorhex = Color::get_color_from_value(getGridNo(i,j,grid_size), 0, grid_size*grid_size-1).rgbToHex();
                    (*fout)<< "X" << i <<"Y"<< j << " [label=\"" << getGridNo(i,j,grid_size) <<  "\", fillcolor=\""<< colorhex <<"\"]"<<endl;
                }
            }
            (*fout) << "}" << endl;

        }

        void draw_grid_value(int grid_size, vector<double>& v_count, vector<double>& v_mean_alpha){

            (*fout)<< "graph "<<GRAPH_NAME<< "{" << endl; // graph Circle {
            (*fout)<< "node [shape=rectangle, width=1, height=1, style=filled];" << endl ; //node [shape=rectangle, width=1, height=1, style=filled];
            
            //print { rank=same; A1; A2; A3; A4; A5; } using for loop
            for (int j = grid_size-1; j >=0; j--) {
                (*fout) << "{ rank=same; ";
                for (int i = 0; i < grid_size; i++) {
                    (*fout) << "X" << j <<"Y"<< i << "; ";
                }
                (*fout) << "}" << endl;
            }


            //X0Y0--X1Y0--X2Y0--X3Y0--X4Y0--X5Y0--X6Y0--X7Y0
            //X0Y1--X1Y1--X2Y1--X3Y1--X4Y1--X5Y1--X6Y1--X7Y1
            for (int j = grid_size-1; j >=0; j--) {
                for (int i = 0; i < grid_size-1; i++) {
                    (*fout) << "X" << i <<"Y"<< j << "--";
                }
                (*fout) << "X" << grid_size-1 <<"Y"<< j ;
                (*fout)  << endl;
            }

            //to align ys node [group="samex0"]  // try to keep same Y value
            // X0Y0 X0Y1 X0Y2 X0Y3 
            for (int j = grid_size-1; j >=0; j--) {
                (*fout) << "node [group=\""<<"samex"<<j<<"\"]"<<endl;
                for (int i = 0; i < grid_size-1; i++) {
                    (*fout) << "X" << j <<"Y"<< i << " ";
                }
                (*fout) << "X" << j <<"Y"<< grid_size-1;
                (*fout)  << endl;
            }



            double min_v =  static_cast<double>(*(std::min_element(v_mean_alpha.begin(), v_mean_alpha.end())));
            double max_v =  static_cast<double>(*(std::max_element(v_mean_alpha.begin(), v_mean_alpha.end())));
            // cout<<min_v<<endl;
            // cout<<max_v<<endl;

            for (int j = grid_size-1; j >=0; j--) {
                for (int i = 0; i < grid_size; i++) {
                    //string colorhex = "#FF0000";
                    if(getGridValue(i,j,grid_size, v_count) ==0){
                        string colorhex = "#ffffff";
                        (*fout)<< "X" << i <<"Y"<< j << " [label=\"" << "" <<  "\", fillcolor=\""<< colorhex <<"\""<< ",pos=\"" << i<<","<<j  <<"!\"]"<<endl;
                    }else{
                        stringstream ss;
                        ss << std::fixed << std::setprecision(0) << getGridValue(i,j,grid_size, v_count) << ","<< std::fixed << std::setprecision(1) << getGridValue(i,j,grid_size, v_mean_alpha);
                        string strprint = ss.str();
                        //= to_string(getGridValue(i,j,grid_size, v_mean_alpha)) +" =N"+ to_string(getGridValue(i,j,grid_size, v_count));
                        string colorhex = Color::get_color_from_value(getGridValue(i,j,grid_size, v_mean_alpha), min_v, max_v).rgbToHex();
                        (*fout)<< "X" << i <<"Y"<< j << " [label=\"" <<  strprint <<  "\", fillcolor=\""<< colorhex <<"\""<< ",pos=\"" << i<<","<<j  <<"!\"]"<<endl;
                    }
                    
                }
            }
            (*fout) << "}" << endl;

        }


        void generateDotFileSample() {
            (*fout)<< "graph Circle {" << endl;
            (*fout) << "  node [shape=circle]" << endl;
            (*fout) << "  1 [pos=\"0,0!\"]" << endl; // Center of the circle
            const int numPoints = 12; // Number of points in the circle
            const double radius = 1.0; // Radius of the circle
            for (int i = 0; i < numPoints; ++i) {
                double angle = 2 * M_PI * i / numPoints;
                double x = radius * cos(angle);
                double y = radius * sin(angle);
                (*fout) << "  " << i+2 << " [pos=\"" << x << "," << y << "!\"]" << endl;
            }
            for (int i = 0; i < numPoints; ++i) {
                (*fout) << "  " << i+1 << " -- " << (i+2 > numPoints ? 1 : i+2) << endl;
            }
            (*fout) << "}" << endl;
        }
};

/*
int main() {
    // Create a Color object with RGB values
    cout<<Color::get_color_from_value(6, 0, 24);

    Graphviz gv;
    
    //gv.generateDotFileSample();
    gv.draw_grid(5);
    gv.draw_legend(0, 24);

    cout << "Dot file generated successfully." << endl;

    // Convert RGB to hexadecimal
    //std::string hexValue = red.rgbToHex();

    // Output the hexadecimal value
    //std::cout << "Hexadecimal representation: " << hexValue << std::endl;

    return 0;
}
*/