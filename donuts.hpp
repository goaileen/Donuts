//
//  donuts.hpp
//  p4-donuts
// IDENTIFIER  = 8729BF7B2234189B8DF8A6DD31770D2B18569C27
//

#ifndef donuts_hpp
#define donuts_hpp

#include <stdio.h>
#include <getopt.h>
#include <vector>
#include <queue>
#include <deque>
#include <limits>
#include "xcode_redirect.hpp"

#endif /* donuts_hpp */

using namespace std;

class Donuts {
private:
    enum class Category:char { USA, Canada, Border };
    // x-y coordinates = vertice, how they are connected = edges
    struct Coordinate {
        int x;
        int y;
        size_t prevCoord;
        double distance = numeric_limits<double>::infinity();
        bool visited = false;
        Category category;
    };
    double weight = 0;
    // For the MST mode, you should print the total weight of the MST you generate by itself on a line; this weight is the sum of the weights of all edges in your MST (in terms of Euclidean distance). You should then print all edges in the MST. All output should be printed to standard output (cout).
    bool MST = 0;
    bool FASTTSP = 0;
    bool OPTTSP = 0;
    size_t num_vertice = 0;

    vector <size_t> path; // when done producing TSP, have to print total and vertices names and their names are their indicies. so when print those, its good if u had a vector that had those to begin with.
    vector <Coordinate> coordinates; // all coords
    double lower_bound;
    
    struct ptcOPTTSP {
        // distance matrix member var
        vector<vector<double>> distances;
        vector<size_t> current_path; // partial path
        vector<size_t> best_path;
        double current_path_cost;
        double best_cost;
    };
    
public:
    void getMode(int argc, char * argv[]);
    void read_data();
    double distance(Coordinate cor1, Coordinate cor2);
    void MST_mode();
    double MST_ptC(size_t permLength, ptcOPTTSP &ptc);
    void FASTTSP_mode();
    double FASTTSP_ptC(ptcOPTTSP &ptc);
    void OPTTSP_mode();
    void genPerms(size_t permLength, ptcOPTTSP &ptc);
    bool promising(size_t permLength, ptcOPTTSP &ptc);
    };
