//
//  donuts.cpp
//  p4-donuts
// IDENTIFIER  = 8729BF7B2234189B8DF8A6DD31770D2B18569C27
//

#include "donuts.hpp"
#include <getopt.h>
#include <string>
#include <vector>
#include <deque>
#include <queue>
#include <cmath>
#include <iomanip>
#include "xcode_redirect.hpp"

using namespace std;

void Donuts::getMode(int argc, char * argv[]) {
    string mode = "";
    opterr = false;
    int choice;
    int index = 0;
    option long_options[] = {{"mode", required_argument, nullptr, 'm'},
        {"help", no_argument, nullptr, 'h'},
        { nullptr, 0, nullptr, '\0' }};
    /*
     Options with required_argument (print) need a colon after the
     char, options with no_argument do not (help).
     */
    while ((choice = getopt_long(argc, argv, "m:h", long_options, &index)) != -1) {
        switch (choice) {
            case 'm':
                mode = optarg;
                if (mode == "MST") {
                    MST = true;
                }
                else if (mode == "FASTTSP") {
                    FASTTSP = true;
                }
                else if (mode == "OPTTSP") {
                    OPTTSP = true;
                }
                if (!MST && !FASTTSP && !OPTTSP) {
                    cerr << "Invalid mode" << endl;
                }
                break;
                
            case 'h':
                cout << "help option" << endl;
                exit(0);

            default:
                if (!MST && !FASTTSP && !OPTTSP) {
                    cerr << "No mode specified" << endl;
                }
                else {
                    cerr << "Invalid command line option" << endl;
                }
                exit(1);
        } // switch
    } // while
} // get_options

double Donuts::distance(Coordinate cor1, Coordinate cor2) {
    // so as helper, you give it two coordinates, it gives you a double. if those coordinates are USA/CAN or CAN/USA say infinity, but anything else, give back a number (double).
    if (MST && ((cor1.category == Category::USA && cor2.category == Category::Canada) || (cor1.category == Category::Canada && cor2.category == Category::USA))) {
        return numeric_limits<double>::infinity();
    }
    double x_sub = (cor1.x) - (cor2.x);
    double y_sub = (cor1.y) - (cor2.y);
    x_sub *= x_sub;
    y_sub *= y_sub;
    double dist = sqrt(x_sub + y_sub);
    return dist;
}

void Donuts::read_data() {
    int x, y;
    bool border_exist = 0;
    cin >> num_vertice;
    Coordinate newCoord;
    size_t i = 0;
    while (i < num_vertice) {
        cin >> x >> y;
        newCoord.x = x;
        newCoord.y = y;
        if (x == 0 && y >= 0) { // vertical border
                newCoord.category = Category::Border;
                border_exist = 1;
        } // if border
        else if (x > 0 && y > 0) { // CANADA (top right)
                newCoord.category = Category::Canada;
        }   // if canada
        else if (x >= 0 && y == 0) { // horizontal border
                newCoord.category = Category::Border;
                border_exist = 1;
        } // if border
        else {
            newCoord.category = Category::USA;
        } // USA
        coordinates.push_back(newCoord);
        i++;
    }
    if (MST && !border_exist) {
        cerr << "Cannot construct MST" << endl;
        exit(1);
    }
} // read_data()

void Donuts::MST_mode() {
    if (MST) {
        size_t index = 0;
        double best_dist = numeric_limits<double>::infinity();
        double cur_dist = 0;
        coordinates[0].distance = 0;
        // num_vertice = outie size
        // dont use count as an index, its just how many times we go thru it
        for (size_t count = 0; count < static_cast<size_t>(num_vertice); ++count) {
                for (size_t v = 0; v < num_vertice; ++v) {
                    if (!coordinates[v].visited) {
                        // step 1: find smallest false
                        // select the vertex v having the smallest tentative distance.
                        cur_dist = coordinates[v].distance;
                        if (cur_dist <= best_dist) {
                            best_dist = cur_dist;
                            index = v;
                        }
                    }
                }
            // step 2: mark as true
            coordinates[index].visited = 1;
            for (size_t w = 0; w < static_cast<size_t>(num_vertice); ++w) {
            // step 3: update the false neighbors
                if (!coordinates[w].visited) {
                    if (coordinates[w].distance > distance(coordinates[index], coordinates[w])) {
                        coordinates[w].distance = distance(coordinates[index], coordinates[w]);
                        coordinates[w].prevCoord = index;
                    }
                }
            }
            weight += best_dist;
            best_dist = numeric_limits<double>::infinity();
        }
        // print
        cout << weight << endl;
        for (size_t i = 1; i < num_vertice; ++i) {
            if (i <= coordinates[i].prevCoord) {
                cout << i << " " << coordinates[i].prevCoord << endl;
            }
            else {
                cout << coordinates[i].prevCoord << " " << i << endl;
            }
        }
    }
} // mst

double Donuts::MST_ptC(size_t permLength, ptcOPTTSP &ptc) {
    size_t index = 0;
    double best_dist = numeric_limits<double>::infinity();
    double cur_dist = 0;
    coordinates[permLength].distance = 0;
    coordinates[permLength].visited = 0;
    
    // num_vertice = outie size
    // dont use count as an index, its just how many times we go thru it
    for (size_t count = permLength; count < coordinates.size(); ++count) {
            for (size_t v = permLength; v < ptc.current_path.size(); ++v) {
                if (!coordinates[v].visited) {
                    // step 1: find smallest false
                    // select the vertex v having the smallest tentative distance.
                    cur_dist = coordinates[v].distance;
                    if (cur_dist <= best_dist) {
                        best_dist = cur_dist;
                        index = v;
                    }
                }
            }
        // step 2: mark as true
        coordinates[index].visited = 1;
        for (size_t w = permLength; w < coordinates.size(); ++w) {
        // step 3: update the false neighbors
            if (!coordinates[w].visited) {
                if (coordinates[w].distance > distance(coordinates[index], coordinates[w])) {
                    coordinates[w].distance = distance(coordinates[index], coordinates[w]);
                    coordinates[w].prevCoord = index;
                }
            }
        }
        weight += best_dist;
        best_dist = numeric_limits<double>::infinity();
    }
    return weight;
} // mst

double Donuts::FASTTSP_ptC(ptcOPTTSP &ptc) {
    // arbitrary insertion

        // @6314
        // A -> B -> A
        // Step 1: Initialize a partial tour with a vertex i, chosen arbitrarily (you can just start with the first vertex available).
        // Step 2: Choose another arbitrary vertex j and set the initial partial tour to i → j → i.
    ptc.best_path.push_back(0);
    ptc.best_path.push_back(1);
    ptc.best_path.push_back(0);
        // 3. Arbitrarily select a vertex k that is currently not in the partial tour
        // int k = 2
        // When you are calculating the costs for ik, kj, and ij edges, make sure you are doing a nested subscript when indexing into your vector of vertices using i and j. What I mean is i and j represent indices to your path vector.
        // path vector stores the actual indices of the vertices you are calculating distance for, so you want to index into your vector of vertices like vertices[path[i]] or vertices[path[j]]. Keep vertices[k] because k is the actual index of the vertex you are trying to insert between i and j.
        // @6411
        for (size_t k = 2; k < num_vertice; ++k) {
            double best_cost = numeric_limits<double>::infinity();
            // edges
            auto insertLoc = ptc.best_path.begin(); // add to this to get where we wanna insert
            for (auto i = ptc.best_path.begin(); i != ptc.best_path.end() - 1; ++i) {
                auto j_index = i + 1;
                //          double c_ik + c_kj - c_ij
                // undef behavior
                //                double edge_cost = distance(coordinates[path[i]], coordinates[k]) + distance(coordinates[k], coordinates[path[j_index]]) - distance(coordinates[path[i]], coordinates[path[j_index]]);
                double edge_cost = distance(coordinates[*i], coordinates[k]) + distance(coordinates[k], coordinates[*j_index]) - distance(coordinates[*i], coordinates[*j_index]);
                
                if (edge_cost < best_cost) {
                    best_cost = edge_cost;
                    insertLoc = i + 1;
                    // j_index + 1 ?
                }
                // Calculating the cost of edges ik, kj, and ij goes inside the inner loop. The outer loop goes until the amount of vertices, and the inner loop goes until the current path size you are calculating in part B minus 1 because you are using index i and index i+1(which is j).
            }
            // insert takes iterator pos
            //            path.insert(<#const_iterator position#>, <#const_reference x#>)
            ptc.best_path.insert(insertLoc, k);
            //            path.insert(path.begin() + j_index, k);
        }
        // @6413
        // go thru path and print
        for (size_t i = 0; i < ptc.best_path.size() - 1; ++i) {
            weight += distance(coordinates[ptc.best_path[i]], coordinates[ptc.best_path[i + 1]]);
        }
        cout << weight << endl;
        for (size_t i = ptc.best_path.size() - 1; i > 0; --i) {
            cout << ptc.best_path[i] << " ";
        }
        cout << endl;
    
    return weight;
}

void Donuts::FASTTSP_mode() {
    // arbitrary insertion
    if (FASTTSP) {
        // @6314
        // A -> B -> A
        // Step 1: Initialize a partial tour with a vertex i, chosen arbitrarily (you can just start with the first vertex available).
        // Step 2: Choose another arbitrary vertex j and set the initial partial tour to i → j → i.
        path.push_back(0);
        path.push_back(1);
        path.push_back(0);
        // 3. Arbitrarily select a vertex k that is currently not in the partial tour
        // int k = 2
        // When you are calculating the costs for ik, kj, and ij edges, make sure you are doing a nested subscript when indexing into your vector of vertices using i and j. What I mean is i and j represent indices to your path vector.
        // path vector stores the actual indices of the vertices you are calculating distance for, so you want to index into your vector of vertices like vertices[path[i]] or vertices[path[j]]. Keep vertices[k] because k is the actual index of the vertex you are trying to insert between i and j.
        // @6411
        for (size_t k = 2; k < num_vertice; ++k) {
            double best_cost = numeric_limits<double>::infinity();
            // edges
            auto insertLoc = path.begin(); // add to this to get where we wanna insert
            for (auto i = path.begin(); i != path.end()-1; ++i) {
                auto j_index = i + 1;
                //          double c_ik + c_kj - c_ij
                // undef behavior
//                double edge_cost = distance(coordinates[path[i]], coordinates[k]) + distance(coordinates[k], coordinates[path[j_index]]) - distance(coordinates[path[i]], coordinates[path[j_index]]);
            double edge_cost = distance(coordinates[*i], coordinates[k]) + distance(coordinates[k], coordinates[*j_index]) - distance(coordinates[*i], coordinates[*j_index]);
                
                if (edge_cost < best_cost) {
                    best_cost = edge_cost;
                    insertLoc = i + 1;
                    // j_index + 1 ?
                }
    // Calculating the cost of edges ik, kj, and ij goes inside the inner loop. The outer loop goes until the amount of vertices, and the inner loop goes until the current path size you are calculating in part B minus 1 because you are using index i and index i+1(which is j).
            }
            // insert takes iterator pos
//            path.insert(<#const_iterator position#>, <#const_reference x#>)
            path.insert(insertLoc, k);
//            path.insert(path.begin() + j_index, k);
        }
    // @6413
        // go thru path and print
        for (size_t i = 0; i < path.size() - 1; ++i) {
            weight += distance(coordinates[path[i]], coordinates[path[i + static_cast<size_t>(1)]]);
        }
        cout << weight << endl;
        for (size_t i = path.size() - 1; i > 0; --i) {
            cout << path[i] << " ";
        }
        cout << endl;
    }
} //  fasttsp mode

// use in Part C opttsp
// Another thing need to do is speed things up a bit. This is the promising function. It is supposed to answer the question: does this look like a good way to start?
// Remember: promising
//arm1 = 0 to closest unvisited, arm2 = last thing in fixed portion path to unvisited
//If we do the estimate, we have to do (curCost + arm1 + arm2 + MST) < bestCost.
//This is what promising answers ^ this question.
// arm1 and arm2 are the cheapest connections from your current path to the MST. TotalEst is the total estimate (your path + arm 1 + arm 2 + mst). Promise is a Boolean that relates this total to whether or not the solution is promising. Think about when this is true or false
bool Donuts::promising(size_t permLength, ptcOPTTSP &ptc) {
    // implement
    // find out if this path is promising
    // : take their part A code, copying it to part C maybe as member fxn, and modifying it so that when it modifies distance it uses a helper that doesnt care about border/Can/USA and makes it able to do an MST of some of the vertices and not all of the vertices.
    double current_min = 0;
    double start_arm = numeric_limits<double>::infinity();
    double end_arm = numeric_limits<double>::infinity();
    double mst_ptC = MST_ptC(permLength, ptc); //  neow we do mst of some of the vertices not all
    
    for(size_t i = permLength; i < ptc.current_path.size(); i ++) {
        current_min = ptc.distances[ptc.current_path[0]][ptc.current_path[i]];
        start_arm = min(start_arm, current_min);
        
        current_min = ptc.distances[ptc.current_path[permLength-1]][ptc.current_path[i]];
        end_arm = min(end_arm, current_min);
    }
    
    lower_bound = ptc.current_path_cost + start_arm + end_arm + mst_ptC;
    
    if (lower_bound < ptc.best_cost) {
        return true;
    }
    else {
        return false;
    }
    return 0;
}

// for partC
void Donuts::genPerms(size_t permLength, ptcOPTTSP &ptc) {
  if (permLength == ptc.current_path.size()) {
    
      double total_cost = ptc.current_path_cost + distance(coordinates[ptc.current_path[ptc.current_path.size() - 1]], coordinates[ptc.current_path[0]]);
      
      if (ptc.current_path_cost < ptc.best_cost) {
          ptc.best_cost = total_cost;
          ptc.best_path = ptc.current_path;
      }
    return;
  }  // if ..complete path

  if (!promising(permLength, ptc)) {
    return;
  }  // if ..not promising

  for (size_t i = permLength; i < ptc.current_path.size(); ++i) {
      // Got to add the code that keeps the currentPathCost updated and with each recursive call to genPerms you should be adding something to the total and after the recursive call you should take it away from the total. Think of preconditions, post conditions,
    swap(ptc.current_path[permLength], ptc.current_path[i]);
      ptc.current_path_cost += ptc.distances[ptc.current_path[permLength]][ptc.current_path[permLength - 1]];
    genPerms(permLength + 1, ptc);
      ptc.current_path_cost -= ptc.distances[ptc.current_path[permLength]][ptc.current_path[permLength - 1]];
    swap(ptc.current_path[permLength], ptc.current_path[i]);
  }  // for ..unpermuted elements
}  // genPerms()


// In part C, we call part B to get our best we’ve ever seen and call part A to get our MST.
void Donuts::OPTTSP_mode() {
    //  you could determine that some branches of the search cannot lead to optimal solutions. For example, you could compute lower bounds on the length of any full tour that can be found in a given branch. If such a lower bound exceeds the cost of a full solution you have found previously, you can skip this branch as hopeless
    // Your logic sounds reasonable to me, make sure though that you are keeping track of the minimum distance between the start of the partial tour and your MST (arm1) AND the minimum distance between the end of the partial tour and your MST (arm2)
    // The first, but we also need to consider the MST calculates from the unvisited nodes for the lower bound in addition to the two weights added from our connecting arms (A to some node and C to some node).
    
    // declare ptcOPTTSP --> pass it throguh
    if (OPTTSP) {
        ptcOPTTSP ptC;
        for(size_t i = 0; i < coordinates.size(); i ++) {
            vector<double> temp;
            for (size_t j = 0; j < coordinates.size(); j ++) {
                double dist = distance(coordinates[i], coordinates[j]);
                if( i == j) temp.push_back(0);
                else temp.push_back(dist);
            }
            
            ptC.distances.push_back(temp);
            
        }
        ptC.best_cost = FASTTSP_ptC(ptC);
        ptC.current_path = ptC.best_path;
        genPerms(1, ptC);

    }
}

int main(int argc, char *argv[]) {
    cout << std::setprecision(2);
    cout << std::fixed;
    
    ios_base::sync_with_stdio(false);
    xcode_redirect(argc, argv);
    
    Donuts donutPath;
    
    donutPath.getMode(argc, argv);
    donutPath.read_data();
    donutPath.MST_mode();
    donutPath.FASTTSP_mode();
    donutPath.OPTTSP_mode();
    
    return 0;
}
