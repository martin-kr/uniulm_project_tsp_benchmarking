#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <limits>
#include <numeric>
#include <cmath>
#include <omp.h>

class TwoOptSolver {
public:
    TwoOptSolver(const std::string& filename) {
        std::ifstream file(filename);
        std::string line;
        bool read_distances = false;

        while (std::getline(file, line)) {
            std::istringstream iss(line);

            if (line.find("EOF") != std::string::npos) {
                break;
            } else if (line.find("NODE_COORD_SECTION") != std::string::npos) {
                read_distances = true;
            } else if (read_distances) {
                float index, x, y;
                iss >> index >> x >> y;
                euc_coordinates.push_back({x, y});
            }
        }

        n = euc_coordinates.size(); 
    }

    int get_cost(int i, int j) {
        float x1 = euc_coordinates[i].first;
        float y1 = euc_coordinates[i].second;
        float x2 = euc_coordinates[j].first;
        float y2 = euc_coordinates[j].second;
        return std::round(std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2))); 
    }

    int calc_path_length() {
        int path_cost = 0;

        for (int i = 0; i < n - 1; ++i) {
            int vertex_i = best_route[i];
            int vertex_j = best_route[i + 1];
            int cost_ij = get_cost(vertex_i, vertex_j);
            path_cost += cost_ij;
        }

        // n -> 0
        int vertex_i = best_route[n - 1];
        int vertex_j = best_route[0];
        int cost_ij = get_cost(vertex_i, vertex_j);
        path_cost += cost_ij;

        return path_cost;
    }

    void swap_edges(int index_i, int index_j) {
        ++index_i;

        while (index_i < index_j) {
            std::swap(best_route[index_i], best_route[index_j]);
            ++index_i;
            --index_j;
        }
    }

    std::pair<int, std::vector<int>> two_opt() {
        best_route.resize(n);
        std::iota(best_route.begin(), best_route.end(), 0);
        best_length = calc_path_length();
        bool improvement_found = true;

        while (improvement_found) {
            improvement_found = false;
            bool break_outer_loop = false;

            #pragma omp parallel for schedule(dynamic) shared(improvement_found, break_outer_loop)
            for (int i = 0; i < n - 1; ++i) {
                if (break_outer_loop) {
                    continue;
                }

                for (int j = i + 2; j < n; ++j) {
                    int j_1 = (j + 1) % n;

                    int vertex_i = best_route[i];
                    int vertex_i1 = best_route[i + 1];
                    int vertex_j = best_route[j];
                    int vertex_j1 = best_route[j_1];

                    int r_i_i1 = get_cost(vertex_i, vertex_i1);
                    int r_j_j1 = get_cost(vertex_j, vertex_j1);
                    int r_i_j = get_cost(vertex_i, vertex_j);
                    int r_i1_j1 = get_cost(vertex_i1, vertex_j1);

                    int minus = r_i_i1 + r_j_j1;
                    int plus = r_i_j + r_i1_j1;

                    int len_delta = plus - minus;

                    if (len_delta < 0 && !break_outer_loop) {
                        #pragma omp critical
                        {
			    if (!break_outer_loop) {
                                improvement_found = true;
                                best_length += len_delta;
                                swap_edges(i, j);
                                break_outer_loop = true;
			    }
                        }
                    }
                }
            }
        }

        return {best_length, best_route};
    }

private:
    std::vector<std::pair<float, float>> euc_coordinates;
    std::vector<int> best_route;
    int n;
    int best_length;
};

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    TwoOptSolver solver(filename);
    auto result = solver.two_opt();

    // prepare output
    int cost = result.first;
    std::vector<int> path = result.second;

    // Construct the output filename
    std::string output_filename = "tc4_mt_" + filename + ".txt";
    
    // Open the output file stream
    std::ofstream output_file(output_filename);

    if (output_file.is_open()) {
        // Write the cost and path in the shortened format
        output_file << cost << "\n[";

        for (size_t i = 0; i < path.size(); ++i) {
            output_file << path[i];
            if (i < path.size() - 1) output_file << ", ";
        }

        output_file << "]\n";
    } else {
        std::cerr << "Failed to open the file.\n";
    }

    // Close the file
    output_file.close();
}
