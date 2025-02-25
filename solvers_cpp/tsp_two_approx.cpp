#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <queue>
#include <algorithm>
#include <set>
#include <tuple>
#include <cmath>

class TwoApproximationSolver {
public:
    TwoApproximationSolver(const std::string& filename) {
        std::ifstream file(filename);
        std::string line;
        bool read_distances = false;

        while (std::getline(file, line)) {
            std::istringstream iss(line);

            if (line.find("EOF") != std::string::npos) {
                break;
            } else if (line.find("EDGE_WEIGHT_TYPE") != std::string::npos) {
                if (line.find("EXPLICIT") != std::string::npos) {
                    explicit_weights = true;
                } 
                else if (line.find("CEIL") != std::string::npos) {
                    explicit_weights = false;
                    ceil = true; 
                }
                else {
                    explicit_weights = false;
                    ceil = false; 
                }
                    
            } else if (line.find("EDGE_WEIGHT_SECTION") != std::string::npos || line.find("NODE_COORD_SECTION") != std::string::npos) {
                read_distances = true;
            } else if (explicit_weights && read_distances) {
                std::vector<int> row;
                int value;
                while (iss >> value) {
                    row.push_back(value);
                }
                distance_matrix.push_back(row);
            } else if (!explicit_weights && read_distances) {
                float index, x, y;
                iss >> index >> x >> y;
                euc_coordinates.push_back({x, y});
            }
        }

        if (explicit_weights) {
            n = distance_matrix.size();
        } else {
            n = euc_coordinates.size();
        }
    }

    int get_cost(int i, int j) {
        if (explicit_weights) {
            return distance_matrix[i][j];
        } else {
            float x1 = euc_coordinates[i].first;
            float y1 = euc_coordinates[i].second;
            float x2 = euc_coordinates[j].first;
            float y2 = euc_coordinates[j].second;

            if (ceil) {
                return std::ceil(std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2))); 
            } 
            else {
                return std::round(std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)));
            }
            
        }
    }

    void prim_mst() {
        mst_edges.clear();
        total_weight = 0;
        std::vector<bool> visited(n, false);
        std::priority_queue<std::tuple<int, int, int>, std::vector<std::tuple<int, int, int>>, std::greater<>> min_heap;
        min_heap.push({0, 0, -1});  // (weight, vertex, parent)

        while (mst_edges.size() < n - 1) {
            auto [weight, u, parent] = min_heap.top();
            min_heap.pop();

            if (visited[u]) {
                continue;
            }

            visited[u] = true;
            if (parent != -1) {
                mst_edges.push_back({parent, u, weight});
                total_weight += weight;
            }

            for (int v = 0; v < n; ++v) {
                if (!visited[v] && get_cost(u, v) > 0) {
                    min_heap.push({get_cost(u, v), v, u});
                }
            }
        }
    }

    std::vector<int> dfs_mst(int start) {
        std::vector<bool> visited(n, false);
        std::vector<int> tour;
        std::vector<int> stack = {start};

        while (!stack.empty()) {
            int u = stack.back();
            stack.pop_back();
            if (visited[u]) {
                continue;
            }

            visited[u] = true;
            tour.push_back(u);

            for (const auto& edge : mst_edges) {
                if (std::get<0>(edge) == u && !visited[std::get<1>(edge)]) {
                    stack.push_back(std::get<1>(edge));
                } else if (std::get<1>(edge) == u && !visited[std::get<0>(edge)]) {
                    stack.push_back(std::get<0>(edge));
                }
            }
        }

        return tour;
    }

    std::pair<int, std::vector<int>> two_approx(int start) {
        prim_mst();
        std::vector<int> tour = dfs_mst(start);

        std::vector<int> unique_tour;
        std::set<int> seen;
        int cost = 0;

        for (size_t i = 0; i < tour.size(); ++i) {
            int vertex = tour[i];
            if (seen.find(vertex) == seen.end()) {
                unique_tour.push_back(vertex);
                seen.insert(vertex);

                if (unique_tour.size() > 1) {
                    int prev_vertex = unique_tour[unique_tour.size() - 2];
                    cost += get_cost(prev_vertex, vertex);
                }
            }
        }

        cost += get_cost(unique_tour.back(), unique_tour[0]);

        return {cost, unique_tour};
    }

private:
    std::vector<std::vector<int>> distance_matrix;
    std::vector<std::pair<float, float>> euc_coordinates;
    std::vector<std::tuple<int, int, int>> mst_edges;
    int n;
    int total_weight;
    bool explicit_weights;
    bool ceil;
};

void write_output(const std::string& filename, int cost, const std::vector<int>& path) {

    std::string output_filename = "solver_2approx_" + filename + ".txt";
    std::ofstream output_file(output_filename);

    if (output_file.is_open()) {
        output_file << cost << "\n[";

        for (size_t i = 0; i < path.size(); ++i) {
            output_file << path[i];
            if (i < path.size() - 1) output_file << ", ";
        }

        output_file << "]\n";
    } else {
        std::cerr << "Failed to open the file.\n";
    }

    output_file.close();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    TwoApproximationSolver solver(filename);
    auto result = solver.two_approx(0);

    // prepare output
    int cost = result.first;
    std::vector<int> path = result.second;
    write_output(filename, cost, path);

    return 0;
}