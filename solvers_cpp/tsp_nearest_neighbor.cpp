#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <ctime>

class NearestNeighborSolver {
public:
    NearestNeighborSolver(const std::string& filename) {
        std::ifstream file(filename);
        std::string line;
        bool read_distances = false;
        bool explicit_weights = true;

        while (std::getline(file, line)) {
            std::istringstream iss(line);

            if (line.find("EOF") != std::string::npos) {
                break;
            } else if (line.find("EDGE_WEIGHT_TYPE") != std::string::npos) {
                if (line.find("EXPLICIT") != std::string::npos) {
                    explicit_weights = true;
                } else {
                    explicit_weights = false;
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
                int index, x, y;
                iss >> index >> x >> y;
                euc_coordinates.push_back({x, y});
            }
        }

        if (explicit_weights) {
            n = distance_matrix.size();
            this->explicit_weights = true;
        } else {
            n = euc_coordinates.size();
            this->explicit_weights = false;
        }
    }

    int get_cost(int i, int j) {
        if (explicit_weights) {
            return distance_matrix[i][j];
        } else {
            int x1 = euc_coordinates[i].first;
            int y1 = euc_coordinates[i].second;
            int x2 = euc_coordinates[j].first;
            int y2 = euc_coordinates[j].second;
            return std::ceil(std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)));
        }
    }

    std::pair<int, std::vector<int>> nearest_neighbor(int start = 0) {
        std::vector<bool> visited(n, false);
        std::vector<int> path;
        path.push_back(start);
        visited[start] = true;
        int total_cost = 0;

        for (int i = 0; i < n - 1; ++i) {
            int last_visited = path.back();
            int nearest = -1;
            int min_dist = std::numeric_limits<int>::max();

            for (int j = 0; j < n; ++j) {
                if (!visited[j] && get_cost(last_visited, j) > 0 && get_cost(last_visited, j) < min_dist) {
                    nearest = j;
                    min_dist = get_cost(last_visited, j);
                }
            }

            path.push_back(nearest);
            visited[nearest] = true;
            total_cost += min_dist;
        }

        total_cost += get_cost(path.back(), path[0]);  // Return to start
        return {total_cost, path};
    }

private:
    std::vector<std::vector<int>> distance_matrix;
    std::vector<std::pair<int, int>> euc_coordinates;
    int n;
    bool explicit_weights;
};

void write_output(const std::string& filename, int cost, const std::vector<int>& path) {
    // Construct the output filename
    std::string output_filename = "christofides_" + filename + ".txt";
    
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

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    srand(time(0));  // Seed for random number generation
    NearestNeighborSolver solver(filename);
    auto result = solver.nearest_neighbor();

    // prepare output
    int cost = result.first;
    std::vector<int> path = result.second;
    write_output(filename, cost, path);

    return 0;
}