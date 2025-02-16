#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
#include <cstdlib>
#include <ctime>

class NearestNeighborSolver {
    public:
        NearestNeighborSolver(const std::string& filename) {
            std::ifstream file(filename);
            std::string line;
            bool read_distances = false;

            while (std::getline(file, line)) {
                line.erase(std::remove(line.begin(), line.end(), ':'), line.end());
                std::istringstream iss(line);

                if (line.find("EOF") != std::string::npos) {
                    break;
                } else if (line.find("EDGE_WEIGHT_SECTION") != std::string::npos) {
                    read_distances = true;
                } else if (read_distances) {
                    std::vector<int> row;
                    int value;
                    while (iss >> value) {
                        row.push_back(value);
                    }
                    distance_matrix.push_back(row);
                }
            }

            n = distance_matrix.size();
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
                    if (!visited[j] && distance_matrix[last_visited][j] > 0 && distance_matrix[last_visited][j] < min_dist) {
                        nearest = j;
                        min_dist = distance_matrix[last_visited][j];
                    }
                }

                path.push_back(nearest);
                visited[nearest] = true;
                total_cost += min_dist;
            }

            total_cost += distance_matrix[path.back()][path[0]];  // Return to start
            return {total_cost, path};
        }

    private:
        std::vector<std::vector<int>> distance_matrix;
        int n;
};

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    srand(time(0));  // Seed for random number generation
    NearestNeighborSolver solver(filename);
    auto result = solver.nearest_neighbor();
    std::cout << "Total cost: " << result.first << "\nPath: ";
    for (int city : result.second) {
        std::cout << city << " ";
    }
    std::cout << std::endl;
    return 0;
}