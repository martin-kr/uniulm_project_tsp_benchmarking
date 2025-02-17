#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <queue>
#include <set>
#include <unordered_map>
#include <deque>
#include <tuple>
#include <limits>
#include <numeric>
#include <iterator>

class ChristofidesSolver {
public:
    ChristofidesSolver(const std::string& filename) {
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

    void prim_mst() {
        mst_edges.clear();
        total_weight = 0;
        std::vector<bool> visited(n, false);
        std::priority_queue<std::tuple<int, int, int>, std::vector<std::tuple<int, int, int>>, std::greater<>> min_heap;
        min_heap.push(std::make_tuple(0, 0, -1));  // (weight, vertex, parent)

        while (mst_edges.size() < n - 1) {
            auto [weight, u, parent] = min_heap.top();
            min_heap.pop();

            if (visited[u]) {
                continue;
            }

            visited[u] = true;
            if (parent != -1) {
                mst_edges.push_back(std::make_tuple(parent, u, weight));
                total_weight += weight;
            }

            for (int v = 0; v < n; ++v) {
                if (!visited[v] && distance_matrix[u][v] > 0) {
                    min_heap.push(std::make_tuple(distance_matrix[u][v], v, u));
                }
            }
        }
    }

    void find_odd_degree_vertices() {
        std::vector<int> degree(n, 0);

        for (const auto& [u, v, _] : mst_edges) {
            degree[u]++;
            degree[v]++;
        }

        odd_vertices.clear();
        for (int v = 0; v < n; ++v) {
            if (degree[v] % 2 == 1) {
                odd_vertices.push_back(v);
            }
        }
    }

    void minimum_weight_perfect_matching() {
        std::vector<std::tuple<int, int, int>> edges;

        for (size_t i = 0; i < odd_vertices.size(); ++i) {
            for (size_t j = i + 1; j < odd_vertices.size(); ++j) {
                int u = odd_vertices[i];
                int v = odd_vertices[j];
                int weight = distance_matrix[u][v];
                edges.push_back(std::make_tuple(weight, u, v));
            }
        }

        std::sort(edges.begin(), edges.end());

        std::set<int> matched;
        matching_edges.clear();

        for (const auto& [weight, u, v] : edges) {
            if (matched.find(u) == matched.end() && matched.find(v) == matched.end()) {
                matching_edges.push_back(std::make_tuple(u, v, weight));
                matched.insert(u);
                matched.insert(v);
            }
        }
    }

    void find_euler_circuit() {
        mst_mwpm_graph = mst_edges;
        mst_mwpm_graph.insert(mst_mwpm_graph.end(), matching_edges.begin(), matching_edges.end());

        std::unordered_map<int, std::deque<int>> adj;

        for (const auto& [u, v, _] : mst_mwpm_graph) {
            adj[u].push_back(v);
            adj[v].push_back(u);
        }

        for (const auto& [node, neighbors] : adj) {
            if (neighbors.size() % 2 != 0) {
                std::cerr << "Vertex " << node << " has an odd degree, so no Euler circuit exists." << std::endl;
                return;
            }
        }

        euler_circuit.clear();
        std::deque<int> current_path = {adj.begin()->first};

        while (!current_path.empty()) {
            int u = current_path.back();

            if (!adj[u].empty()) {
                int v = adj[u].front();
                adj[u].pop_front();
                adj[v].erase(std::remove(adj[v].begin(), adj[v].end(), u), adj[v].end());
                current_path.push_back(v);
            } else {
                euler_circuit.push_back(current_path.back());
                current_path.pop_back();
            }
        }

        std::reverse(euler_circuit.begin(), euler_circuit.end());
    }

    std::pair<int, std::vector<int>> christofides() {
        prim_mst();
        find_odd_degree_vertices();
        minimum_weight_perfect_matching();
        find_euler_circuit();

        std::vector<int> unique_tour;
        std::set<int> seen;
        int cost = 0;

        for (int vertex : euler_circuit) {
            if (seen.find(vertex) == seen.end()) {
                unique_tour.push_back(vertex);
                seen.insert(vertex);

                if (unique_tour.size() > 1) {
                    int prev_vertex = unique_tour[unique_tour.size() - 2];
                    cost += distance_matrix[prev_vertex][vertex];
                }
            }
        }

        cost += distance_matrix[unique_tour.back()][unique_tour[0]];

        return {cost, unique_tour};
    }

private:
    std::vector<std::vector<int>> distance_matrix;
    std::vector<std::tuple<int, int, int>> mst_edges;
    std::vector<int> odd_vertices;
    std::vector<std::tuple<int, int, int>> matching_edges;
    std::vector<std::tuple<int, int, int>> mst_mwpm_graph;
    std::vector<int> euler_circuit;
    int n;
    int total_weight;
};

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    ChristofidesSolver solver(filename);
    auto result = solver.christofides();
    std::cout << "Total cost: " << result.first << "\nPath: ";
    for (int city : result.second) {
        std::cout << city << " ";
    }
    std::cout << std::endl;
    return 0;
}


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <queue>
#include <set>
#include <unordered_map>
#include <deque>
#include <tuple>
#include <limits>
#include <cmath>

class ChristofidesSolver {
public:
    ChristofidesSolver(const std::string& filename) {
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

    void find_odd_degree_vertices() {
        std::vector<int> degree(n, 0);

        for (const auto& [u, v, _] : mst_edges) {
            degree[u]++;
            degree[v]++;
        }

        odd_vertices.clear();
        for (int v = 0; v < n; ++v) {
            if (degree[v] % 2 == 1) {
                odd_vertices.push_back(v);
            }
        }
    }

    void minimum_weight_perfect_matching() {
        std::vector<std::tuple<int, int, int>> edges;

        for (size_t i = 0; i < odd_vertices.size(); ++i) {
            for (size_t j = i + 1; j < odd_vertices.size(); ++j) {
                int u = odd_vertices[i];
                int v = odd_vertices[j];
                int weight = get_cost(u, v);
                edges.push_back({weight, u, v});
            }
        }

        std::sort(edges.begin(), edges.end());

        std::set<int> matched;
        matching_edges.clear();

        for (const auto& [weight, u, v] : edges) {
            if (matched.find(u) == matched.end() && matched.find(v) == matched.end()) {
                matching_edges.push_back(std::make_tuple(u, v, weight));
                matched.insert(u);
                matched.insert(v);
            }
        }
    }

    void find_euler_circuit() {
        mst_mwpm_graph = mst_edges;
        mst_mwpm_graph.insert(mst_mwpm_graph.end(), matching_edges.begin(), matching_edges.end());

        std::unordered_map<int, std::deque<int>> adj;

        for (const auto& [u, v, _] : mst_mwpm_graph) {
            adj[u].push_back(v);
            adj[v].push_back(u);
        }

        for (const auto& [node, neighbors] : adj) {
            if (neighbors.size() % 2 != 0) {
                std::cerr << "Vertex " << node << " has an odd degree, so no Euler circuit exists." << std::endl;
                return;
            }
        }

        euler_circuit.clear();
        std::deque<int> current_path = {adj.begin()->first};

        while (!current_path.empty()) {
            int u = current_path.back();

            if (!adj[u].empty()) {
                int v = adj[u].front();
                adj[u].pop_front();
                adj[v].erase(std::remove(adj[v].begin(), adj[v].end(), u), adj[v].end());
                current_path.push_back(v);
            } else {
                euler_circuit.push_back(current_path.back());
                current_path.pop_back();
            }
        }

        std::reverse(euler_circuit.begin(), euler_circuit.end());
    }

    std::pair<int, std::vector<int>> christofides() {
        prim_mst();
        find_odd_degree_vertices();
        minimum_weight_perfect_matching();
        find_euler_circuit();

        std::vector<int> unique_tour;
        std::set<int> seen;
        int cost = 0;

        for (int vertex : euler_circuit) {
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
    std::vector<std::pair<int, int>> euc_coordinates;
    std::vector<std::tuple<int, int, int>> mst_edges;
    std::vector<int> odd_vertices;
    std::vector<std::tuple<int, int, int>> matching_edges;
    std::vector<std::tuple<int, int, int>> mst_mwpm_graph;
    std::vector<int> euler_circuit;
    int n;
    int total_weight;
    bool explicit_weights;
};

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    ChristofidesSolver solver(filename);
    auto result = solver.christofides();
    std::cout << "Total cost: " << result.first << "\nPath: ";
    for (int city : result.second) {
        std::cout << city << " ";
    }
    std::cout << std::endl;
    return 0;
}