#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <random>
#include <set>
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <fstream>

using namespace std;
using namespace chrono;

class Graph {
private:
    int numVertices;
    bool directed;
    vector<vector<int>> adjacencyList;
    vector<pair<int, int>> edgeList;

public:
    Graph(int n, bool dir) : numVertices(n), directed(dir) {
        adjacencyList.resize(n);
    }

    void addEdge(int u, int v) {
        edgeList.emplace_back(u, v);
        adjacencyList[u].push_back(v);
        if (!directed) {
            adjacencyList[v].push_back(u);
        }
    }

    vector<vector<int>> getAdjacencyMatrix() const {
        vector<vector<int>> matrix(numVertices, vector<int>(numVertices, 0));
        for (const auto& edge : edgeList) {
            int u = edge.first;
            int v = edge.second;
            matrix[u][v] = 1;
            if (!directed) {
                matrix[v][u] = 1;
            }
        }
        return matrix;
    }

    vector<vector<int>> getIncidenceMatrix() const {
        int numEdges = edgeList.size();
        vector<vector<int>> matrix(numVertices, vector<int>(numEdges, 0));
        for (int i = 0; i < numEdges; ++i) {
            int u = edgeList[i].first;
            int v = edgeList[i].second;
            if (directed) {
                matrix[u][i] = 1;
                matrix[v][i] = -1;
            }
            else {
                matrix[u][i] = 1;
                matrix[v][i] = 1;
            }
        }
        return matrix;
    }

    vector<vector<int>> getAdjacencyList() const {
        return adjacencyList;
    }

    vector<pair<int, int>> getEdgeList() const {
        return edgeList;
    }

    int getNumVertices() const { return numVertices; }
    bool isDirected() const { return directed; }
};

class GraphGenerator {
public:
    struct Parameters {
        int min_nodes;
        int max_nodes;
        int min_edges;
        int max_edges;
        int max_edges_per_vertex;
        bool directed;
        int max_in_degree;
        int max_out_degree;
    };

    Graph generate(const Parameters& params) {
        random_device rd;
        mt19937 gen(rd());

        int num_nodes = uniform_int_distribution<>(params.min_nodes, params.max_nodes)(gen);

        int max_possible_edges;
        if (params.directed) {
            max_possible_edges = num_nodes * (num_nodes - 1);
        }
        else {
            max_possible_edges = num_nodes * (num_nodes - 1) / 2;
        }

        if (params.directed) {
            int max_out_total = num_nodes * params.max_out_degree;
            int max_in_total = num_nodes * params.max_in_degree;
            max_possible_edges = min(max_possible_edges, min(max_out_total, max_in_total));
        }
        else {
            int max_edges_total = (num_nodes * params.max_edges_per_vertex) / 2;
            max_possible_edges = min(max_possible_edges, max_edges_total);
        }

        int min_edges = max(params.min_edges, 0);
        int max_edges = min(params.max_edges, max_possible_edges);
        if (min_edges > max_edges) {
            throw invalid_argument("Invalid edge parameters");
        }
        int num_edges = uniform_int_distribution<>(min_edges, max_edges)(gen);

        Graph graph(num_nodes, params.directed);
        set<pair<int, int>> existing_edges;
        vector<int> out_degree(num_nodes, 0);
        vector<int> in_degree(num_nodes, 0);

        for (int e = 0; e < num_edges; ) {
            int u = uniform_int_distribution<>(0, num_nodes - 1)(gen);

            if (params.directed) {
                if (out_degree[u] >= params.max_out_degree) continue;
            }
            else {
                if (out_degree[u] >= params.max_edges_per_vertex) continue;
            }

            int v;
            do {
                v = uniform_int_distribution<>(0, num_nodes - 1)(gen);
            } while (v == u);

            if (!params.directed && u > v) {
                swap(u, v);
            }

            if (params.directed) {
                if (in_degree[v] >= params.max_in_degree) continue;
            }
            else {
                if (out_degree[v] >= params.max_edges_per_vertex) continue;
            }

            pair<int, int> edge(u, v);
            if (existing_edges.count(edge)) continue;

            graph.addEdge(u, v);
            existing_edges.insert(edge);
            out_degree[u]++;
            in_degree[v]++;
            if (params.directed) {
                // pass
            }
            else {
                out_degree[v]++;
                in_degree[u]++;
            }
            e++;
        }

        return graph;
    }
};

vector<int> BFS(const Graph& graph, int start, int end, double& duration) {
    auto start_time = high_resolution_clock::now();
    int n = graph.getNumVertices();
    vector<bool> visited(n, false);
    vector<int> parent(n, -1);
    queue<int> q;

    visited[start] = true;
    q.push(start);

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        if (u == end) break;

        auto adjList = graph.getAdjacencyList();
        for (int v : adjList[u]) {
            if (!visited[v]) {
                visited[v] = true;
                parent[v] = u;
                q.push(v);
            }
        }
    }

    vector<int> path;
    if (!visited[end]) {
        duration = duration_cast<microseconds>(
            high_resolution_clock::now() - start_time).count();
        return path;
    }

    int current = end;
    while (current != -1) {
        path.push_back(current);
        current = parent[current];
    }
    reverse(path.begin(), path.end());

    duration = duration_cast<microseconds>(
        high_resolution_clock::now() - start_time).count();
    return path;
}

vector<int> DFS(const Graph& graph, int start, int end, double& duration) {
    auto start_time = high_resolution_clock::now();
    int n = graph.getNumVertices();
    vector<bool> visited(n, false);
    vector<int> parent(n, -1);
    stack<int> stack;

    stack.push(start);
    visited[start] = true;
    bool found = false;

    while (!stack.empty()) {
        int u = stack.top();
        stack.pop();

        if (u == end) {
            found = true;
            break;
        }

        auto adjList = graph.getAdjacencyList();
        for (int v : adjList[u]) {
            if (!visited[v]) {
                visited[v] = true;
                parent[v] = u;
                stack.push(v);
            }
        }
    }

    vector<int> path;
    if (!found) {
        duration = duration_cast<microseconds>(
            high_resolution_clock::now() - start_time).count();
        return path;
    }

    int current = end;
    while (current != -1) {
        path.push_back(current);
        current = parent[current];
    }
    reverse(path.begin(), path.end());

    duration = duration_cast<microseconds>(
        high_resolution_clock::now() - start_time).count();
    return path;
}

int main() {
    GraphGenerator::Parameters params;
    params.min_nodes = 5;
    params.max_nodes = 5;
    params.min_edges = 4;
    params.max_edges = 6;
    params.max_edges_per_vertex = 3;
    params.directed = false;
    params.max_in_degree = 3;
    params.max_out_degree = 3;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> directed_dist(0, 1);
    params.directed = directed_dist(gen) == 1;

    GraphGenerator generator;

    // Открываем файл для записи результатов
    ofstream outputFile("results4.csv");
    outputFile << "Graph,Vertices,Edges,BFS Time (microseconds),DFS Time (microseconds),BFS Path,DFS Path\n";

    for (int i = 1; i < 21; ++i) {
        params.min_nodes = 1 + i * 5;
        params.max_nodes = 1 + i * 5;
        params.min_edges = 1 + i * 7;
        params.max_edges = 2 + i * 10;

        try {
            Graph graph = generator.generate(params);

            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<> dis(0, graph.getNumVertices() - 1);
            int A = dis(gen);
            int B = dis(gen);

            double bfs_time, dfs_time;
            auto bfs_path = BFS(graph, A, B, bfs_time);
            auto dfs_path = DFS(graph, A, B, dfs_time);

            cout << "Graph " << i << ":\n";
            cout << "Nodes: " << graph.getNumVertices() << "\n";
            cout << "Edges: " << graph.getEdgeList().size() << "\n";
            cout << "A: " << A << ", B: " << B << "\n";

            if (i <= 1) {
                cout << "Adjacency Matrix:\n";
                auto adjMatrix = graph.getAdjacencyMatrix();
                for (const auto& row : adjMatrix) {
                    for (int val : row) {
                        cout << val << " ";
                    }
                    cout << "\n";
                }

                cout << "Incidence Matrix:\n";
                auto inciMatrix = graph.getIncidenceMatrix();
                for (const auto& row : inciMatrix) {
                    for (int val : row) {
                        cout << val << " ";
                    }
                    cout << "\n";
                }

                cout << "Adjacency List:\n";
                auto adjList = graph.getAdjacencyList();
                for (size_t j = 0; j < adjList.size(); ++j) {
                    cout << j << ": ";
                    for (int neighbor : adjList[j]) {
                        cout << neighbor << " ";
                    }
                    cout << "\n";
                }

                cout << "Edge List:\n";
                for (const auto& edge : graph.getEdgeList()) {
                    cout << edge.first << " - " << edge.second << "\n";
                }
            }

            cout << "BFS Time: " << bfs_time << " microseconds\n";
            cout << "DFS Time: " << dfs_time << " microseconds\n";

            // Вывод путей
            cout << "BFS Path: ";
            if (bfs_path.empty()) {
                cout << "No path found";
            }
            else {
                for (int node : bfs_path) {
                    cout << node << " ";
                }
            }
            cout << "\n";

            cout << "DFS Path: ";
            if (dfs_path.empty()) {
                cout << "No path found";
            }
            else {
                for (int node : dfs_path) {
                    cout << node << " ";
                }
            }
            cout << "\n";

            cout << "-------------------------\n";

            // Записываем результаты в CSV-файл
            outputFile << i << "," << graph.getNumVertices() << "," << graph.getEdgeList().size() << ","
                << bfs_time << "," << dfs_time << ",";

            // Записываем BFS Path
            if (bfs_path.empty()) {
                outputFile << "No path";
            }
            else {
                for (int node : bfs_path) {
                    outputFile << node << " ";
                }
            }
            outputFile << ",";

            // Записываем DFS Path
            if (dfs_path.empty()) {
                outputFile << "No path";
            }
            else {
                for (int node : dfs_path) {
                    outputFile << node << " ";
                }
            }
            outputFile << "\n";
        }
        catch (const exception& e) {
            cerr << "Error generating graph " << i << ": " << e.what() << "\n";
        }
    }

    outputFile.close();
    cout << "Results saved to results.csv\n";

    return 0;
}
