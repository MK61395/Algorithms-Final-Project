// I21-0851 Muhammad Kashif
// I21-0556 Muhammad Usman Nazeer
// I21-0464 Bilal Saleem
// I21-2460 Muhammad Huzaifa
// P
// CS-2009 Design and Analysis of Algorithm
// Final Project

#include <iostream>
#include <string>
#include <queue>
#include <random>
#include <stack>
#include <limits>
#include <chrono>
#include <fstream>
#include <sstream>
using namespace std;


const int MAX_NODES = 50;
const int MAX_EVENTS = 3;


class Node {
public:
    string name; // Name of the node
    string role; // Role of the node
    string department; // Department of the node
    int academicImpact; // Academic impact level of the node
    int leadershipRoles; // Number of leadership roles held by the node
    int involvement; // Level of involvement in events

    // Array to track participation in events, true if participated, false otherwise
    bool isEventParticipant[MAX_EVENTS];

    // Constructor to initialize node properties
    Node(const string& n, const string& r, const string& d, int impact, int leadership, int events)
        : name(n), role(r), department(d), academicImpact(impact), leadershipRoles(leadership), involvement(events) {
        // Initialize the participation in events array to false for all events
        for (int i = 0; i < MAX_EVENTS; i++) {
            isEventParticipant[i] = false;
        }
    }
};


class Graph {
private:
    Node* nodes[MAX_NODES]; // Array to hold pointers to nodes
    double** adjacencyMatrix; // Matrix representing connections between nodes based on edge weights
    int** relationshipMatrix; // Matrix representing relationships between nodes
    int event_count; // Count of events in the graph
    string* event_names; // Array to store names of events
    int nodeCount; // Count of nodes in the graph

public:
    // Constructor to initialize the graph
    Graph() : nodeCount(0), event_count(0) {
        adjacencyMatrix = new double* [MAX_NODES]; // Allocate memory for adjacency matrix
        relationshipMatrix = new int* [MAX_NODES]; // Allocate memory for relationship matrix
        this->event_names = new string[MAX_EVENTS]; // Allocate memory for event names array

        int i = 0;
        while (i < MAX_NODES) {
            adjacencyMatrix[i] = new double[MAX_NODES]; // Initialize adjacency matrix
            relationshipMatrix[i] = new int[MAX_NODES]; // Initialize relationship matrix

            int j = 0;
            while (j < MAX_NODES) {
                adjacencyMatrix[i][j] = 0.0; // Set default value for adjacency matrix
                relationshipMatrix[i][j] = 0; // Set default value for relationship matrix
                j++;
            }
            i++;
        }
    }

    // Destructor to free dynamically allocated memory
    ~Graph() {
        for (int i = 0; i < MAX_NODES; ++i) {
            delete[] adjacencyMatrix[i]; // Free memory for each row of adjacency matrix
            delete[] relationshipMatrix[i]; // Free memory for each row of relationship matrix
        }
        delete[] adjacencyMatrix; // Free memory allocated for adjacency matrix
        delete[] relationshipMatrix; // Free memory allocated for relationship matrix
    }



    // This function adds a new node to our graph.
    void addNode(Node* newNode) {
        // We put the new node into our list of nodes and increase the count.
        nodes[nodeCount++] = newNode;
    }

    // This function adds a relationship (an edge) between two nodes in the graph.
    // It sets a random weight and a type of relationship between them.
    void addEdge(Node* src, Node* dest, double minWeight, double maxWeight, int relationshipType) {
        // We find the positions (indices) of the source and destination nodes in our node list.
        int source = getNodeIndex(src);
        int destination = getNodeIndex(dest);

        // If both nodes exist in our list...
        if (source != -1 && destination != -1) {
            // ...we assign a random weight between a minimum and maximum value to the edge.
            double weight = ((double)rand() / RAND_MAX) * (maxWeight - minWeight) + minWeight;
            // We update the adjacency matrix with this weight for both directions (source to destination and vice versa).
            adjacencyMatrix[source][destination] = weight;
            adjacencyMatrix[destination][source] = weight;

            // We set the type of relationship in the relationship matrix for both directions as well.
            relationshipMatrix[source][destination] = relationshipType;
            relationshipMatrix[destination][source] = relationshipType;
        }
    }

    // This function prints out the connections between nodes in the graph.
    void printGraph() {
        // For each node in our list...
        for (int i = 0; i < nodeCount; ++i) {
            // ...we display its name and connections.
            cout << nodes[i]->name << " is connected to: " << endl;
            // Check each other node to see if there's a connection.
            for (int j = 0; j < nodeCount; ++j) {
                // If there's a non-zero value in the adjacency matrix, there's a connection.
                if (adjacencyMatrix[i][j] != 0) {
                    // Print the connected node's name, weight of connection, and type of relationship.
                    cout << nodes[j]->name << " (Weight: " << adjacencyMatrix[i][j] << ", Relationship: ";
                    // Translate the relationship type to a readable text.
                    switch (relationshipMatrix[i][j]) {
                    case 1:
                        cout << "Academic Collaboration";
                        break;
                    case 2:
                        cout << "Mentorship";
                        break;
                    case 3:
                        cout << "Extracurricular Involvement";
                        break;
                    default:
                        cout << "Unknown Relationship";
                        break;
                    }
                    cout << ")" << endl;
                }
            }
            cout << endl;
        }
    }


    // This function prints the adjacency matrix representing connections between nodes in the graph.
    void printAdjacencyMatrix() {
        cout << "Adjacency Matrix:" << endl;
        int i = 0;
        while (i < nodeCount) {
            int j = 0;
            while (j < nodeCount) {
                // Output the values of the adjacency matrix.
                cout << adjacencyMatrix[i][j] << "  ";
                j++;
            }
            cout << endl;
            i++;
        }
    }

    // This function performs a breadth-first search traversal starting from a given node.
    void bfs(Node* startNode) {
        bool visited[MAX_NODES] = { false }; // Track visited nodes.
        queue<Node*> q; // Queue for BFS traversal.

        int startIndex = getNodeIndex(startNode);
        if (startIndex == -1) {
            cout << "Start node not found in the graph." << endl;
            return;
        }

        visited[startIndex] = true;
        q.push(startNode);

        cout << "BFS Traversal from Node " << startNode->name << ":" << endl;

        // BFS algorithm implementation.
        while (!q.empty()) {
            Node* currentNode = q.front();
            q.pop();
            cout << currentNode->name << " ";

            int currentIndex = getNodeIndex(currentNode);
            int i = 0;
            while (i < nodeCount) {
                if (adjacencyMatrix[currentIndex][i] != 0 && !visited[i]) {
                    visited[i] = true;
                    q.push(nodes[i]);
                }
                i++;
            }
        }
        cout << endl;
    }

    // This function performs a depth-first search traversal starting from a given node.
    void dfs(Node* startNode) {
        bool visited[MAX_NODES] = { false }; // Track visited nodes.
        stack<Node*> s; // Stack for DFS traversal.

        int startIndex = getNodeIndex(startNode);
        if (startIndex == -1) {
            cout << "Start node not found in the graph." << endl;
            return;
        }

        visited[startIndex] = true;
        s.push(startNode);

        cout << "DFS Traversal from Node " << startNode->name << ":" << endl;

        // DFS algorithm implementation.
        while (!s.empty()) {
            Node* currentNode = s.top();
            s.pop();
            if (!visited[getNodeIndex(currentNode)]) {
                cout << currentNode->name << " ";
                visited[getNodeIndex(currentNode)] = true;
            }

            int i = 0;
            while (i < nodeCount) {
                if (adjacencyMatrix[getNodeIndex(currentNode)][i] != 0 && !visited[i]) {
                    s.push(nodes[i]);
                }
                i++;
            }
        }
        cout << endl;
    }

    // This function performs topological sorting of the nodes in the graph.
    void topologicalSort() {
        stack<Node*> stack; // Stack for topological sorting.
        bool visited[MAX_NODES] = { false }; // Track visited nodes.

        // Visit each node and perform topological sort.
        int i = 0;
        while (i < nodeCount) {
            if (!visited[i]) {
                topologicalSortUtil(i, visited, stack);
            }
            i++;
        }

        // Output the topologically sorted nodes.
        cout << "Topological Sorting: ";
        while (!stack.empty()) {
            cout << stack.top()->name << " ";
            stack.pop();
        }
        cout << endl;
    }


    // This function finds the Minimum Spanning Tree (MST) using Prim's algorithm.
    void primMST() {
        // Initialization
        bool inMST[MAX_NODES] = { false }; // Tracks nodes in the MST.
        double key[MAX_NODES]; // Holds minimum weights for nodes.
        int parent[MAX_NODES]; // Stores the parent nodes.

        // Initializing keys and parent nodes.
        int i = 0;
        while (i < nodeCount) {
            key[i] = numeric_limits<double>::max(); // Initialize keys with maximum values.
            parent[i] = -1; // Initialize parent nodes as -1.
            i++;
        }

        key[0] = 0; // Starting with the first node.
        parent[0] = -1; // No parent for the starting node.

        // Finding MST
        int count = 0;
        double totalCost = 0;
        while (count < nodeCount - 1) {
            int u = minKey(key, inMST); // Find the minimum key node not in MST.

            if (u == -1) {
                cout << "Error: Graph is not connected." << endl;
                return;
            }

            // Include the selected node in the MST.
            inMST[u] = true;
            int v = 0;
            while (v < nodeCount) {
                // Update keys and parent nodes for adjacent vertices.
                if (adjacencyMatrix[u][v] != 0 && !inMST[v] && adjacencyMatrix[u][v] < key[v]) {
                    parent[v] = u;
                    key[v] = adjacencyMatrix[u][v];
                }
                v++;
            }
            count++;
        }

        // Printing the MST edges and their total cost.
        cout << "Minimum Spanning Tree (using Prim's algorithm):" << endl;
        i = 1;
        while (i < nodeCount) {
            if (parent[i] != -1) {
                cout << nodes[parent[i]]->name << " - " << nodes[i]->name << " (Weight: " << adjacencyMatrix[i][parent[i]] << ")" << endl;
                totalCost += adjacencyMatrix[i][parent[i]];
            }
            i++;
        }
        cout << "Total Cost of Minimum Spanning Tree: " << totalCost << endl;
    }

    // This function calculates the degree centrality for each node in the graph.
    void calculateDegreeCentrality() {
        cout << "Degree Centrality for each node:" << endl;
        int i = 0;
        while (i < nodeCount) {
            // Count connections to calculate degree centrality.
            int j = 0;
            int connections = 0;
            while (j < nodeCount) {
                if (adjacencyMatrix[i][j] != 0) {
                    connections++;
                }
                j++;
            }
            // Calculate and display degree centrality for each node.
            double degreeCentrality = connections / static_cast<double>(nodeCount - 1);
            cout << nodes[i]->name << ": " << degreeCentrality << endl;
            i++;
        }
    }

    // This function finds and prints nodes with a specified degree centrality.
    void printNodesWithSameDegreeCentrality(const double targetCentrality, const Node* const* nodes, const int nodeCount) {
        cout << "Nodes with the same degree centrality (" << targetCentrality << "):" << endl;
        bool found = false;
        int i = 0;
        while (i < nodeCount) {
            // Calculate degree centrality for each node.
            int j = 0;
            int connections = 0;
            while (j < nodeCount) {
                if (adjacencyMatrix[i][j] != 0) {
                    connections++;
                }
                j++;
            }
            double degreeCentrality = connections / static_cast<double>(nodeCount - 1);

            // Check and print nodes with the specified degree centrality.
            if (degreeCentrality == targetCentrality) {
                found = true;
                cout << nodes[i]->name << endl;
            }
            i++;
        }

        if (!found) {
            cout << "No other nodes found with the same degree centrality." << endl;
        }
    }

    // This function finds nodes with the highest degree centrality in the graph.
    void findNodesWithHighestDegreeCentrality() {
        double maxDegreeCentrality = 0.0;
        Node* nodeWithMaxDegreeCentrality = nullptr;

        cout << "Nodes with the highest degree centrality:" << endl;
        int i = 0;
        while (i < nodeCount) {
            // Calculate degree centrality for each node.
            int connections = 0;
            int j = 0;
            while (j < nodeCount) {
                if (adjacencyMatrix[i][j] != 0) {
                    connections++;
                }
                j++;
            }
            double degreeCentrality = connections / static_cast<double>(nodeCount - 1);

            // Find nodes with the highest degree centrality.
            if (degreeCentrality > maxDegreeCentrality) {
                maxDegreeCentrality = degreeCentrality;
                nodeWithMaxDegreeCentrality = nodes[i];
            }
            else if (degreeCentrality == maxDegreeCentrality) {
                // Print nodes with the same degree centrality if encountered.
                if (nodeWithMaxDegreeCentrality != nullptr) {
                    cout << nodeWithMaxDegreeCentrality->name << endl;
                }
                printNodesWithSameDegreeCentrality(maxDegreeCentrality, nodes, nodeCount);
                return;
            }
            i++;
        }

        // Print node(s) with the highest degree centrality.
        if (nodeWithMaxDegreeCentrality != nullptr) {
            cout << nodeWithMaxDegreeCentrality->name << endl;
        }
    }


    // This function finds academic collaborations within the same department between faculty and students.
    void findAcademicCollaborations() {
        cout << "Academic Collaborations within the same department:" << endl;
        int i = 0;
        while (i < nodeCount) {
            int j = i + 1;
            while (j < nodeCount) {
                if (nodes[i]->role == "Faculty" && nodes[j]->role == "Student" && nodes[i]->department == nodes[j]->department && adjacencyMatrix[i][j] != 0) {
                    // Display academic collaborations between faculty and students in the same department.
                    cout << nodes[i]->name << " (" << nodes[i]->role << ", " << nodes[i]->department << ") - " << nodes[j]->name << " (" << nodes[j]->role << ", " << nodes[j]->department << ")" << endl;
                }
                j++;
            }
            i++;
        }
    }

    // This function finds interdisciplinary collaborations across different departments.
    void findInterdisciplinaryCollaborations() {
        cout << "Interdisciplinary Collaborations across different departments:" << endl;
        int i = 0;
        while (i < nodeCount) {
            int j = i + 1;
            while (j < nodeCount) {
                if (nodes[i]->department != nodes[j]->department && adjacencyMatrix[i][j] != 0) {
                    // Display collaborations between nodes from different departments.
                    cout << nodes[i]->name << " (" << nodes[i]->role << ", " << nodes[i]->department << ") - " << nodes[j]->name << " (" << nodes[j]->role << ", " << nodes[j]->department << ")" << endl;
                }
                j++;
            }
            i++;
        }
    }

    // This function finds influential individuals based on academic impact, leadership roles, and event involvement.
    void findInfluentialIndividuals() {
        Node* influentialNodes[MAX_NODES]; // Array to store influential nodes.
        int influentialCount = 0;

        int i = 0;
        while (i < nodeCount) {
            if (nodes[i]->academicImpact > 80 && nodes[i]->leadershipRoles > 5 && nodes[i]->involvement > 20) {
                // Check if a node is influential based on specified criteria.
                influentialNodes[influentialCount++] = nodes[i];
            }
            i++;
        }

        if (influentialCount == 0) {
            cout << "No influential individuals found based on the criteria." << endl;
        }
        else {
            // Display influential individuals based on the specified criteria.
            cout << "Influential Individuals based on academic impact, leadership roles, and involvement in events:" << endl;
            i = 0;
            while (i < influentialCount) {
                cout << "Name: " << influentialNodes[i]->name << ", Role: " << influentialNodes[i]->role << ", Department: " << influentialNodes[i]->department << endl;
                i++; // Moved the increment inside the loop
            }
        }
    }

    // This function finds the Minimum Spanning Tree (MST) using Kruskal's algorithm.
    void kruskalMST() {
        // Initialization
        double edgeWeights[MAX_NODES * MAX_NODES];
        int edges[MAX_NODES * MAX_NODES][2];
        int edgeCount = 0;

        // Store edges and their weights
        int i = 0;
        while (i < nodeCount) {
            int j = i + 1;
            while (j < nodeCount) {
                if (adjacencyMatrix[i][j] != 0) {
                    edgeWeights[edgeCount] = adjacencyMatrix[i][j];
                    edges[edgeCount][0] = i;
                    edges[edgeCount][1] = j;
                    edgeCount++;
                }
                j++;
            }
            i++;
        }

        // Sort edges based on their weights
        i = 0;
        while (i < edgeCount - 1) {
            int j = 0;
            while (j < edgeCount - i - 1) {
                if (edgeWeights[j] > edgeWeights[j + 1]) {
                    swap(edgeWeights[j], edgeWeights[j + 1]);
                    swap(edges[j][0], edges[j + 1][0]);
                    swap(edges[j][1], edges[j + 1][1]);
                }
                j++;
            }
            i++;
        }

        // Initialize sets for union-find
        int parent[MAX_NODES];
        i = 0;
        while (i < nodeCount) {
            parent[i] = i;
            i++;
        }

        // Find MST using Kruskal's algorithm
        cout << "Minimum Spanning Tree (using Kruskal's algorithm):" << endl;
        i = 0;
        while (i < edgeCount) {
            int u = edges[i][0];
            int v = edges[i][1];
            int setU = findSet(parent, u);
            int setV = findSet(parent, v);

            if (setU != setV) {
                cout << nodes[u]->name << " - " << nodes[v]->name << " (Weight: " << edgeWeights[i] << ")" << endl;
                unionSets(parent, setU, setV);
            }
            i++;
        }
    }


    // Finds the shortest path distances from a source node using Dijkstra's algorithm.
    void dijkstra(Node* sourceNode) {
        // Initialization
        int sourceIndex = getNodeIndex(sourceNode);

        if (sourceIndex == -1) {
            cout << "Source node not found in the graph." << endl;
            return;
        }

        double distance[MAX_NODES];
        for (int i = 0; i < MAX_NODES; ++i) {
            distance[i] = std::numeric_limits<double>::max(); // Initialize distances to maximum value.
        }
        distance[sourceIndex] = 0.0; // Set the distance to the source node as 0.

        queue<int> pq;
        pq.push(sourceIndex); // Push the source node index into the priority queue.

        // Dijkstra's algorithm to compute shortest paths
        while (!pq.empty()) {
            int u = pq.front();
            pq.pop();
            int v = 0;
            while (v < nodeCount) {
                if (adjacencyMatrix[u][v] != 0) {
                    double weight = adjacencyMatrix[u][v];

                    // Relax edges and update distances if a shorter path is found.
                    if (distance[u] + weight < distance[v]) {
                        distance[v] = distance[u] + weight;
                        pq.push(v);
                    }
                }
                v++;
            }
        }

        // Display shortest distances from the source node to all other nodes.
        cout << "Shortest Distances from Node " << sourceNode->name << " using Dijkstra's algorithm:" << endl;
        for (int i = 0; i < nodeCount; ++i) {
            cout << sourceNode->name << " to " << nodes[i]->name << ": " << distance[i] << endl;
        }
    }

    // Adds an event to the event list.
    void addEvent(const string& e) {
        if (event_count < MAX_EVENTS) {
            event_names[event_count++] = e;
        }
        else {
            cout << "Number of events has reached its maximum limit." << endl;
        }
    }

    // Marks a participant's attendance for a specific event.
    void Attendance(Node* participant, const string& eventName) {
        int eventIndex = -1;
        // Find the index of the event in the event list.
        int i = 0;
        while ( i < event_count) {
            if (event_names[i] == eventName) {
                eventIndex = i;
                break;
            }
            i++;
        }
        if (eventIndex != -1) {
            participant->isEventParticipant[eventIndex] = true;
        }
        else {
            cout << "Event not found." << endl;
        }
    }

    // Identifies and displays popular events based on attendance.
    void popular_events() {
        // Calculate attendance per event.
        int attendancePerEvent[MAX_EVENTS] = { 0 };
        int i = 0;
        while (i < event_count) {
            int j = 0;
            while (j < nodeCount) {
                if (nodes[j]->isEventParticipant[i]) {
                    attendancePerEvent[i]++;
                }
                j++;
            }
            i++;
        }

        // Display events with higher attendance than a threshold.
        cout << "Popular Events based on Attendance:" << endl;
        int k = 0;
        while (k < event_count) {
            if (attendancePerEvent[k] > (nodeCount / 50)) {
                cout << event_names[k] << " - Attendance: " << attendancePerEvent[k] << " participants" << endl;
            }
            k++;
        }
    }

    // Identifies and displays frequent participants in events.
    void Frequent_participants() {
        int Count[MAX_NODES] = { 0 }; // Array to count event participation.

        // Count event participation for each node.
        int i = 0;
        while (i < nodeCount) {
            int j = 0;
            while (j < event_count) {
                if (nodes[i]->isEventParticipant[j]) {
                    Count[i]++;
                }
                j++;
            }
            i++;
        }

        // Display participants who participated in more than half of the events.
        cout << "Frequent Participants in Events:" << endl;
        int k = 0;
        while (k < nodeCount) {
            if (Count[k] > (event_count / 2)) {
                cout << nodes[k]->name << " - Participated in " << Count[k] << " events" << endl;
            }
            k++;
        }
    }

    // Removes a node from the graph and updates connections accordingly.
    void removeNode(Node* nodeToRemove) {
        // Find the index of the node to remove.
        int index = getNodeIndex(nodeToRemove);
        if (index == -1) {
            cout << "Node does not exist." << endl;
            return;
        }

        // Shift elements to remove the node from the nodes array.
        int i = index;
        while (i < nodeCount - 1) {
            nodes[i] = nodes[i + 1];
            int j = 0;
            while ( j < nodeCount) {
                adjacencyMatrix[i][j] = adjacencyMatrix[i + 1][j];
                adjacencyMatrix[j][i] = adjacencyMatrix[j][i + 1];
                relationshipMatrix[i][j] = relationshipMatrix[i + 1][j];
                relationshipMatrix[j][i] = relationshipMatrix[j][i + 1];
                j++;
            }
            ++i;
        }

        nodeCount--;

        delete nodes[nodeCount];
        nodes[nodeCount] = nullptr;

        // Clear connections related to the removed node.
        for (int i = 0; i < nodeCount; ++i) {
            adjacencyMatrix[nodeCount][i] = 0;
            adjacencyMatrix[i][nodeCount] = 0;
            relationshipMatrix[nodeCount][i] = 0;
            relationshipMatrix[i][nodeCount] = 0;
        }
        cout << "Node has been removed." << endl;
    }

    // Adds a new node to the graph.
    void add_new_Node(Node* newNode) {
        if (nodeCount < MAX_NODES) {
            nodes[nodeCount++] = newNode;
            cout << "Node has been added." << endl;
        }
        else {
            cout << "Maximum node limit reached." << endl;
        }
    }

    // Finds the shortest path distances from a source node using Bellman-Ford algorithm.
    void bellmanFord(Node* sourceNode) {
        // Initialization
        int sourceIndex = getNodeIndex(sourceNode);

        if (sourceIndex == -1) {
            cout << "Source node not found in the graph." << endl;
            return;
        }

        double distance[MAX_NODES];
        int i = 0;
        while (i < MAX_NODES) {
            distance[i] = numeric_limits<double>::max();
            i++;
        }
        distance[sourceIndex] = 0.0;

        // Bellman-Ford algorithm to compute shortest paths
        i = 0;
        while (i < nodeCount - 1) {
            int u = 0;
            while (u < nodeCount) {
                int v = 0;
                while (v < nodeCount) {
                    if (adjacencyMatrix[u][v] != 0) {
                        double weight = adjacencyMatrix[u][v];

                        // Relax edges and update distances if a shorter path is found.
                        if (distance[u] + weight < distance[v]) {
                            distance[v] = distance[u] + weight;
                        }
                    }
                    v++;
                }
                u++;
            }
            i++;
        }

        // Check for negative weight cycles
        int u = 0;
        while (u < nodeCount) {
            int v = 0;
            while (v < nodeCount) {
                if (adjacencyMatrix[u][v] != 0) {
                    double weight = adjacencyMatrix[u][v];

                    if (distance[u] + weight < distance[v]) {
                        cout << "Graph contains a negative weight cycle. Bellman-Ford cannot find shortest paths." << endl;
                        return;
                    }
                }
                v++;
            }
            u++;
        }

        // Display shortest distances from the source node to all other nodes.
        cout << "Shortest Distances from Node " << sourceNode->name << " using Bellman-Ford algorithm:" << endl;
        i = 0;
        while (i < nodeCount) {
            cout << sourceNode->name << " to " << nodes[i]->name << ": " << distance[i] << endl;
            i++;
        }
    }

    void detectCommunities() {
        bool visited[MAX_NODES] = { false }; // Track visited nodes
        cout << "Communities within the university:" << endl;
        int i = 0;
        while( i < nodeCount) {
            if (!visited[i]) {
                cout << "Community " << i + 1 << ": ";
                detectCommunityDFS(i, visited);
                cout << endl;
            }
            i++;
        }
    }
private:
    // Function to get the index of a node in the 'nodes' array
    int getNodeIndex(Node* node) {
        int i = 0;
        while ( i < nodeCount) {
            if (nodes[i] == node) {
                return i;
            }
            i++;
        }
        return -1; // Return -1 if the node is not found
    }

    // Utility function for topological sorting using DFS
    void topologicalSortUtil(int v, bool visited[], stack<Node*>& stack) {
        visited[v] = true;

        // Recursively traverse the adjacent nodes
        int i = 0;
        while ( i < nodeCount) {
            if (adjacencyMatrix[v][i] != 0 && !visited[i]) {
                topologicalSortUtil(i, visited, stack);
            }
            i++;
        }

        // Push the current node onto the stack
        stack.push(nodes[v]);
    }

    // Function to find the node with the minimum key value from the set of vertices not yet included in MST
    int minKey(double key[], bool inMST[]) {
        double min = numeric_limits<double>::max();
        int min_index = -1;

        // Traverse all the vertices to find the minimum key
        int v = 0;
        while ( v < nodeCount) {
            if (!inMST[v] && key[v] < min) {
                min = key[v];
                min_index = v;
            }
            v++;
        }

        return min_index; // Return the index of the minimum key vertex
    }

    // Function to find the subset of an element 'i' in the disjoint set
    int findSet(int parent[], int i) {
        // Path compression technique to find the representative element of the set
        if (parent[i] != i) {
            parent[i] = findSet(parent, parent[i]);
        }
        return parent[i]; // Return the representative (root) of the set
    }

    // Function to perform union of two subsets x and y
    void unionSets(int parent[], int x, int y) {
        // Find the root nodes (representatives) of x and y
        int rootX = findSet(parent, x);
        int rootY = findSet(parent, y);

        // Union by rank: attach the smaller depth tree under the root of the deeper tree
        parent[rootX] = rootY;
    }

    void detectCommunityDFS(int nodeIndex, bool visited[]) {
        visited[nodeIndex] = true;
        cout << nodes[nodeIndex]->name << " ";

        // Traverse adjacent nodes for the current node
        int j = 0;
        while(j < nodeCount) {
            if (adjacencyMatrix[nodeIndex][j] != 0 && !visited[j]) {
                detectCommunityDFS(j, visited);
            }
            j++;
        }
    }
};

// Appends operation and time complexity data to a CSV file
void appendTimeComplexityToCSV(const string& fileName, const string& operation, double timeTaken) {
    ofstream file;
    file.open(fileName, ios_base::app); // Open the file in append mode

    // Check if the file is successfully opened
    if (file.is_open()) {
        // Write operation and time taken as a new line in CSV format
        file << operation << "," << timeTaken << "\n";
        file.close(); // Close the file after writing
    }
    else {
        cout << "Unable to open file to append data." << endl; // Display error if file opening fails
    }
}


int main()
{
    // Creating a graph object
    Graph universityGraph;

    // Array to store node pointers and counting the nodes
    Node* nodes[MAX_NODES];
    int nodeCount = 0;

    // Array to store node pointers and counting the nodes
    ifstream file("data.csv");
    if (!file.is_open()) {
        cout << "Error opening the file." << endl;
        return 1;
    }

    // Read node information from the file
    // Extracting data for nodes and populating the graph
    // Skipping the header line
    string line;
    getline(file, line); // Skip header

    while (getline(file, line) && nodeCount < MAX_NODES) {
        stringstream ss(line);// Parse the line and extract node attributes
        string nodeName, nodeRole, nodeDepartment;
        int academicImpact, leadershipRoles, involvement;

        // Create a new node and add it to the graph
        getline(ss, nodeName, ',');
        getline(ss, nodeRole, ',');
        getline(ss, nodeDepartment, ',');
        ss >> academicImpact;
        ss.ignore();
        ss >> leadershipRoles;
        ss.ignore();
        ss >> involvement;

        Node* newNode = new Node(nodeName, nodeRole, nodeDepartment, academicImpact, leadershipRoles, involvement);  // Storing node pointers in the nodes array
        universityGraph.addNode(newNode);
        nodes[nodeCount++] = newNode;// Incrementing the node count
    }

    file.close();

    // Set minimum and maximum edge weights
    // Populate graph edges based on certain conditions
    // Randomly generate relationship types for edges
    double minWeight = 0.5;
    double maxWeight = 1.5;

    for (int i = 0; i < nodeCount; ++i) {
        for (int j = i + 1; j < nodeCount; ++j) {
            if (i != j) {
                int relationshipType = rand() % 3 + 1;
                if (nodes[i]->role == "Faculty" && nodes[j]->role == "Student" && nodes[i]->department == nodes[j]->department) {
                    universityGraph.addEdge(nodes[i], nodes[j], minWeight, maxWeight * 1.5, relationshipType);
                }
                else if (nodes[i]->department != nodes[j]->department) {
                    universityGraph.addEdge(nodes[i], nodes[j], minWeight, maxWeight, relationshipType);
                }
            }
        }
    }
    
    universityGraph.addEvent("Event 1");
    universityGraph.addEvent("Event 2");
    universityGraph.addEvent("Event 3");


    int choice;
    // Menu-driven interface for various graph operations
    do {
        // Display menu options for different graph operations

        cout << "University Graph Operations Menu:" << endl;
        cout << "1. Print Graph" << endl;
        cout << "2. Print Adjacency Matrix" << endl;
        cout << "3. Breadth-First Search (BFS)" << endl;
        cout << "4. Depth-First Search (DFS)" << endl;
        cout << "5. Topological Sort" << endl;
        cout << "6. Minimum Spanning Tree (Prim's algorithm)" << endl;
        cout << "7. Calculate Degree Centrality" << endl;
        cout << "7. Find Nodes with Highest Degree Centrality" << endl;
        cout << "8. Find Academic Collaborations" << endl;
        cout << "9. Find Interdisciplinary Collaborations" << endl;
        cout << "10. Find Influential Individuals" << endl;
        cout << "11.Minimum Spanning Tree (Kruskal's algorithm)" << endl;
        cout << "12.Shortest Path  (Dijkstra's algorithm)" << endl;
        cout << "13.Event Attendance Analysis" << endl;
        cout << "14.To Add New Node(Individuals)" << endl;
        cout << "15.To Remove Node(Individuals)" << endl;
        cout << "16.Bellman's Shortest Path" << endl;
        cout << "17.Detect Communities" << endl;
        string n, r, d;

        //Take user input for the desired operation
        cout << "0. Exit" << endl;
        int option = 0;
        cout << "Enter your choice: ";
        cin >> choice;

        switch (choice)
        {
        // Perform specific operations based on the user's choice
        // Measure time taken for certain operations using chrono library
        // Display the time taken for each operation
        // Append the time complexities to a CSV file using 'appendTimeComplexityToCSV'
        // Switch-case block to handle different user choices for graph operations
        // Each case corresponds to a specific operation using methods from the 'universityGraph' object
        // Measure and record time taken for operations and write to 'time_complexities.csv'
        // Continue until the user chooses to exit
        case 1:
            universityGraph.printGraph();
            break;
        case 2:
            universityGraph.printAdjacencyMatrix();
            break;
        case 3: {
            auto start = chrono::high_resolution_clock::now();
            int index = 0;
            cout << "Enter starting Index" << endl;
            cin >> index;
            universityGraph.bfs(nodes[0]);

            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> duration = end - start;

            cout << "Time taken for BFS: " << duration.count() << " seconds" << endl;
            appendTimeComplexityToCSV("time_complexities.csv", "BFS", duration.count());
            break;
        }
        case 4: {
            auto start = chrono::high_resolution_clock::now();
            int index = 0;
            cout << "Enter starting Index" << endl;
            cin >> index;
            universityGraph.dfs(nodes[index]);

            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> duration = end - start;

            cout << "Time taken for DFS: " << duration.count() << " seconds" << endl;
            appendTimeComplexityToCSV("time_complexities.csv", "DFS", duration.count());
            break;
        }

        case 5: {
            auto start = chrono::high_resolution_clock::now();
            universityGraph.topologicalSort();

            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> duration = end - start;

            cout << "Time taken for Topological Sort: " << duration.count() << " seconds" << endl;
            appendTimeComplexityToCSV("time_complexities.csv", "Topological Sort", duration.count());
            break;
        }

        case 6: {
            auto start = chrono::high_resolution_clock::now();
            universityGraph.primMST();

            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> duration = end - start;

            cout << "Time taken for Prim's MST: " << duration.count() << " seconds" << endl;
            appendTimeComplexityToCSV("time_complexities.csv", "primMST", duration.count());
            break;
        }
        case 7: {
            auto startDegreeCentrality = chrono::high_resolution_clock::now();
            universityGraph.calculateDegreeCentrality();
            universityGraph.findNodesWithHighestDegreeCentrality();

            auto endDegreeCentrality = chrono::high_resolution_clock::now();
            chrono::duration<double> durationDegreeCentrality = endDegreeCentrality - startDegreeCentrality;

            cout << "Time taken for Degree Centrality and finding Nodes with Highest Degree Centrality: "
                << durationDegreeCentrality.count() << " seconds" << endl;
            appendTimeComplexityToCSV("time_complexities.csv", "DegreeCentrality", durationDegreeCentrality.count());
            break;
        }

        case 8: {
            auto startAcademicCollaborations = chrono::high_resolution_clock::now();
            universityGraph.findAcademicCollaborations();

            auto endAcademicCollaborations = chrono::high_resolution_clock::now();
            chrono::duration<double> durationAcademicCollaborations = endAcademicCollaborations - startAcademicCollaborations;

            cout << "Time taken for Finding Academic Collaborations: "
                << durationAcademicCollaborations.count() << " seconds" << endl;
            appendTimeComplexityToCSV("time_complexities.csv", "AcademicCollaborations", durationAcademicCollaborations.count());
            break;
        }

        case 9: {
            auto startInterdisciplinaryCollaborations = chrono::high_resolution_clock::now();
            universityGraph.findInterdisciplinaryCollaborations();

            auto endInterdisciplinaryCollaborations = chrono::high_resolution_clock::now();
            chrono::duration<double> durationInterdisciplinaryCollaborations = endInterdisciplinaryCollaborations - startInterdisciplinaryCollaborations;

            cout << "Time taken for Finding Interdisciplinary Collaborations: "
                << durationInterdisciplinaryCollaborations.count() << " seconds" << endl;
            appendTimeComplexityToCSV("time_complexities.csv", "InterdisciplinaryCollaborations", durationInterdisciplinaryCollaborations.count());
            break;
        }

        case 10: {
            auto startInfluentialIndividuals = chrono::high_resolution_clock::now();
            universityGraph.findInfluentialIndividuals();

            auto endInfluentialIndividuals = chrono::high_resolution_clock::now();
            chrono::duration<double> durationInfluentialIndividuals = endInfluentialIndividuals - startInfluentialIndividuals;

            cout << "Time taken for Finding Influential Individuals: "
                << durationInfluentialIndividuals.count() << " seconds" << endl;
            appendTimeComplexityToCSV("time_complexities.csv", "InfluentialIndividual", durationInfluentialIndividuals.count());
            break;
        }

        case 11: {
            auto startKruskals = chrono::high_resolution_clock::now();
            universityGraph.kruskalMST();

            auto endKruskals = chrono::high_resolution_clock::now();
            chrono::duration<double> durationInfluentialIndividuals = endKruskals - startKruskals;

            cout << "Time taken for Finding Influential Individuals: "
                << durationInfluentialIndividuals.count() << " seconds" << endl;
            appendTimeComplexityToCSV("time_complexities.csv", "Kruskals", durationInfluentialIndividuals.count());
            break;
        }
        case 12: {
            auto startDijkStra = chrono::high_resolution_clock::now();
            int index = 0;
            cout << "Enter starting Index" << endl;
            cin >> index;
            universityGraph.dijkstra(nodes[index]);

            auto endDijkStra = chrono::high_resolution_clock::now();
            chrono::duration<double> durationInfluentialIndividuals = endDijkStra - startDijkStra;

            cout << "Time taken for Finding Influential Individuals: "
                << durationInfluentialIndividuals.count() << " seconds" << endl;
            appendTimeComplexityToCSV("time_complexities.csv", "Dijkstra", durationInfluentialIndividuals.count());
            break;
        }

        case 13:
        {
            auto startDijkStra = chrono::high_resolution_clock::now();
            do {
                int node_index, event_index;
                cout << "Enter Node index 1-" << MAX_NODES << ": ";
                cin >> node_index;
                cout << "Enter Event index 1-" << MAX_EVENTS << ": ";
                cin >> event_index;

                if (node_index > 0 && node_index <= MAX_NODES && event_index > 0 && event_index <= MAX_EVENTS) {
                    universityGraph.Attendance(nodes[node_index - 1], "Event " + to_string(event_index));
                    cout << "Attendance marked for Node " << node_index << " in Event " << event_index << endl;
                }
                else {
                    cout << "Invalid Node or Event index." << endl;
                }
                option++;
            } while (option < 5);
            cout << "Event Attendance Analysis:" << endl;
            universityGraph.popular_events();
            cout << endl;
            universityGraph.Frequent_participants();
            auto endDijkStra = chrono::high_resolution_clock::now();
            chrono::duration<double> durationInfluentialIndividuals = endDijkStra - startDijkStra;

            cout << "Event Attendance Analysis: "
                << durationInfluentialIndividuals.count() << " seconds" << endl;
            appendTimeComplexityToCSV("time_complexities.csv", "Attendence Analysis", durationInfluentialIndividuals.count());
            break;
        }
        case 14:
            Node * n1;

            int academicImpact;
            int leadershipRoles;
            int involvement;
            cout << "Enter name" << endl;
            cin >> n;
            cout << "Enter role" << endl;
            cin >> r;
            cout << "Enter department" << endl;
            cin >> d;
            cout << "Enter academic impact" << endl;
            cin >> academicImpact;
            cout << "Enter leadership role" << endl;
            cin >> leadershipRoles;
            cout << "Enter involvemrnt" << endl;
            cin >> involvement;
            n1 = new Node(n, r, d, academicImpact, leadershipRoles, involvement);
            universityGraph.add_new_Node(n1);

            break;
        case 15:
            int c;
            cout << "Enter node" << endl;
            cin >> c;
            if (c >= 0 && c < MAX_NODES)
            {
                Node* n1 = nodes[c];
                universityGraph.removeNode(n1);
            }

            else
                cout << "Invalid Index" << endl;
            break;
        case 16:
        {
            auto startDijkStra = chrono::high_resolution_clock::now();
            int index = 0;
            cout << "Enter starting Index" << endl;
            cin >> index;
            universityGraph.bellmanFord(nodes[index]);

            auto endDijkStra = chrono::high_resolution_clock::now();
            chrono::duration<double> durationInfluentialIndividuals = endDijkStra - startDijkStra;

            cout << "Time taken Bellman's Ford: "
                << durationInfluentialIndividuals.count() << " seconds" << endl;
            appendTimeComplexityToCSV("time_complexities.csv", "Bellman's", durationInfluentialIndividuals.count());
            break;
        }
        case 17:
        {
            auto startDijkStra = chrono::high_resolution_clock::now();
            int index = 0;
            universityGraph.detectCommunities();

            auto endDijkStra = chrono::high_resolution_clock::now();
            chrono::duration<double> durationInfluentialIndividuals = endDijkStra - startDijkStra;

            cout << "Time taken Community Detection "
                << durationInfluentialIndividuals.count() << " seconds" << endl;
            appendTimeComplexityToCSV("time_complexities.csv", "Community Detection", durationInfluentialIndividuals.count());
            break;
        }
        default:
            cout << "Invalid choice. Please enter a number from the menu." << endl;
        }
    } while (choice != 0);



    return 0; 
}
