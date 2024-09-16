#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>

using namespace std;

template<typename VertexT, typename WeightT>
class graph {
private:
    map<VertexT, vector<pair<VertexT, WeightT>>> adjacencyList;

public:
    //Default constructor
    graph(){}

    //Returns the number of verticies in the graph, considering each unique vertex
    int NumVertices() const{
        return adjacencyList.size();
    }

    //Computes the total numbers of edges in the graph
    //Returns the total number of edges
    int NumEdges() const {
        int totalEdges = 0;
        //iterates through each vertex in the adjacency list
        //counts the number of edges connected to each vertex
        for(const auto& vertex : adjacencyList){
            totalEdges += vertex.second.size();
        }
        return totalEdges;
    }

    //Adds a new vertex to the graph represented by the adjacency list
    //Parameter includes the vertex to be added
    //Returns true if the vertex is added successfully, false otherwise
    bool addVertex(VertexT vertex) {
        //Searches for the vertex in the adjacency list by having a pointer at the last position 
        if(adjacencyList.find(vertex) != adjacencyList.end()){
            return false;
        }
        //Adds a new vertex to the graph by creating an entry in the adjacency list for that vertex with an empty set
        adjacencyList[vertex] = {};
        return true;
    }


    //Adds a new edge to the vertex
    //Parameters include 'from' is the source vertex, 'to' is the destination vertex, 'weight' is the output parameter for the edge weight
    //Returns true if the edge is added successfule, false otherwise
    bool addEdge(VertexT from, VertexT to, WeightT weight) {
        //Check if both vertices exist in the graph
        if(adjacencyList.find(from) == adjacencyList.end() || adjacencyList.find(to) == adjacencyList.end()){
            return false;
        }

        //Gets edges associates with 'from'
        auto& edges = adjacencyList[from];

        //Checks if an edge (from, to) exists within the vertex's edges
        auto edgeIt = edges.end();
        for(auto it = edges.begin(); it != edges.end(); ++it){
            if(it->first == to){
                edgeIt = it;
                break;
            }
        }

        if(edgeIt != edges.end()){
            //Edge already exists, update the weight
            edgeIt->second = weight;
        }else{
            //Edge doesn't exist, add the new edge
            edges.push_back({ to, weight });
        }
        return true;
    }

    //Retrieves the weight of an edge between two vertices if they exist
    //Parameters include 'from' is the source vertex, 'to' is the destination vertex, 'weight' is the output parameter for the edge weight
    //Returns true if the edge exists and its weight is stored in 'weight', false otherwise
    bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
        //Find the source vertex in the adjacency list
        auto iter = adjacencyList.find(from);
        //If the source vertex exists
        if(iter != adjacencyList.end()){
            //Iterate through its edges
            for (const auto& edge : iter->second){
                //If an edge to the destination vertex is found, store the weight
                if(edge.first == to){
                    weight = edge.second;
                    return true;
                }
            }
        }
        return false;
    }

    //Retrieves a set containing all neighbors of a given vertex
    //Parameters includes the vertex for which neighbors are to be found
    //Returns a set containing all neighbots vertices of vertex
    set<VertexT> neighbors(VertexT vertex) const {
        //Initializes an empty set to store neighbors
        set<VertexT> neighbors;
        //Checks if the vertex exists in the adjacency list
        if(adjacencyList.find(vertex) != adjacencyList.end()){
            //Iterates through the edges of vertex 
            for (const auto& neighbor : adjacencyList.at(vertex)){
                //If it is a neighbor, insert the vertex into the set
                neighbors.insert(neighbor.first);
            }
        }
        return neighbors;
    }

    //Retrieves all vertices present in the graph
    //Returns a vector containing all the vertices in the graph
    vector<VertexT> getVertices() const{
        //Initializes an empty vector to store verticies
        vector<VertexT> vertices;
        //Iterate through the adjacency list
        for(const auto& vertex : adjacencyList){
            //add each vertex to the vector
            vertices.push_back(vertex.first);
        }
        //Return the vector containing all verticies
        return vertices;
    }

    //Dumps the internal state of the graph
    //Outputs the number of vertices and lists each vertex with its corresponding edges
    //Displays the number of vertices and their connected edges
    void dump(ostream& output) const {
        output << "**************************************************" << endl;
        output << "********************* GRAPH ***********************" << endl;
        output << "**Num vertices: " << adjacencyList.size() << endl;
        output << "**Vertices:" << endl;

        //Outputs each vertex in the graph
        for (const auto& vertex : adjacencyList) {
            output << " " << vertex.first << endl;
        }
        output << "**Edges:" << endl;
        //Outputs each edge for each vertex
        for (const auto& vertex : adjacencyList) {
            output << " " << vertex.first << ": ";
            for (const auto& neighbor : vertex.second) {
                output << "(" << neighbor.first << "," << neighbor.second << ") ";
            }
            output << endl;
        }
        output << "**************************************************" << endl;
    }
};


