// application.cpp <Starter Code>
// <Riya Patel>
//
//
// Adam T Koehler, PhD
// University of Illinois Chicago
// CS 251, Fall 2023
//
// Project Original Variartion By:
// Joe Hummel, PhD
// University of Illinois at Chicago
//
// 
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <queue>
#include "tinyxml2.h"
#include "dist.h"
#include "graph.h"
#include "osm.h"
#include <stack>




using namespace std;
using namespace tinyxml2;


const double INF = std::numeric_limits<double>::max();

//Searchs for a building within a vector of buildings based on a query string
//Parameters include buildings, a vector containing all the buildings and a query which is an abbreviation
//Return am empty building info if no buildings are found
BuildingInfo searchBuilding(
  const vector<BuildingInfo>& buildings, const string& query){

  //iterate through each building matching its abbreviation
  for(const auto& building : buildings){
      if(building.Abbrev == query){
          return building;
      }
  }

  //iterate through each building matching its full name
  for(const auto& building : buildings){
    if(building.Fullname.find(query) != string::npos){
      return building;
    }
  }
  return BuildingInfo(); 
}

//Finds the building in a provided vector closest to the midpoint betweem two buildings
//parameters include, the two buildings you want to find the middle building of and a vector of all the buildings
//Returns the buildingInfo object representing the nearest building
BuildingInfo findNearestBuildingToMidpoint(
    const BuildingInfo& building1, const BuildingInfo& building2, const vector<BuildingInfo>& Buildings){
    
    //calculate the coordinates of the midpoint between the buildings
    Coordinates midpoint = centerBetween2Points(building1.Coords.Lat, building1.Coords.Lon, building2.Coords.Lat, building2.Coords.Lon);
  
    BuildingInfo destinationBuilding;
    double minDistance = INF;

    //iterate through each building in the vector
    for(const auto& building : Buildings){
      //calculate the distance between the current building and the midpoint
      double distance = distBetween2Points(midpoint.Lat, midpoint.Lon, building.Coords.Lat, building.Coords.Lon);
      //if the distance is smaller, replcate the building
      if(distance < minDistance){
        minDistance = distance;
        destinationBuilding = building;
        
      }
    }
    return destinationBuilding;
}

//Finds the ID of the nearest building node to a given coordinate among footways
//Parameters include a coordinate object, a vector containing FootwayInfo objects, a map containing node IDs
//Returns a long long ID of the nearest building node to the given coordinate
long long findNearestBuilding(
  const Coordinates& building, const vector<FootwayInfo>& Footways, const map<long long, Coordinates>& Nodes){

  double minDistance = INF;
  long long nearestBuildingID = -1; 

  //iterate through each footway
  for(const auto& footway : Footways){
      //iterate through each footway
      for (const long long& nodeId : footway.Nodes){
        //retrieve the coordinates
          const Coordinates& loc = Nodes.at(nodeId);
          double distance = distBetween2Points(building.Lat, building.Lon, loc.Lat, loc.Lon);

          //check if the distance to the current node is smaller than the minimum distance
          if(distance < minDistance){
              minDistance = distance;
              nearestBuildingID = nodeId; 
          }
      }
  }
    return nearestBuildingID;
  }


//Prints the information of a building using its ID
//A long long ID of the building ndoe, a map containing node IDs mapped to their coordinates
void printBuilding(
  long long buildingID, map<long long, Coordinates>& Nodes){

  //attempt to find the node corresponding to the building ID
  auto nodeIt = Nodes.find(buildingID);

  //Check if the node corresponding to the building ID is found
  if(nodeIt != Nodes.end()){
      const Coordinates& buildingCoords = nodeIt->second;
      //print out the information of the building ID
      cout << " " << buildingID << endl;
      cout << " (" << buildingCoords.Lat << ", " << buildingCoords.Lon << ")" << endl;
  }
}

//Finds the shortest path in a weighted graph froma starting vertex to all other vertices 
//Parameters include graph, verticies, starting/ending vertex. distance map, pathNode map
template<typename VertexT, typename WeightT>
double Dijkstra(
  const graph<VertexT, WeightT>& G, VertexT startV, VertexT endV, map<VertexT, double>& distance, map<VertexT, VertexT>& pathNode){
  
    map<VertexT, bool> visited;

    //Initialize distances and pathNode for each vertex in the graph
    for(const VertexT& v : G.getVertices()){
        distance[v] = INF;
        pathNode[v] = VertexT{}; //default vertexT value
        visited[v] = false;
    }

    //set the distance of starting vertex
    distance[startV] = 0;
    while(true){
        VertexT curr;
        double minDist = INF;

        //Find the vertex with the smallest distance among unvisited vertices
        for(const VertexT& vertex : G.getVertices()){
            if(!visited[vertex] && distance[vertex] < minDist){
              minDist = distance[vertex];
              curr = vertex;
            }
        }

        //if all vertices have been visited, break the loop
        if(minDist == INF){
          break;
        }

        visited[curr] = true;

        //update distances to neighboring vertices through the current vertex
        for(const VertexT& neighbor : G.neighbors(curr)){
            WeightT weight;
            G.getWeight(curr, neighbor, weight);

            //calculate the new distance to neighbor via curr, update the distance if a shorter path is found
            double newNeighbor = distance[curr] + weight;
            if(newNeighbor < distance[neighbor]){
                distance[neighbor] = newNeighbor;
                pathNode[neighbor] = curr;
            }
        }

        //If the end vertex is reached, return its distance
        if(curr == endV){
          return distance[endV];
        }
    } 
    return INF;
}

//Reconstructs and prints the path between two vertices based on the pathNode
//Parameters include a map, startVertex and endVertex
void getPath(
  map<long long, long long>& pathNode, long long startVertex, long long endVertex){
  
  //If the start and end vertex are the same, there is no path, so print the single vertex
  if(startVertex == endVertex){
      cout << "Path: " << startVertex << endl;
      return;
  }

  //initialize data structures for storing and printing the path
  vector<long long> path;
  stack<long long> vertexStack;

  //push the endVertex onto stack
  vertexStack.push(endVertex);

  //Reconstruct the path by traversing the pathNode map backwards
  while (pathNode[endVertex] != startVertex && pathNode[endVertex] != 0){
      //move the previous vertex and psuh the vertex onto stacl
      endVertex = pathNode[endVertex];
      vertexStack.push(endVertex);
  }

  //push the start vertex on to stack
  vertexStack.push(startVertex);

  cout << "Path: ";
  bool isFirst = true;
  //print the reconstructed vertices from the stack
  while (!vertexStack.empty()){
      if(isFirst){
          cout << vertexStack.top();
          isFirst = false;
      } else {
          cout << "->" << vertexStack.top();
      }
      vertexStack.pop();
  }
  cout << endl;
}

//Runs the map and lets user input different locations
//Parameters include map Nodes, vector of footways, vector of buildings, graph
void application(
    map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways,
    vector<BuildingInfo>& Buildings, graph<long long, double>& G){
    string person1Building, person2Building;

/*-----------------------------------------------------------------------------------------*/
    
    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);

    while (person1Building != "#"){

      cout << endl;
      cout << "Enter person 2's building (partial name or abbreviation)> ";
      getline(cin, person2Building);

      //searches the building that the user inputs
      BuildingInfo building1 = searchBuilding(Buildings, person1Building);
      //if building info is not found, ask for another input
      if(building1.Fullname.empty()){
          cout << "Person 1's building not found" << endl;
          cout << "Enter person 1's building (partial name or abbreviation), or #> ";
          getline(cin, person1Building);
          continue;
      }
        
      //searches for building 2, if building info is not found, ask for another input
      BuildingInfo building2 = searchBuilding(Buildings, person2Building);
      if(building2.Fullname.empty()){
          cout << "Person 2's building not found" << endl;
          cout << "Enter person 1's building (partial name or abbreviation), or #> ";
          getline(cin, person1Building);
          continue;
      }
      //checks if they want to end the program
      if(person1Building == "#" && person2Building == "#"){
        break;
      }

/*-----------------------------------------------------------------------------------------*/

      //prints the building name, latitude and longitude 
      cout << endl;
      cout << "Person 1's point:" << endl;
      cout << " " << building1.Fullname  << endl;
      cout << " (" << building1.Coords.Lat << ", " << building1.Coords.Lon << ")" << endl;

      cout << "Person 2's point:" << endl;
      cout << " " << building2.Fullname <<  endl;
      cout << " (" << building2.Coords.Lat << ", " << building2.Coords.Lon << ")" << endl;

      BuildingInfo destinationBuilding = findNearestBuildingToMidpoint(building1, building2, Buildings);
     
      cout << "Destination Building:" << endl;
      cout << " " << destinationBuilding.Fullname  << endl;
      cout << " (" << destinationBuilding.Coords.Lat << ", " << destinationBuilding.Coords.Lon << ")" << endl;
  

/*-----------------------------------------------------------------------------------------*/  
      

      // Find nearest nodes to each building
      cout << endl;
      long long nearestToBuilding1 = findNearestBuilding(building1.Coords, Footways, Nodes);
      cout << "Nearest P1 node:" << endl;
      printBuilding(nearestToBuilding1, Nodes);

      long long nearestToBuilding2 = findNearestBuilding(building2.Coords, Footways, Nodes);
      cout << "Nearest P2 node:" << endl;
      printBuilding(nearestToBuilding2, Nodes);
      long long nearestToDestinationBuilding = findNearestBuilding(destinationBuilding.Coords, Footways, Nodes);
      
      //If the building is not found, then report that they were unable to reach the destination
      //Print the nearest destionation is found
      if(nearestToDestinationBuilding < INF){
        cout << "Nearest destination node:" << endl;
        printBuilding(nearestToDestinationBuilding, Nodes);
        cout << endl;
      }else{
        cout << "At least one person was unable to reach the destination building. Finding next closest building..." << endl;
      }


/*-----------------------------------------------------------------------------------------*/   
   
    map<long long, double> distance1, distance2;
    map<long long, long long> pathNode1, pathNode2;

    // Call Dijkstra's algorithm
    double shortDistance1 = Dijkstra(G, nearestToBuilding1, nearestToDestinationBuilding, distance1, pathNode1);
    double shortDistance2 = Dijkstra(G, nearestToBuilding2, nearestToDestinationBuilding, distance2, pathNode2);


    //If the distance is INF, report that the destionation is unreachable
    if(shortDistance1 == INF || shortDistance2 == INF){
        cout << "Sorry, destination unreachable." << endl;
    }else{
        //print the path if there is one
        cout << "Person 1's distance to dest: " << shortDistance1 << " miles" << endl;
        getPath(pathNode1, nearestToBuilding1, nearestToDestinationBuilding);
        cout << endl;

        cout << "Person 2's distance to dest: " << shortDistance2 << " miles" << endl;
        getPath(pathNode2, nearestToBuilding2, nearestToDestinationBuilding);
        cout << endl;
    }

    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  }
}


int main() {
  graph<long long, double> G;

  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates>  Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo>          Footways;
  // info about each building, in no particular order
  vector<BuildingInfo>         Buildings;
  XMLDocument                  xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;



  //
  // TO DO: build the graph, output stats:
  //
  for(const auto& node : Nodes){
    G.addVertex(node.first);
  }

 
  for(const auto& footway : Footways){
    const auto& nodes = footway.Nodes;
          
    for(size_t i = 0; i < nodes.size() - 1; i++){
        long long fromID = nodes[i];
        long long toID = nodes[i + 1];

        const auto& fromCoord = Nodes[fromID];
        const auto& toCoord = Nodes[toID];

        double distance = distBetween2Points(fromCoord.Lat, fromCoord.Lon, toCoord.Lat, toCoord.Lon); 
      
        G.addEdge(fromID, toID, distance);
        G.addEdge(toID, fromID, distance);
    }
  }

  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  // Execute Application
  application(Nodes, Footways, Buildings, G);

  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}
