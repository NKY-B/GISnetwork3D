// [[Rcpp::depends(BH)]]
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <Rcpp.h>
#include <vector>
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using namespace Rcpp;


// Define data type
typedef bg::model::point<double, 2, bg::cs::cartesian> POINT2D;
typedef bg::model::linestring<POINT2D> LINESTRING2D;
typedef bg::model::polygon<POINT2D> POLYGON2D;
typedef bg::model::multi_polygon<POLYGON2D> MULTIPOLYGON2D;

//_______________________________________________________________________________________________________________________________
// ---> tilePoint_cpp - Find the points within tile

// [[Rcpp::export]] 
DataFrame tilePoint_cpp(const DataFrame tile_df, const DataFrame pt) {
  
  // Define a value type that stores a point and its ID
  typedef std::pair<POINT2D, int> value;
  
  // Extract coordinates from tile
  NumericVector tile_x = tile_df["X"]; // Extract the X coordinate from the input tile
  NumericVector tile_y = tile_df["Y"]; // Extract the Y coordinate from the input tile
  NumericVector tile_L2 = tile_df["L2"]; // Extract the L2 from tile
  
  // Create a map to store tiles by L2
  std::map<int, POLYGON2D> tiles;
  
  // Loop through the points and group them by L2
  for (size_t i = 0; i < tile_x.size(); i++) {
    int L2 = tile_L2[i];
    double x = tile_x[i];
    double y = tile_y[i];
    
    // Add points to the corresponding tile
    bg::append(tiles[L2].outer(), POINT2D(x, y));
  }
  
  // Extract coordinates from pt
  NumericVector ptx = pt["X"]; // Extract the X coordinate from the input pt
  NumericVector pty = pt["Y"]; // Extract the Y coordinate from the input pt
  
  // Build the r-tree with points from pt
  bgi::rtree<value, bgi::quadratic<16>> rtree;
  for (int i = 0; i < ptx.size(); ++i) {
    POINT2D point(ptx[i], pty[i]);
    rtree.insert(std::make_pair(point, i + 1)); // Store ID as row + 1
  }
  
  // Prepare the containers to store result
  std::vector<int> id1; // ID of the tile
  std::vector<int> id2; // ID of the point
  
  // Find all points within each tile
  for (int i = 0; i < tiles.size(); ++i) {
    // Container to store the neighbours
    std::vector<value> neighbours;
    
    // Query the r-tree for all point within the buffer
    rtree.query(bgi::intersects(tiles[i+1]), std::back_inserter(neighbours));
    
    // Loop through each neighbour to store the id
    for (int j = 0; j < neighbours.size(); ++j){
      id1.push_back(i + 1);
      id2.push_back(neighbours[j].second); // Access the ID from the pair
    }
    
    if (neighbours.size() == 0) {
      id1.push_back(i + 1);
      id2.push_back(-999);
    }
  }
  
  // Return the result as a DataFrame
  return DataFrame::create(
    Named("gridID") = id1,
    Named("id") = id2);
}




//_______________________________________________________________________________________________________________________________
// ---> tileNode_cpp - Find the points within tile

// [[Rcpp::export]] 
DataFrame tileNode_cpp(const DataFrame tile_df, const DataFrame pt, const std::vector<std::string>& nodeID) {
  
  // Define a value type that stores a point, its ID and the nodeID
  typedef std::tuple<POINT2D, int, std::string> value;
  
  // Extract coordinates from tile
  NumericVector tile_x = tile_df["X"]; // Extract the X coordinate from the input tile
  NumericVector tile_y = tile_df["Y"]; // Extract the Y coordinate from the input tile
  NumericVector tile_L2 = tile_df["L2"]; // Extract the L2 from tile
  
  // Create a map to store tiles by L2
  std::map<int, POLYGON2D> tiles;
  
  // Loop through the points and group them by L2
  for (size_t i = 0; i < tile_x.size(); i++) {
    int L2 = tile_L2[i];
    double x = tile_x[i];
    double y = tile_y[i];
    
    // Add points to the corresponding tile
    bg::append(tiles[L2].outer(), POINT2D(x, y));
  }
  
  // Extract coordinates from pt
  NumericVector ptx = pt["X"]; // Extract the X coordinate from the input pt
  NumericVector pty = pt["Y"]; // Extract the Y coordinate from the input pt
  
  // Build the r-tree with points from pt
  bgi::rtree<value, bgi::quadratic<16>> rtree;
  for (int i = 0; i < ptx.size(); ++i) {
    POINT2D point(ptx[i], pty[i]);
    rtree.insert(std::make_tuple(point, i + 1, nodeID[i])); // Store ID as row + 1, and store the nodeID
  }
  
  // Prepare the containers to store result
  std::vector<int> id1; // ID of the tile
  std::vector<int> id2; // ID of the point
  std::vector<std::string> nID; // ID of the point
  
  // Find all points within each tile
  for (int i = 0; i < tiles.size(); ++i) {
    // Container to store the neighbours
    std::vector<value> neighbours;
    
    // Query the r-tree for all point within the buffer
    rtree.query(bgi::intersects(tiles[i+1]), std::back_inserter(neighbours));
    
    // Loop through each neighbour to store the id
    for (int j = 0; j < neighbours.size(); ++j){
      id1.push_back(i+1);
      id2.push_back(std::get<1>(neighbours[j]));
      nID.push_back(std::get<2>(neighbours[j]));
    }
    
    if(neighbours.size() == 0){
      id1.push_back(i+1);
      id2.push_back(-999);
      nID.push_back("-999");
    }
    
  }
  
  // Return the result as a DataFrame
  return DataFrame::create(
    Named("gridID") = id1,
    Named("id") = id2,
    Named("nID") = nID);
}




