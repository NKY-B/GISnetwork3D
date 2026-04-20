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
typedef bg::model::point<double, 3, bg::cs::cartesian> POINT3D;
typedef bg::model::linestring<POINT3D> LINESTRING3D;
typedef bg::model::point<double, 2, bg::cs::cartesian> POINT2D;
typedef bg::model::linestring<POINT2D> LINESTRING2D;
typedef bg::model::polygon<POINT2D> POLYGON2D;
typedef bg::model::multi_polygon<POLYGON2D> MULTIPOLYGON2D;

//_______________________________________________________________________________________________________________________________
// ---> withinDist3d_p2p_cpp - Find the points within a predefined distance from points

// [[Rcpp::export]] 
DataFrame withinDist3d_p2p_cpp(const DataFrame pt1, const DataFrame pt2, const double dist) {
  
  // Define a value type that stores a point and its ID
  typedef std::pair<POINT2D, int> value;
  
  // Extract coordinates from pt1
  NumericVector x1 = pt1["X"]; // Extract the X coordinate from the input pt
  NumericVector y1 = pt1["Y"]; // Extract the Y coordinate from the input pt
  NumericVector z1 = pt1["Z"]; // Extract the Z coordinate from the input pt
  int n1 = x1.size(); // Get the number of input pt1
  
  // Extract coordinates from pt2
  NumericVector x2 = pt2["X"]; // Extract the X coordinate from the input pt
  NumericVector y2 = pt2["Y"]; // Extract the Y coordinate from the input pt
  NumericVector z2 = pt2["Z"]; // Extract the Z coordinate from the input pt
  int n2 = x2.size(); // Get the number of input pt1
  
  // Build the r-tree with points from pt2
  bgi::rtree<value, bgi::quadratic<16>> rtree;
  for (int i = 0; i < n2; ++i) {
    POINT2D point(x2[i], y2[i]);
    rtree.insert(std::make_pair(point, i + 1)); // Store ID as row + 1
  }
  
  // Prepare the containers to store result
  std::vector<int> id1; // ID of the from point
  std::vector<int> id2; // ID of the to point
  std::vector<double> distances; // Extract the distance
  
  // Declare strategies
  const double buffer_distance = dist * 1.1; // "OGC defines within as completely within and not on the border." So the buffer is expanded to include the on-border points. 
  const int points_per_circle = 36;
  boost::geometry::strategy::buffer::distance_symmetric<double> distance_strategy(buffer_distance);
  boost::geometry::strategy::buffer::join_round join_strategy(points_per_circle);
  boost::geometry::strategy::buffer::end_round end_strategy(points_per_circle);
  boost::geometry::strategy::buffer::point_circle circle_strategy(points_per_circle);
  boost::geometry::strategy::buffer::side_straight side_strategy;
  
  // Perform within-distance search for each point in pt
  for (int i = 0; i < n1; ++i) {
    // Create point buffer
    POINT2D query_point(x1[i], y1[i]); // Define the query point
    boost::geometry::model::multi_polygon<POLYGON2D> Buffer; // Create a container for storing the buffer polygon
    boost::geometry::buffer(query_point, Buffer,
                            distance_strategy, side_strategy,
                            join_strategy, end_strategy, circle_strategy); // Create buffer
    
    // Container to store the neighbours
    std::vector<value> neighbours;
    
    // Query the r-tree for all point within the buffer
    rtree.query(bgi::within(Buffer), std::back_inserter(neighbours));
    
    // Loop through each query point and calculate its distance to each neighbour
    for (int j = 0; j < neighbours.size(); ++j){
      
      // Calculate distance
      POINT3D query_point3D(x1[i], y1[i], z1[i]);
      int neighbourIndex = std::get<1>(neighbours[j]) - 1;
      POINT3D neighbour_point3D(x2[neighbourIndex], y2[neighbourIndex], z2[neighbourIndex]);
      double distance = bg::distance(query_point3D, neighbour_point3D);
      
      // Export data if it is within the distance
      if(distance <= dist){
        id1.push_back(i+1);
        id2.push_back(std::get<1>(neighbours[j]));  
        distances.push_back(distance);}
    }
  }
  // Return the result as a DataFrame
  return DataFrame::create(
    Named("id1") = id1,
    Named("id2") = id2,
    Named("dist") = distances);
}


//_______________________________________________________________________________________________________________________________
// ---> withinDist3d_p_self_cpp - Find the points within a predefined distance

// [[Rcpp::export]] 
DataFrame withinDist3d_p_self_cpp(const DataFrame pt, const double dist) {
  
  // Define a value type that stores a point and its ID
  typedef std::pair<POINT2D, int> value;
  
  // Extract coordinates from pt
  NumericVector x = pt["X"]; // Extract the X coordinate from the input pt
  NumericVector y = pt["Y"]; // Extract the Y coordinate from the input pt
  NumericVector z = pt["Z"]; // Extract the Z coordinate from the input pt
  int n = x.size(); // Get the number of input point
  
  // Build the r-tree with points from pt
  bgi::rtree<value, bgi::quadratic<16>> rtree;
  for (int i = 0; i < n; ++i) {
    POINT2D point(x[i], y[i]);
    rtree.insert(std::make_pair(point, i + 1)); // Store ID as row + 1
  }
  
  // Prepare the containers to store result
  std::vector<int> id1; // ID of the from point
  std::vector<int> id2; // ID of the to point
  std::vector<double> distances; // Extract the distance
  
  // Declare strategies
  const double buffer_distance = dist * 1.1; // "OGC defines within as completely within and not on the border." So the buffer is expanded to include the on-border points. 
  const int points_per_circle = 36;
  boost::geometry::strategy::buffer::distance_symmetric<double> distance_strategy(buffer_distance);
  boost::geometry::strategy::buffer::join_round join_strategy(points_per_circle);
  boost::geometry::strategy::buffer::end_round end_strategy(points_per_circle);
  boost::geometry::strategy::buffer::point_circle circle_strategy(points_per_circle);
  boost::geometry::strategy::buffer::side_straight side_strategy;
  
  // Perform within-distance search for each point in pt
  for (int i = 0; i < n; ++i) {
    // Create point buffer
    POINT2D query_point(x[i], y[i]); // Define the query point
    boost::geometry::model::multi_polygon<POLYGON2D> Buffer; // Create a container for storing the buffer polygon
    boost::geometry::buffer(query_point, Buffer,
                            distance_strategy, side_strategy,
                            join_strategy, end_strategy, circle_strategy); // Create buffer
    
    // Container to store the neighbours
    std::vector<value> neighbours;
    
    // Query the r-tree for all point within the buffer
    rtree.query(bgi::within(Buffer), std::back_inserter(neighbours));
    
    // Loop through each query point and calculate its distance to each neighbour
    for (int j = 0; j < neighbours.size(); ++j){
      
      // Calculate distance
      POINT3D query_point3D(x[i], y[i], z[i]);
      int neighbourIndex = std::get<1>(neighbours[j]) - 1;
      POINT3D neighbour_point3D(x[neighbourIndex], y[neighbourIndex], z[neighbourIndex]);
      double distance = bg::distance(query_point3D, neighbour_point3D);
      
      // Export data if it is within the distance
      if(distance <= dist && (i+1) != (neighbourIndex + 1)){
        id1.push_back(i+1);
        id2.push_back(std::get<1>(neighbours[j]));  
        distances.push_back(distance);}
      }
    }
  // Return the result as a DataFrame
  return DataFrame::create(
    Named("id1") = id1,
    Named("id2") = id2,
    Named("dist") = distances);
}


//_______________________________________________________________________________________________________________________________
// ---> withinDist3d_l2p_cpp - Find the points within a predefined distance from linestring

// [[Rcpp::export]] 
DataFrame withinDist3d_l2p_cpp(const DataFrame lines_df, const DataFrame pt_df, const double dist) {
  
  // Define a value type that stores a point and its ID
  typedef std::pair<POINT2D, int> value;
  
  // Extract coordinates from linestring
  NumericVector lx = lines_df["X"]; // Extract the X coordinate from the input linestring
  NumericVector ly = lines_df["Y"]; // Extract the Y coordinate from the input linestring
  NumericVector lz = lines_df["Z"]; // Extract the Z coordinate from the input linestring
  NumericVector lL1 = lines_df["L1"]; // Extract the L1 from linestring
  
  //___________________________________
  // ----->>> Create 2D linestring
  // Create a map to store 2D linestrings by ID
  std::map<int, LINESTRING2D> linestrings2D;
  
  // Loop through the points and group them by ID
  for (size_t i = 0; i < lL1.size(); i++) {
    int line_id = lL1[i];
    double x_coord = lx[i];
    double y_coord = ly[i];
    
    // Add points to the corresponding 2D linestring by line ID
    linestrings2D[line_id].push_back(POINT2D(x_coord, y_coord));
  }
  
  //___________________________________
  // ----->>> Create 3D linestring
  // Create a map to store 3D linestrings by ID
  std::map<int, LINESTRING3D> linestrings3D;
  
  // Loop through the points and group them by ID
  for (size_t i = 0; i < lL1.size(); i++) {
    int line_id = lL1[i];
    double x_coord = lx[i];
    double y_coord = ly[i];
    double z_coord = lz[i];
    
    // Add points to the corresponding 2D linestring by line ID
    linestrings3D[line_id].push_back(POINT3D(x_coord, y_coord, z_coord));
  }
  
  //___________________________________
  // ----->>> Build Spatial Index for the input point
  
  // Extract coordinates from pt
  NumericVector px = pt_df["X"]; // Extract the X coordinate from the input pt
  NumericVector py = pt_df["Y"]; // Extract the Y coordinate from the input pt
  NumericVector pz = pt_df["Z"]; // Extract the Z coordinate from the input pt
  int n = px.size(); // Get the number of input pt
  
  // Build the r-tree with points from pt
  bgi::rtree<value, bgi::quadratic<16>> rtree;
  for (int i = 0; i < n; ++i) {
    POINT2D point(px[i], py[i]);
    rtree.insert(std::make_pair(point, i + 1)); // Store ID as row + 1
  }
  
  //___________________________________
  // ----->>> Search point for each linestring
  
  // Prepare the containers to store result
  std::vector<int> L1; // ID of the linestring
  std::vector<int> id; // ID of the point
  std::vector<double> distances; // Extract the distance
  
  // Declare strategies
  const double buffer_distance = dist * 1.1; // "OGC defines within as completely within and not on the border." So the buffer is expanded to include the on-border points. 
  const int points_per_circle = 36;
  boost::geometry::strategy::buffer::distance_symmetric<double> distance_strategy(buffer_distance);
  boost::geometry::strategy::buffer::join_round join_strategy(points_per_circle);
  boost::geometry::strategy::buffer::end_round end_strategy(points_per_circle);
  boost::geometry::strategy::buffer::point_circle circle_strategy(points_per_circle);
  boost::geometry::strategy::buffer::side_straight side_strategy;
  
  // Perform within-distance search for each line
  for (int i = 0; i < linestrings2D.size(); ++i) {
    // Create line buffer
    LINESTRING2D query_line = linestrings2D[i+1];
    boost::geometry::model::multi_polygon<POLYGON2D> Buffer; // Create a container for storing the buffer polygon
    boost::geometry::buffer(query_line, Buffer,
                            distance_strategy, side_strategy,
                            join_strategy, end_strategy, circle_strategy); // Create buffer
    
    // Container to store the neighbours
    std::vector<value> neighbours;
    
    // Query the r-tree for all point within the buffer
    rtree.query(bgi::within(Buffer), std::back_inserter(neighbours));
    
    // Loop through each query line and calculate its distance to each neighbour
    for (int j = 0; j < neighbours.size(); ++j){
      
      // Calculate distance
      LINESTRING3D query_line3D = linestrings3D[i+1]; // Index line by the lingstring ID
      int neighbourIndex = std::get<1>(neighbours[j]) - 1; // Get the positional index for the point
      POINT3D neighbour_point3D(px[neighbourIndex], py[neighbourIndex], pz[neighbourIndex]); // Create the destination point
      double distance = bg::distance(query_line3D, neighbour_point3D); // Calculate distance
      
      // Export data if it is within the distance
      if(distance <= dist){
        L1.push_back(i+1);
        id.push_back(std::get<1>(neighbours[j]));  
        distances.push_back(distance);}
    }
  }
  // Return the result as a DataFrame
  return DataFrame::create(
    Named("L1") = L1,
    Named("PTid") = id,
    Named("dist") = distances);
}

