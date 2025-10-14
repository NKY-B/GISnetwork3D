// [[Rcpp::depends(BH)]]
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <Rcpp.h>
#include <vector>
#include <iostream>
#include <boost/foreach.hpp>
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using namespace Rcpp;


// Define data type: <POINT3D> and <LINESTRING3D>
typedef bg::model::point<double, 3, bg::cs::cartesian> POINT3D;
typedef bg::model::linestring<POINT3D> LINESTRING3D;
typedef bg::model::segment<POINT3D> SEGMENT3D;

//_______________________________________________________________________________________________________________________________
// ---> snap3d_p2p_cpp - Change the coordinate of point x to the nearest coordinate of point y if they are within a threshold distance

// [[Rcpp::export]] 
DataFrame snap3d_p2p_cpp(const DataFrame pt1, const DataFrame pt2, const double tolerance) {
  
  // Define a value type that stores a point and its ID
  typedef std::pair<POINT3D, int> value;
  
  // Extract coordinates from pt2
  NumericVector x2 = pt2["X"];
  NumericVector y2 = pt2["Y"];
  NumericVector z2 = pt2["Z"];
  int n2 = x2.size();
  
  // Build the r-tree with points from pt2
  bgi::rtree<value, bgi::quadratic<16>> rtree;
  for (int i = 0; i < n2; ++i) {
    POINT3D point(x2[i], y2[i], z2[i]);
    rtree.insert(std::make_pair(point, i + 1)); // Store ID as row + 1
  }
  
  // Extract coordinates from pt1
  NumericVector x1 = pt1["X"];
  NumericVector y1 = pt1["Y"];
  NumericVector z1 = pt1["Z"];
  
  // Find the row number of the points from pt1
  int n1 = x1.size();
  
  // Prepare the containers to store result
  IntegerVector result(n1);
  std::vector<double> X; // To store the final X coordinate
  std::vector<double> Y; // To store the final Y coordinate
  std::vector<double> Z; // To store the final Z coordinate
  std::vector<double> dist; // Extract the distance
  std::vector<int> did; // Extract nearest pt ID
  
  // Perform KNN search for each point in pt1
  for (int i = 0; i < n1; ++i) {
    // Define parameters
    POINT3D query_point(x1[i], y1[i], z1[i]); // Create query point
    std::vector<value> neighbours; // Create neighbours to store the KNN result
    
    // Query the r-tree for k nearest neighbours
    rtree.query(bgi::nearest(query_point, 1), std::back_inserter(neighbours)); // Find the nearest point
    POINT3D nearest_point = std::get<0>(neighbours[0]); // Extract the nearest point geometry
    int nearest_point_id = std::get<1>(neighbours[0]); // Extract the nearest point ID
    
    // Calculate the distance between the query point and the nearest point
    double distance = bg::distance(query_point, nearest_point);
    
    if(distance <= tolerance){
      X.push_back(nearest_point.get<0>());
      Y.push_back(nearest_point.get<1>());
      Z.push_back(nearest_point.get<2>());
      dist.push_back(distance);
      did.push_back(nearest_point_id);
    } else{
      X.push_back(x1[i]);
      Y.push_back(y1[i]);
      Z.push_back(z1[i]);
      dist.push_back(distance);
      did.push_back(nearest_point_id);
    }
    
  }
  
  // Return the result as a DataFrame
  return DataFrame::create(
    Named("X") = X,
    Named("Y") = Y,
    Named("Z") = Z,
    Named("did") = did,
    Named("dist") = dist);
}


//_______________________________________________________________________________________________________________________________
// ---> snap3d_p2l_cpp - Change the coordinate of point x to the nearest coordinate on line segment y if they are within a threshold distance

// [[Rcpp::export]] 
DataFrame snap3d_p2l_cpp(const DataFrame pt, const DataFrame segments, const double tolerance) {
  
  // Define a value type that stores a segment, its line ID, and its segment ID
  typedef std::tuple<SEGMENT3D, int, int> value; // SEGMENT3D, line ID, segment ID
  
  // Extract coordinates and line IDs from segment
  NumericVector x1 = segments["from_x"]; // X coordinate of the start node
  NumericVector y1 = segments["from_y"]; // Y coordinate of the start node
  NumericVector z1 = segments["from_z"]; // Z coordinate of the start node
  NumericVector x2 = segments["to_x"]; // X coordinate of the end node
  NumericVector y2 = segments["to_y"]; // Y coordinate of the end node
  NumericVector z2 = segments["to_z"]; // Z coordinate of the end node
  IntegerVector L1 = segments["L1"]; // Linestring ID of the segment
  IntegerVector L2 = segments["L2"]; // Segment ID of the segment
  int n_segments = x1.size(); // Get the total number of segments
  
  // Build the r-tree with segments
  bgi::rtree<value, bgi::quadratic<16>> rtree;
  for (int i = 0; i < n_segments; ++i) {
    POINT3D p1(x1[i], y1[i], z1[i]);
    POINT3D p2(x2[i], y2[i], z2[i]);
    SEGMENT3D segment(p1, p2);
    rtree.insert(std::make_tuple(segment, L1[i], L2[i])); // Store L1 and L2
  }
  
  // Extract coordinates from query points
  NumericVector px = pt["X"]; // X coordinate of the pt
  NumericVector py = pt["Y"]; // Y coordinate of the pt
  NumericVector pz = pt["Z"]; // Z coordinate of the pt
  int n_points = px.size(); // Get the total number of points
  
  // Prepare the containers to store result
  std::vector<double> X; // To store the final X coordinate
  std::vector<double> Y; // To store the final Y coordinate
  std::vector<double> Z; // To store the final Z coordinate
  std::vector<double> dist; // Extract the distance
  std::vector<int> dL1; // Extract nearest linestring ID
  std::vector<int> dL2; // Extract nearest segment ID
  
  // Perform KNN search for each point in pt
  for (int i = 0; i < n_points; ++i) {
    // Define parameters
    POINT3D query_point(px[i], py[i], pz[i]); // Create query point
    std::vector<value> neighbours; // Create neighbours to store the KNN result
    
    // Query the r-tree for k nearest neighbours
    rtree.query(bgi::nearest(query_point, 1), std::back_inserter(neighbours)); // Find the nearest segment
    SEGMENT3D nearest_segment = std::get<0>(neighbours[0]); // Extract the nearest segment geometry
    int nearest_linestring_id = std::get<1>(neighbours[0]); // Extract the nearest linestring ID
    int nearest_segment_id = std::get<2>(neighbours[0]); // Extract the nearest segment ID
    
    // Project point to line
    // Use Boost.Geometry to calculate the closest segment linking the point and the linestring
    bg::model::segment<POINT3D> nearest_seg; // Define segment container
    bg::closest_points(query_point, nearest_segment, nearest_seg); // Generate segment
    POINT3D nearest_point = nearest_seg.second; // Extract the closest point along the segment
    double distance = bg::distance(query_point, nearest_point); // Calculate the distance between the query point and the nearest point
    
    if(distance <= tolerance){
      // If the distance is within the tolerance distance, replace the coordinates of the input pt to the projected coordinates along segment
      X.push_back(nearest_point.get<0>());
      Y.push_back(nearest_point.get<1>());
      Z.push_back(nearest_point.get<2>());
      dist.push_back(distance);
      dL1.push_back(nearest_linestring_id);
      dL2.push_back(nearest_segment_id);
      
    } else{
      // If the distance is NOT within the tolerance distance, retain the original coordinates of the input point
      X.push_back(px[i]);
      Y.push_back(py[i]);
      Z.push_back(pz[i]);
      dist.push_back(distance);
      dL1.push_back(nearest_linestring_id);
      dL2.push_back(nearest_segment_id);
    }
    
  }
  
  // Return the result as a DataFrame
  return DataFrame::create(
    Named("X") = X,
    Named("Y") = Y,
    Named("Z") = Z,
    Named("L1") = dL1,
    Named("L2") = dL2,
    Named("dist") = dist);
}

