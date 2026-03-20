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
// ---> nearest_p2p_cpp - Nearest object for point to point

// [[Rcpp::export]] 
IntegerVector nearest3d_p2p_cpp(const DataFrame pt1, const DataFrame pt2) {
  
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
  
  // Prepare the result list
  IntegerVector result(n1);
  
  // Perform KNN search for each point in pt1
  for (int i = 0; i < n1; ++i) {
    POINT3D query_point(x1[i], y1[i], z1[i]);
    std::vector<value> neighbours;
    
    // Query the r-tree for k nearest neighbours
    rtree.query(bgi::nearest(query_point, 1), std::back_inserter(neighbours));
    
    // Collect the ID of the neighbours
    result[i] = neighbours[0].second; // Extract the ID of the nearest neighbour
    
  }
  
  return result;
}


//_______________________________________________________________________________________________________________________________
// ---> nearest_p2l_cpp - Nearest object for point to line

// [[Rcpp::export]] 
IntegerVector nearest3d_p2l_cpp(const DataFrame pt, const DataFrame segments) {
  
  // Define a value type that stores a segment, its line ID, and its segment ID
  typedef std::tuple<SEGMENT3D, int, int> value; // SEGMENT3D, line ID, segment ID
  
  // Extract coordinates and line IDs from segment
  NumericVector x1 = segments["from_x"];
  NumericVector y1 = segments["from_y"];
  NumericVector z1 = segments["from_z"];
  NumericVector x2 = segments["to_x"];
  NumericVector y2 = segments["to_y"];
  NumericVector z2 = segments["to_z"];
  IntegerVector L1 = segments["L1"];
  IntegerVector L2 = segments["L2"];
  int n_segments = x1.size();
  
  // Build the r-tree with segments
  bgi::rtree<value, bgi::quadratic<16>> rtree;
  for (int i = 0; i < n_segments; ++i) {
    POINT3D p1(x1[i], y1[i], z1[i]);
    POINT3D p2(x2[i], y2[i], z2[i]);
    SEGMENT3D segment(p1, p2);
    rtree.insert(std::make_tuple(segment, L1[i], L2[i])); // Store L1 and SegmentID
  }
  
  // Extract coordinates from query points
  NumericVector px = pt["X"];
  NumericVector py = pt["Y"];
  NumericVector pz = pt["Z"];
  int n_points = px.size();
  
  // Prepare vector for the result
  IntegerVector result(n_points);
  
  // Perform nearest segment search for each query point
  for (int i = 0; i < n_points; ++i) {
    POINT3D FromPoint(px[i], py[i], pz[i]);
    std::vector<value> neighbours;
    
    // Query the r-tree for the nearest neighbour (segment)
    rtree.query(bgi::nearest(FromPoint, 1), std::back_inserter(neighbours));
    
    // Extract the nearest segment's line ID (L1)
    int nearest_L1 = std::get<1>(neighbours[0]);
    
    // Store results
    result[i] = nearest_L1; // Line ID of the nearest segment
  }
  
  // Return results as a DataFrame
  return result;
}





