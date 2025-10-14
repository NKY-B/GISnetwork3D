// [[Rcpp::depends(BH)]]
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
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

//_______________________________________________________________________________________________________________________________
// ---> knn3d_p2p_cpp - Return the K-nearest neighbouring points for each point

// [[Rcpp::export]] 
DataFrame knn3d_p2p_cpp(const DataFrame pt1, const DataFrame pt2, const int k) {
  
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
  // Prepare the containers to store result
  std::vector<int> id1; // ID of the from point
  std::vector<int> id2; // ID of the to point
  std::vector<double> distances; // Extract the distance
  
  // Perform KNN search for each point in pt1
  for (int i = 0; i < n1; ++i) {
    POINT3D query_point(x1[i], y1[i], z1[i]);
    std::vector<value> neighbours;
    
    // Query the r-tree for k nearest neighbours
    rtree.query(bgi::nearest(query_point, k), std::back_inserter(neighbours));
    
    for (int j = 0; j < neighbours.size(); ++j){
      double distance = bg::distance(query_point, std::get<0>(neighbours[j]));
      id1.push_back(i+1);
      id2.push_back(std::get<1>(neighbours[j]));
      distances.push_back(distance);
      }
  }
  // Return the result as a DataFrame
  return DataFrame::create(
    Named("id1") = id1,
    Named("id2") = id2,
    Named("dist") = distances);
}

