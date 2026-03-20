// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;
#include <boost/geometry.hpp>
namespace bg = boost::geometry;

// Define data type: <POINT3D> and <LINESTRING3D>
typedef bg::model::point<double, 3, bg::cs::cartesian> POINT3D;
typedef bg::model::linestring<POINT3D> LINESTRING3D;

//_______________________________________________________________________________________________________________________________
// ---> dist3d_pt2pt_cpp - Point to Point Calculation
// [[Rcpp::export]]
NumericMatrix dist3d_pt2pt_cpp(NumericMatrix points1, NumericMatrix points2) {
  int n1 = points1.nrow();
  int n2 = points2.nrow();
  NumericMatrix distances(n1, n2);
  
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      POINT3D p1(points1(i, 0), points1(i, 1), points1(i, 2));
      POINT3D p2(points2(j, 0), points2(j, 1), points2(j, 2));
      distances(i, j) = bg::distance(p1, p2);
    }
  }
  
  return distances;
}

//_______________________________________________________________________________________________________________________________
// ---> dist3d_pt2line_cpp - Point to Line Calculation
// [[Rcpp::export]]
NumericMatrix dist3d_pt2line_cpp(NumericMatrix points, DataFrame lines_df) {
  int n_points = points.nrow();
  
  // Extract columns from the lines DataFrame
  NumericVector x = lines_df["X"];
  NumericVector y = lines_df["Y"];
  NumericVector z = lines_df["Z"];
  IntegerVector id = lines_df["L1"];
  
  // Create a map to store linestrings by ID
  std::map<int, LINESTRING3D> linestrings;
  
  // Loop through the points and group them by ID
  for (size_t i = 0; i < id.size(); i++) {
    int line_id = id[i];
    double x_coord = x[i];
    double y_coord = y[i];
    double z_coord = z[i];
    
    // Add points to the corresponding linestring
    linestrings[line_id].push_back(POINT3D(x_coord, y_coord, z_coord));
  }
  
  // Prepare the distance matrix
  NumericMatrix distances(n_points, linestrings.size());
  int line_index = 0;
  
  // Calculate distances from each point to each linestring
  for (const auto& line_pair : linestrings) {
    const auto& linestring = line_pair.second;
    
    for (int j = 0; j < n_points; ++j) {
      POINT3D point(points(j, 0), points(j, 1), points(j, 2));
      distances(j, line_index) = bg::distance(point, linestring);
    }
    line_index++;
  }
  
  return distances;
}

//_______________________________________________________________________________________________________________________________
// ---> dist3d_pt2pt_pair_cpp - Point to Point Calculation (Pairwise)
// [[Rcpp::export]]
NumericVector dist3d_pt2pt_pair_cpp(NumericMatrix points1, NumericMatrix points2) {
  int n = points1.nrow();
  NumericVector distances(n);
  
  for (int i = 0; i < n; ++i) {
    POINT3D p1(points1(i, 0), points1(i, 1), points1(i, 2));
    POINT3D p2(points2(i, 0), points2(i, 1), points2(i, 2));
    distances[i] = bg::distance(p1, p2);
  }
  
  return distances;
}


//_______________________________________________________________________________________________________________________________
// ---> dist3d_pt2line_pair_cpp - Point to Line Calculation (Pairwise)

// [[Rcpp::export]]
NumericVector dist3d_pt2line_pair_cpp(NumericMatrix points, DataFrame lines_df) {
  
  // Extract columns from the lines DataFrame
  NumericVector x = lines_df["X"];
  NumericVector y = lines_df["Y"];
  NumericVector z = lines_df["Z"];
  IntegerVector id = lines_df["L1"];
  
  // Create a map to store linestrings by ID
  std::map<int, LINESTRING3D> linestrings;
  
  // Loop through the points and group them by ID
  for (size_t i = 0; i < id.size(); i++) {
    int line_id = id[i];
    double x_coord = x[i];
    double y_coord = y[i];
    double z_coord = z[i];
    
    // Add points to the corresponding linestring
    linestrings[line_id].push_back(POINT3D(x_coord, y_coord, z_coord));
  }
  
  // Prepare the distance matrix
  int n = points.nrow();
  NumericVector distances(n);
  int line_index = 0;
  
  // Calculate distances from each point to each linestring
  for (const auto& line_pair : linestrings) {
    const auto& linestring = line_pair.second;
    POINT3D point(points(line_index, 0), points(line_index, 1), points(line_index, 2));
    distances(line_index) = bg::distance(point, linestring);
    line_index++;
  }
  
  return distances;
}
