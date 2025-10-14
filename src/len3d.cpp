// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;
#include <boost/geometry.hpp>
namespace bg = boost::geometry;

// Define data type: <POINT3D> and <LINESTRING3D>
typedef bg::model::point<double, 3, bg::cs::cartesian> POINT3D;
typedef bg::model::linestring<POINT3D> LINESTRING3D;

//_______________________________________________________________________________________________________________________________
// ---> len3d_cpp - Calculate the 3D length for linestring
// [[Rcpp::export]]
NumericVector len3d_cpp(DataFrame lines_df) {
  
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
    
    // Add point to each linestring by looping through ID in the map object
    linestrings[line_id].push_back(POINT3D(x_coord, y_coord, z_coord));
  }
  
  // Prepare the len output
  NumericVector len(linestrings.size());
  int line_index = 0;
  
  // Calculate length
  for (const auto& line_pair : linestrings) {
    const auto& linestring = line_pair.second;
    len(line_index) = bg::length(linestring);
    line_index++;
  }
  
  return len;
}

