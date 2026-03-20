// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;
#include <boost/geometry.hpp>
namespace bg = boost::geometry;

// Define data type: <POINT3D> and <LINESTRING3D>
typedef bg::model::point<double, 3, bg::cs::cartesian> POINT3D;
typedef bg::model::linestring<POINT3D> LINESTRING3D;

// This function computes the nearest point for each pair of input point and linestring
// [[Rcpp::export]]
DataFrame closestPT3d_cpp(DataFrame point_df, DataFrame lines_df) {
  // Extract columns from the point DataFrame
  NumericVector PTx = point_df["X"];
  NumericVector PTy = point_df["Y"];
  NumericVector PTz = point_df["Z"];
  
  // Extract columns from the lines DataFrame
  NumericVector x = lines_df["X"];
  NumericVector y = lines_df["Y"];
  NumericVector z = lines_df["Z"];
  IntegerVector id = lines_df["L1"];
  
  // Create a map to store linestrings by ID
  std::map<int, LINESTRING3D> linestrings;
  
  // Populate the linestrings map with points grouped by line ID
  for (size_t i = 0; i < id.size(); i++) {
    int line_id = id[i];
    double x_coord = x[i];
    double y_coord = y[i];
    double z_coord = z[i];
    
    // Add points to the corresponding linestring
    linestrings[line_id].push_back(POINT3D(x_coord, y_coord, z_coord));
  }
  
  // Define containers to store the output
  std::vector<double> nearest_X;
  std::vector<double> nearest_Y;
  std::vector<double> nearest_Z;
  
  // Perform pairwise computation
  for (size_t i = 0; i < PTx.size(); i++) {
    POINT3D input_point(PTx[i], PTy[i], PTz[i]); // Current input point
    LINESTRING3D linestring = linestrings[i+1]; // Corresponding linestring (ID matches row index + 1)
    
    // Use Boost.Geometry to calculate the closest segment linking the point and the linestring
    bg::model::segment<POINT3D> closest_seg;
    bg::closest_points(input_point, linestring, closest_seg);
    
    // Extract the closest point from the segment
    POINT3D closest_point = closest_seg.second;
    
    // Store the results
    nearest_X.push_back(bg::get<0>(closest_point));
    nearest_Y.push_back(bg::get<1>(closest_point));
    nearest_Z.push_back(bg::get<2>(closest_point));
  }
  
  // Return results as a DataFrame
  return DataFrame::create(
    Named("X") = nearest_X,
    Named("Y") = nearest_Y,
    Named("Z") = nearest_Z
  );
}

