#include <Rcpp.h>
#include <vector>
#include <iostream>
using namespace Rcpp;

//_______________________________________________________________________________________________________________________________
// ---> getFirstNode - subtract the first node of each linestring

// [[Rcpp::export]] 
NumericMatrix getFirstNode_cpp(const DataFrame line_df) {
  
  // Extract coordinates from line_df
  NumericVector x = line_df["X"];
  NumericVector y = line_df["Y"];
  NumericVector z = line_df["Z"];
  NumericVector L1 = line_df["L1"];
  
  // Create a map to store linestrings by ID
  std::map<int, std::vector<std::array<double, 3>>> lines;
  
  // Loop through the points and group them by L1
  for (size_t i = 0; i < L1.size(); i++) {
    int line_L1 = L1[i];
    double x_coord = x[i];
    double y_coord = y[i];
    double z_coord = z[i];
    
    // Group points by L1
    lines[line_L1].push_back({x_coord, y_coord, z_coord});}
  
  // Prepare the output
  std::vector<std::array<double, 3>> result;
  
  // Subset the first row
  for (const auto& line_pair : lines) {
    const auto& linestring = line_pair.second;
    result.push_back(linestring[0]);
  }
  // Convert the result to a NumericMatrix
  int num_lines = result.size();
  NumericMatrix output(num_lines, 3); // 3 columns for X, Y, Z
  
  for (int i = 0; i < num_lines; i++) {
    output(i, 0) = result[i][0]; // X
    output(i, 1) = result[i][1]; // Y
    output(i, 2) = result[i][2]; // Z
  }
  
  return output;
}



//_______________________________________________________________________________________________________________________________
// ---> getLastNode - subtract the last node of each linestring

// [[Rcpp::export]] 
NumericMatrix getLastNode_cpp(const DataFrame line_df) {
  
  // Extract coordinates from line_df
  NumericVector x = line_df["X"];
  NumericVector y = line_df["Y"];
  NumericVector z = line_df["Z"];
  NumericVector L1 = line_df["L1"];
  
  // Create a map to store linestrings by L1
  std::map<int, std::vector<std::array<double, 3>>> lines;
  
  // Loop through the points and group them by L1
  for (size_t i = 0; i < L1.size(); i++) {
    int line_L1 = L1[i];
    double x_coord = x[i];
    double y_coord = y[i];
    double z_coord = z[i];
    
    // Group points by L1
    lines[line_L1].push_back({x_coord, y_coord, z_coord});}
  
  // Prepare the output
  std::vector<std::array<double, 3>> result;
  
  // Subet the last row
  for (const auto& line_pair : lines) {
    const auto& linestring = line_pair.second;
    result.push_back(linestring.back());
  }
  // Convert the result to a NumericMatrix
  int num_lines = result.size();
  NumericMatrix output(num_lines, 3); // 3 columns for X, Y, Z
  
  for (int i = 0; i < num_lines; i++) {
    output(i, 0) = result[i][0]; // X
    output(i, 1) = result[i][1]; // Y
    output(i, 2) = result[i][2]; // Z
  }
  
  return output;
}



//_______________________________________________________________________________________________________________________________
// ---> getNodeIndex_cpp - For a dataframe of linestring, return the index of the first and last nodes

// [[Rcpp::export]] 
IntegerVector getNodeIndex_cpp(const DataFrame line_df) {
  
  // Extract L1 from line_df
  NumericVector L1 = line_df["L1"];
  
  // Create a map to store linestrings by ID
  std::map<int, std::vector<int>> lines;
  
  // Loop through the L1 and group them by L1
  for (size_t i = 0; i < L1.size(); i++) {
    int line_L1 = L1[i];
    
    // Group points by L1
    lines[line_L1].push_back(i);}
  
  // Prepare the output
  std::vector<int> result;
  
  // Subset the first row
  for (int i = 0; i < lines.size(); ++i) {
    std::vector<int> subline = lines[i+1]; // Extract the line by L1
    
    for (int j = 0; j < subline.size(); ++j){
      if(j == 0){ result.push_back(1); } // For the first node, return 1
      else if(j == (subline.size() - 1)){ result.push_back(2);} // For the last node, return 2
      else{ result.push_back(0); } // For the vertiecs, return 0
      }
  }
  
  return wrap(result);
}
