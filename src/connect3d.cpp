#include <Rcpp.h>
#include <vector>
#include <iostream>
using namespace Rcpp;

//_______________________________________________________________________________________________________________________________
// ---> connect3d_pair_cpp - Create data.frame of segments that connect each pair of points

// [[Rcpp::export]] 
DataFrame connect3d_pair_cpp(const DataFrame pt1, const DataFrame pt2) {
  
  // Extract coordinates from pt1
  NumericVector x1 = pt1["X"];
  NumericVector y1 = pt1["Y"];
  NumericVector z1 = pt1["Z"];
  
  // Extract coordinates from pt2
  NumericVector x2 = pt2["X"];
  NumericVector y2 = pt2["Y"];
  NumericVector z2 = pt2["Z"];
  
  // Define containers to store the output
  std::vector<double> X;
  std::vector<double> Y;
  std::vector<double> Z;
  std::vector<int> id;
  
  // Loop through the points and create line segment
  for (size_t i = 0; i < x1.size(); i++) {
    X.push_back(x1[i]);
    X.push_back(x2[i]);
    Y.push_back(y1[i]);
    Y.push_back(y2[i]);
    Z.push_back(z1[i]);
    Z.push_back(z2[i]);
    id.push_back(i+1);
    id.push_back(i+1);
    }
  
  // Return results as a DataFrame
  return DataFrame::create(
    Named("X") = X,
    Named("Y") = Y,
    Named("Z") = Z,
    Named("id") = id);
}


//_______________________________________________________________________________________________________________________________
// ---> connect3d_ManyToMany_cpp - Create data.frame of segments that connect all input points

// [[Rcpp::export]] 
DataFrame connect3d_ManyToMany_cpp(const DataFrame pt1, const DataFrame pt2) {
  
  // Extract coordinates from pt1
  NumericVector x1 = pt1["X"];
  NumericVector y1 = pt1["Y"];
  NumericVector z1 = pt1["Z"];
  
  // Extract coordinates from pt2
  NumericVector x2 = pt2["X"];
  NumericVector y2 = pt2["Y"];
  NumericVector z2 = pt2["Z"];
  
  // Define containers to store the output
  std::vector<double> X;
  std::vector<double> Y;
  std::vector<double> Z;
  std::vector<int> id1;
  std::vector<int> id2;
  
  // Loop through the points and create line segment
  for (size_t i = 0; i < x1.size(); i++) {
    for (size_t j = 0; j < x2.size(); j++){
      X.push_back(x1[i]);
      X.push_back(x2[j]);
      Y.push_back(y1[i]);
      Y.push_back(y2[j]);
      Z.push_back(z1[i]);
      Z.push_back(z2[j]);
      
      id1.push_back(i+1);
      id1.push_back(i+1);
      id2.push_back(j+1);
      id2.push_back(j+1);
      }
    }
  
  // Return results as a DataFrame
  return DataFrame::create(
    Named("X") = X,
    Named("Y") = Y,
    Named("Z") = Z,
    Named("id1") = id1,
    Named("id2") = id2);
}

