#include <Rcpp.h>
#include <vector>
#include <iostream>
using namespace Rcpp;

//_______________________________________________________________________________________________________________________________
// ---> segmentise_cpp - Segmentise linestring

// [[Rcpp::export]] 
List segmentise_cpp(const DataFrame lines_df) {
  
  // Extract coordinates from line_df
  NumericVector x = lines_df["X"];
  NumericVector y = lines_df["Y"];
  NumericVector z = lines_df["Z"];
  NumericVector id = lines_df["L1"];
  
  
  // Define container to store outputs
  std::vector<double> X;
  std::vector<double> Y;
  std::vector<double> Z;
  std::vector<int> L1;
  std::vector<int> L2;
  
  // Create recorder to index the L2
  int n = 0;
  
  // Loop through the points and group them by ID
  for (size_t i = 1; i < id.size(); i++) {
    
    if(id[i] == id[i-1]){
      n++;
      
      X.push_back(x[i-1]);
      X.push_back(x[i]);
      
      Y.push_back(y[i-1]);
      Y.push_back(y[i]);
      
      Z.push_back(z[i-1]);
      Z.push_back(z[i]);
      
      L1.push_back(id[i-1]);
      L1.push_back(id[i]);
      
      L2.push_back(n);
      L2.push_back(n);
    }
    }
    
    // Convert the containers to Rcpp vectors and return as a list
    return List::create(
      _["X"] = X,
      _["Y"] = Y,
      _["Z"] = Z,
      _["L1"] = L1,
      _["L2"] = L2
    );
  }
