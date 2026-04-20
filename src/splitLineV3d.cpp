#include <Rcpp.h>
using namespace Rcpp;

//_______________________________________________________________________________________________________________________________
// ---> splitLineV3d_cpp - Split line by vertices

// [[Rcpp::export]]
DataFrame splitLineV3d_cpp(NumericVector x, 
                           NumericVector y,
                           NumericVector z,
                           IntegerVector L1,
                           IntegerVector split) {
  
  // Create containers for the for loop
  std::vector<double> X; // Extract the X coordinate of each vertices from the subseted line
  std::vector<double> Y; // Extract the Y coordinate of each vertices from the subseted line
  std::vector<double> Z; // Extract the Z coordinate of each vertices from the subseted line
  std::vector<int> newL1; // For storing the new linestring index of each vertices from the subseted line
  std::vector<int> oL1; // For storing the old linestring index of each vertices from the subseted line
  
  // Initialize variables
  int currentNewL1 = 1; // Start new linestring index
  int currentOriginalL1 = L1[0]; // Track original linestring index
  
  X.push_back(x[0]);
  Y.push_back(y[0]);
  Z.push_back(z[0]);
  oL1.push_back(L1[0]);
  newL1.push_back(currentNewL1);
  
  for (size_t i = 1; i < x.size() - 1; ++i) {
    if (L1[i - 1] == L1[i] && L1[i] == L1[i + 1] && split[i] == 1) {
      // If the current vertex is a split point, and the vertices before and after belong to the same linestring
      X.push_back(x[i]);
      Y.push_back(y[i]);
      Z.push_back(z[i]);
      oL1.push_back(L1[i]);
      newL1.push_back(currentNewL1);
      
      currentNewL1++; // For the new added rows, the first vertex belongs to the previous linestring, and the second vertex to next linestring
      
      X.push_back(x[i]);
      Y.push_back(y[i]);
      Z.push_back(z[i]);
      oL1.push_back(L1[i]);
      newL1.push_back(currentNewL1);
    } else {
      if (L1[i] != L1[i - 1]) {
        // If the current L1 is greater than the previous one, start a new linestring
        currentNewL1++;
      }
      X.push_back(x[i]);
      Y.push_back(y[i]);
      Z.push_back(z[i]);
      oL1.push_back(L1[i]);
      newL1.push_back(currentNewL1);
    }
  }
  
  // Add the last vertex
  X.push_back(x[x.size() - 1]);
  Y.push_back(y[y.size() - 1]);
  Z.push_back(z[z.size() - 1]);
  oL1.push_back(L1[L1.size() - 1]);
  newL1.push_back(currentNewL1);
  
  // Return the result as a DataFrame
  return DataFrame::create(
    Named("L1") = oL1,
    Named("NewL1") = newL1,
    Named("X") = X,
    Named("Y") = Y,
    Named("Z") = Z
  );
}
