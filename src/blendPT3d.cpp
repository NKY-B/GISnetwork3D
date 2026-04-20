// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;
#include <boost/geometry.hpp>
namespace bg = boost::geometry;

// Define data type: <POINT3D> and <LINESTRING3D>
typedef bg::model::point<double, 3, bg::cs::cartesian> POINT3D;
typedef bg::model::linestring<POINT3D> LINESTRING3D;

//_______________________________________________________________________________________________________________________________
// ---> blendPT3d_cpp - Blend point to line
// [[Rcpp::export]]
DataFrame blendPT3d_cpp(const DataFrame nearestPt, const DataFrame firstNode, const DataFrame lastNode, IntegerVector segID) {
  // ---> Define parameters
  // Point to blend
  NumericVector ptx = nearestPt["X"];
  NumericVector pty = nearestPt["Y"];
  NumericVector ptz = nearestPt["Z"];
  IntegerVector pt_segID = nearestPt["segID"];
  
  // From node
  NumericVector fromx = firstNode["X"];
  NumericVector fromy = firstNode["Y"];
  NumericVector fromz = firstNode["Z"];
  
  // Last node
  NumericVector Tox = lastNode["X"];
  NumericVector Toy = lastNode["Y"];
  NumericVector Toz = lastNode["Z"];
  
  // Define containers to store the output
  std::vector<double> X, Y, Z, dist;
  std::vector<int> id;
  
  for (int i = 0; i < segID.size(); ++i) {
    // Add first node
    X.push_back(fromx[i]);
    Y.push_back(fromy[i]);
    Z.push_back(fromz[i]);
    id.push_back(segID[i]);
    dist.push_back(0);
    
    for (int j = 0; j < pt_segID.size(); ++j) {
      if (segID[i] == pt_segID[j] &&
          (fromx[i] != ptx[j] || fromy[i] != pty[j] || fromz[i] != ptz[j]) &&
          (Tox[i] != ptx[j] || Toy[i] != pty[j] || Toz[i] != ptz[j])) {
        
        X.push_back(ptx[j]);
        Y.push_back(pty[j]);
        Z.push_back(ptz[j]);
        id.push_back(pt_segID[j]);
        
        // Calculate distance
        POINT3D p1(fromx[i], fromy[i], fromz[i]);
        POINT3D p2(ptx[j], pty[j], ptz[j]);
        dist.push_back(bg::distance(p1, p2));
      }
      }
    // Add last node
    X.push_back(Tox[i]);
    Y.push_back(Toy[i]);
    Z.push_back(Toz[i]);
    id.push_back(segID[i]);
    
    // Calculate distance
    POINT3D p1(fromx[i], fromy[i], fromz[i]);
    POINT3D p2(Tox[i], Toy[i], Toz[i]);
    dist.push_back(bg::distance(p1, p2));
    }
  
  // Return results as a DataFrame
  return DataFrame::create(
    Named("X") = X,
    Named("Y") = Y,
    Named("Z") = Z,
    Named("dist") = dist,
    Named("sid") = id);
  }
