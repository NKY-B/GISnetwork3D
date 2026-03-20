// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <vector>
#include <iostream>
#include <boost/geometry.hpp>
#include <cmath>
using namespace Rcpp;
namespace bg = boost::geometry;

// Define data type: <POINT2D>
typedef bg::model::point<double, 2, bg::cs::cartesian> POINT2D;
typedef bg::model::point<double, 3, bg::cs::cartesian> POINT3D;


// Function for tobr_hf.pace
double tobr_hf_pace(double gradientSlope) {
  double p = 0.6 * exp(3.5 * fabs(gradientSlope + 0.05));
  return p;
}

// Function for campbell.pace
double campbell_pace(double degreeSlope, double a, double b, double c, double d, double e) {
  double r = c * (1 / (M_PI * b * (1 + pow((degreeSlope - a) / b, 2)))) + d + e * degreeSlope;
  return 1 / r;
}

// Function for energy expenditure
double EE(double speed, double percentSlope) {
  // Calculate energy expenditure in W/kg
  double result = 1.44 + 1.94 * pow(speed, 0.43) + 0.24 * pow(speed, 4) +
    0.34 * speed * percentSlope * (1 - pow(1.05, 1 - pow(1.1, percentSlope + 32)));
  // Convert from W/kg to kcal/kg/hr
  result *= 0.860421;
  // Convert from kcal/hr to kcal/kg/s
  result /= 3600;
  return result;
}

// Function for calculating slopes in gradient, percent and degree
double calslope(double hdist, double vdist, const std::string& output_type) {
  // Calculate the slope values
  double percent_slope = (vdist / hdist) * 100;
  double degree_slope = atan(vdist / hdist) * (180 / M_PI);
  double gradient_slope = vdist / hdist;

  // Return the requested output
  if (output_type == "percent") {
    return percent_slope;
  } else if (output_type == "degrees") {
    return degree_slope;
  } else if (output_type == "gradient") {
    return gradient_slope;
  } else {
    std::cerr << "Invalid output type. Please choose 'percent', 'degrees', or 'gradient'." << std::endl;
    return NAN; // Return NaN for invalid input
  }
}

//_______________________________________________________________________________________________________________________________
// ---> calCost_net3d_THF_cpp - Calculate the travel cost for each linestring
// This function uses the Tobler's hiking function to calculate the travel time and use the new LCDA energy expenditure function proposed by Looney et al (2019) to calculate the travel energy.

// [[Rcpp::export]]
DataFrame calCost_net3d_THF_cpp(const DataFrame lines_df, const IntegerVector uniqueL1, double pSlope) {

  // Extract coordinates and index from lines_df
  IntegerVector L1 = lines_df["L1"]; // Line identifiers
  NumericVector x1 = lines_df["from_x"]; // X-coordinates of the start points
  NumericVector y1 = lines_df["from_y"]; // Y-coordinates of the start points
  NumericVector z1 = lines_df["from_z"]; // Z-coordinates of the start points
  NumericVector x2 = lines_df["to_x"]; // X-coordinates of the end points
  NumericVector y2 = lines_df["to_y"]; // Y-coordinates of the end points
  NumericVector z2 = lines_df["to_z"]; // Z-coordinates of the end points

  // Define containers to store the output values by segment
  std::vector<double> timeList; // To store calculated travel times
  std::vector<double> energyList; // To store calculated energy expenditures
  std::vector<double> len3dList; // To store 3D lengths
  std::vector<double> slopeList; // To store slope indicators

  // Loop through the lines to calculate the travel costs
  for (size_t i = 0; i < x1.size(); i++) {

    // Calculate horizontal distance between points
    POINT2D fromPT2D(x1[i], y1[i]); // Start point in 2D
    POINT2D toPT2D(x2[i], y2[i]); // End point in 2D
    double hdist = bg::distance(fromPT2D, toPT2D); // Horizontal distance

    // Calculate surface distance between points in 3D
    POINT3D fromPT3D(x1[i], y1[i], z1[i]); // Start point in 3D
    POINT3D toPT3D(x2[i], y2[i], z2[i]); // End point in 3D
    double sdist = bg::distance(fromPT3D, toPT3D); // Surface distance
    len3dList.push_back(sdist); // Store 3D length

    // Calculate vertical distance
    double vdist = z2[i] - z1[i]; // Vertical distance

    // Calculate travel time using Tobler's hiking function
    double pace = tobr_hf_pace(calslope(hdist, vdist, "gradient")); // Pace based on slope
    double travelTime = pace * sdist; // Total travel time
    if(hdist == 0){travelTime = 0;} // Set travel time to 0 if horizontal distance is 0
    timeList.push_back(travelTime); // Store calculated travel time

    // Calculate travel energy expenditure
    double energyExpenditure = EE(1/pace, calslope(hdist, vdist, "percent")); // Energy expenditure based on pace
    if(hdist == 0){
      energyList.push_back(0); // Store 0 energy if horizontal distance is 0
    } else {
      energyList.push_back(travelTime * energyExpenditure); // Store calculated energy expenditure
    }

    // Calculate slope
    double slope = calslope(hdist, vdist, "percent"); // Calculate slope percentage
    double slopeTRUE = 0; // Initialise slope indicator
    if(std::abs(slope) >= std::abs(pSlope)){ // Check if slope exceeds the threshold
      slopeTRUE = 1; // Set indicator to 1 if slope exceeds threshold
    }
    slopeList.push_back(slopeTRUE); // Store slope indicator
  }

  // Sum costs by L1
  // Create containers to store the aggregated results
  std::vector<double> time; // Total travel time for each line
  std::vector<double> energy; // Total energy expenditure for each line
  std::vector<double> len3d; // Total 3D length for each line
  std::vector<double> slope_len3d; // Total 3D length for slopes
  std::vector<double> slope_time; // Total travel time for slopes
  std::vector<double> slope_energy; // Total energy for slopes

  // Initialise containers for unique line identifiers
  for (size_t i = 0; i < uniqueL1.size(); i++) {
    time.push_back(0); // Initialise travel time for each unique line
    energy.push_back(0); // Initialise energy for each unique line
    len3d.push_back(0); // Initialise 3D length for each unique line
    slope_len3d.push_back(0); // Initialise 3D slope length for each unique line
    slope_time.push_back(0); // Initialise slope travel time for each unique line
    slope_energy.push_back(0); // Initialise slope energy for each unique line
  }

  // Aggregate results for each linestring
  for (size_t i = 0; i < x1.size(); i++) {
    int index = L1[i] - 1; // Get the index for the current line
    time[index] = time[index] + timeList[i]; // Accumulate travel time for the line
    energy[index] = energy[index] + energyList[i]; // Accumulate energy expenditure for the line
    len3d[index] = len3d[index] + len3dList[i]; // Accumulate 3D length for the line

    slope_time[index] = slope_time[index] + timeList[i] * slopeList[i]; // Count time for traveling slopes
    slope_energy[index] = slope_energy[index] + energyList[i] * slopeList[i]; // Count energy for traveling slopes
    slope_len3d[index] = slope_len3d[index] + len3dList[i] * slopeList[i]; // Count length of traveling slopes
    }

  // Return results as a DataFrame
  return DataFrame::create(
    Named("L1") = uniqueL1, // Line identifiers
    Named("time") = time, // Total travel time
    Named("energy") = energy, // Total energy expenditure
    Named("len3d") = len3d, // Total 3D lengths
    Named("slope.time") = slope_time, // Total time for slopes
    Named("slope.energy") = slope_energy, // Total energy for slopes
    Named("slope.len3d") = slope_len3d // Total 3D lengths for slopes
  );
}


//_______________________________________________________________________________________________________________________________
// ---> calCost_net3d_CHF_cpp - Calculate the travel cost for each linestring
// This function uses the hiking function proposed by Campbell et al (2022) to calculate the travel time and use the new LCDA energy expenditure function proposed by Looney et al (2019) to calculate the travel energy.

// [[Rcpp::export]]
DataFrame calCost_net3d_CHF_cpp(const DataFrame lines_df, const IntegerVector uniqueL1, double pSlope,
                                double a, double b, double c, double d, double e) {

  // Extract coordinates and index from lines_df
  IntegerVector L1 = lines_df["L1"]; // Line identifiers
  NumericVector x1 = lines_df["from_x"]; // X-coordinates of the start points
  NumericVector y1 = lines_df["from_y"]; // Y-coordinates of the start points
  NumericVector z1 = lines_df["from_z"]; // Z-coordinates of the start points
  NumericVector x2 = lines_df["to_x"]; // X-coordinates of the end points
  NumericVector y2 = lines_df["to_y"]; // Y-coordinates of the end points
  NumericVector z2 = lines_df["to_z"]; // Z-coordinates of the end points

  // Define containers to store the output values by segment
  std::vector<double> timeList; // To store calculated travel times
  std::vector<double> energyList; // To store calculated energy expenditures
  std::vector<double> len3dList; // To store 3D lengths
  std::vector<double> slopeList; // To store slope indicators

  // Loop through the lines to calculate the travel costs
  for (size_t i = 0; i < x1.size(); i++) {

    // Calculate horizontal distance between points
    POINT2D fromPT2D(x1[i], y1[i]); // Start point in 2D
    POINT2D toPT2D(x2[i], y2[i]); // End point in 2D
    double hdist = bg::distance(fromPT2D, toPT2D); // Horizontal distance

    // Calculate surface distance between points in 3D
    POINT3D fromPT3D(x1[i], y1[i], z1[i]); // Start point in 3D
    POINT3D toPT3D(x2[i], y2[i], z2[i]); // End point in 3D
    double sdist = bg::distance(fromPT3D, toPT3D); // Surface distance
    len3dList.push_back(sdist); // Store 3D length

    // Calculate vertical distance
    double vdist = z2[i] - z1[i]; // Vertical distance

    // Calculate travel time using the Campbell's function
    double pace = campbell_pace(calslope(hdist, vdist, "degrees"), a, b, c, d, e); // Pace based on slope
    double travelTime = pace * sdist; // Total travel time
    if(hdist == 0){travelTime = 0;} // Set travel time to 0 if horizontal distance is 0
    timeList.push_back(travelTime); // Store calculated travel time

    // Calculate travel energy expenditure
    double energyExpenditure = EE(1/pace, calslope(hdist, vdist, "percent")); // Energy expenditure based on pace
    if(hdist == 0){
      energyList.push_back(0); // Store 0 energy if horizontal distance is 0
    } else {
      energyList.push_back(travelTime * energyExpenditure); // Store calculated energy expenditure
    }

    // Calculate slope
    double slope = calslope(hdist, vdist, "percent"); // Calculate slope percentage
    double slopeTRUE = 0; // Initialise slope indicator
    if(std::abs(slope) >= std::abs(pSlope)){ // Check if slope exceeds the threshold
      slopeTRUE = 1; // Set indicator to 1 if slope exceeds threshold
    }
    slopeList.push_back(slopeTRUE); // Store slope indicator
  }

  // Sum costs by L1
  // Create containers to store the aggregated results
  std::vector<double> time; // Total travel time for each line
  std::vector<double> energy; // Total energy expenditure for each line
  std::vector<double> len3d; // Total 3D length for each line
  std::vector<double> slope_len3d; // Total 3D length for slopes
  std::vector<double> slope_time; // Total travel time for slopes
  std::vector<double> slope_energy; // Total energy for slopes

  // Initialise containers for unique line identifiers
  for (size_t i = 0; i < uniqueL1.size(); i++) {
    time.push_back(0); // Initialise travel time for each unique line
    energy.push_back(0); // Initialise energy for each unique line
    len3d.push_back(0); // Initialise 3D length for each unique line
    slope_len3d.push_back(0); // Initialise 3D slope length for each unique line
    slope_time.push_back(0); // Initialise slope travel time for each unique line
    slope_energy.push_back(0); // Initialise slope energy for each unique line
  }

  // Aggregate results for each linestring
  for (size_t i = 0; i < x1.size(); i++) {
    int index = L1[i] - 1; // Get the index for the current line
    time[index] = time[index] + timeList[i]; // Accumulate travel time for the line
    energy[index] = energy[index] + energyList[i]; // Accumulate energy expenditure for the line
    len3d[index] = len3d[index] + len3dList[i]; // Accumulate 3D length for the line

    slope_time[index] = slope_time[index] + timeList[i] * slopeList[i]; // Count time for traveling slopes
    slope_energy[index] = slope_energy[index] + energyList[i] * slopeList[i]; // Count energy for traveling slopes
    slope_len3d[index] = slope_len3d[index] + len3dList[i] * slopeList[i]; // Count length of traveling slopes
  }

  // Return results as a DataFrame
  return DataFrame::create(
    Named("L1") = uniqueL1, // Line identifiers
    Named("time") = time, // Total travel time
    Named("energy") = energy, // Total energy expenditure
    Named("len3d") = len3d, // Total 3D lengths
    Named("slope.time") = slope_time, // Total time for slopes
    Named("slope.energy") = slope_energy, // Total energy for slopes
    Named("slope.len3d") = slope_len3d // Total 3D lengths for slopes
  );
}




//_______________________________________________________________________________________________________________________________
// ---> calSlope_net3d_cpp - Calculate the slope for each segment

// [[Rcpp::export]]
NumericVector calSlope_net3d_cpp(const DataFrame fromNode_df, const DataFrame toNode_df, const IntegerVector uniqueL1, const std::string& unit) {

  // Extract coordinates and Index from lines_df
  IntegerVector L1 = uniqueL1;
  NumericVector x1 = fromNode_df["X"];
  NumericVector y1 = fromNode_df["Y"];
  NumericVector z1 = fromNode_df["Z"];
  NumericVector x2 = toNode_df["X"];
  NumericVector y2 = toNode_df["Y"];
  NumericVector z2 = toNode_df["Z"];

  // Define containers to store the output
  NumericVector result(L1.size());

  // Loop through the lines to calculate the travel costs
  for (size_t i = 0; i < x1.size(); i++) {

    // Horizontal distance
    POINT2D fromPT2D(x1[i], y1[i]);
    POINT2D toPT2D(x2[i], y2[i]);
    double hdist = bg::distance(fromPT2D, toPT2D);

    // Vertical distance
    double vdist = z2[i] - z1[i];

    // Slope
    double slope = calslope(hdist, vdist, unit);
    result[i] = slope;
    }
  return result;
  }





