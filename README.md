# **`GISnetwork3D` − An open-sourced 3D GIS and 3D spatial access measurement tool in R and QGIS**

A package for 3D GIS operation and 3D network analysis.

## Brief Introduction

Traditional map-making and GIS software have long adhered to a planar conception of space. Contrary to proprietary GIS software, open-source alternatives like QGIS and R currently lack comparable capabilities to perform even simple 3D GIS operations and 3D network analysis. Therefore, to bridge this divide, we developed an innovative R package that empowers users to perform 3D GIS operations and advanced 3D network analyses, bringing the power of 3D modelling into open-source environments. Additionally, by linking QGIS with this R package, we seek to create a user-friendly, code-free platform that democratises access to advanced 3D GIS functionalities, fostering new insights and applications across disciplines.

The `GISnetwork3D` has the following objectives:

1.  **Supports 3D geometry operations:** Provides functionality for measuring the geometry of 3D objects and querying spatial and topological relationships between 3D geometries. Examples include 3D length measurement, K-nearest neighbour analysis, distance calculation, and snapping functions.

2.  **Supports automated/convenient 3D path cost calculation:** Supports calculating the path cost for each segment (pre-calculated walking time and energy consumption), and aggregates the results without requiring segmentation of the entire network.

3.  **Supports slope and direction-aware 3D routing:** This routing tool considers the following factors simultaneously: (1) path direction; (2) a single linestring with dual weight by directions; (3) a network containing both one-way and two-way paths (some paths have dual weights, while others have only one).

4.  **Supports user-friendly, comprehensive 3D accessibility analysis:** Provides user-friendly functionality for: (1) single-purpose and multi-purpose accessibility analysis, (2) path attribute summary and (3) route geometry generation in one single code/function.

# Cheatsheets

## 1. Geometry Operation

![Cheatsheet for Geometry Operation](images/Figure. GIS Operation.png)

## 2. Network Construction

![Cheatsheet for Network Construction](images/Figure. Network Construction.png)

## 3. Network Analysis

![Cheatsheet for Network Construction](images/Figure. Network Analysis.png)

# Installing R Package

The R package is hosted on GitHub. The installation requires the `devtools` package for this process. Please run the following code only if you do not have `devtools` installed.

```{install.packages("devtools")}
```

Then, you can install the `GISnetwork3D` package directly from the GitHub. This package relies on the following R packages: `Rcpp`, `sf`, `dplyr`, `igraph`, `raster`, `sfheaders`, `purrr`, `furrr`, `future` and `BH`. Their availability will be detected and be installed, in case missing, automatically during the installation of the `GISnetwork3D` package.

```{library(devtools)}
install_github("NKY-B/GISnetwork3D", auth_token = "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")
```
