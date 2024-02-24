##GISnetwork3D=group
##line=vector line
##condition=string "NA"
##revLine=output vector

library(sf)
library(dplyr)
library(GISnetwork3D)

if(condition == "NA"){
condition = -999
} else {
condition = line[[condition]]
}

result = reverseLine_by_condition(line = line, condition = condition)

revLine = result


