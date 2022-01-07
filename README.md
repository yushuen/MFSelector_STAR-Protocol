# A parallel computing version of MFSelector based on doSNOW package
The original version of MFSelector (Monotonic Feature Selector) can be download via its [official site](http://microarray.ym.edu.tw:8080/tools/module/MFSelector/index.jsp?mode=home). We made minimum modifications of the original R code to enable it can perform parallel computing based on doSNOW package. We propose this version as a solution for parallel MFSelector on Windows platform.

## Usage
Assume you have required inputs ready, here is an example for executing MFSelector in parallel on Windows with our solution.

```
# Load required package
library(foreach)
library(doSNOW)

# To ensure the original MFSelector has been detached
detach("package:MFSelector", unload=TRUE)

# Read R code of Monotonic Feature Selector
source("MFSelector_doSNOW.r")

# Detect computing cores in the current machine
total_cores <- detectCores()

# Half of total cores
half_cores <- total_cores/2

# Run MFSelector with half of total cores
mfselector(data, nsc, stageord = F, stagename = F, type = 1, nline = T, dline = T, pdf = 1:100, cmp = 0, permut = 100, svdenoise = 0.03, svdetimes = 4, cores = half_cores)

```

## Note
* For details about the usage of MFSelector, please check the manual provided at [the official site](http://microarray.ym.edu.tw:8080/tools/module/MFSelector/index.jsp?mode=support).
* Modifications we made are all tracked.
