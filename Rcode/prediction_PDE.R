# 11/11/25 prediction uncertainty
# 11/16/25 prediction as sample from ZIBinom for 1000 MCMC samples

library(fields, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(RColorBrewer, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(mvtnorm, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
# library(gridExtra)
library(fBasics, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
# library(coda, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(dplyr, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
library(egg, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
# library(patchwork)
# library(ggplot2, lib.loc = "/home/yawenguan/jiaqichen/test100/library")
# library(cowplot)
# library(grid)

# setwd("/home/yawenguan/jiaqichen/Region88/paper/simulation/prediction uncertainty")
source("/home/yawenguan/jiaqichen/test100/Rcode/functions_v6.70.R")

N_sample <- 1000

sink(file = "test.txt")
cat("test: ")
N_sample
sink(file = NULL)
