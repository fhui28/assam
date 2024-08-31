rm(list = ls())
library(tidyverse)
library(sdmTMB)
library(mgcv)

set.seed(2) ## simulate some data... 
dat <- gamSim(1,
              n = 1000,
              dist = "normal",
              scale = 1)

gamfit <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3),
         family = gaussian(),
         data = dat,
         method = "ML")
summary(gamfit)


sdmfit <- sdmTMB(y ~ s(x0) + s(x1) + s(x2) + s(x3),
              family = gaussian(),
              data = dat,
              reml = FALSE,
              spatial = FALSE)
summary(sdmfit)


gam.vcomp(gamfit)
sdmfit$sd_report$par.fixed[6:10] %>% exp

coef(gamfit)
sdmfit$sd_report$par.fixed


par(mfrow = c(2,2))
plot(gamfit, rug = TRUE)

par(mfrow = c(2,2))
plot_smooth(sdmfit, select = 1)
plot_smooth(sdmfit, select = 2)
plot_smooth(sdmfit, select = 3)
plot_smooth(sdmfit, select = 4)
