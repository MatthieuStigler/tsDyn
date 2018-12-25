library(tidyverse)
library(usethis)

path_orig <-  "/home/matifou/pCloudDrive/Documents/Papers/MyPapers/Articles/Threshold Cointegration Handbook/example_applied/zeroyld1"

path_z <-  "/home/matifou/pCloudDrive/Documents/Papers/MyPapers/Articles/Threshold Cointegration Handbook/example_applied/joe_02r/zeroyld.dat"
a <- read_csv(path_z)
dat <- read_table(path_z, col_names = FALSE)
dat_orig <-  read_table(path_orig, skip = 1, col_names = FALSE)

## mats
dat_orig[1:7,]
mats <- c(t(as.matrix(dat_orig[1:7,])))
mats_c <-  mats[!is.na(mats)]
length(mats_c) == 7* 8

## maturity per Hansen:
rs_c <- c(seq(0,18,1),21,24,30, seq(36,(36+7*12),12))

## compare
mats_c2 <-round(12*mats_c)

rs_c
all(rs_c== mats_c2[1:length(rs_c)] )

colnames(dat) <-  c("Year",  "Month", "NMATS", "NOBS", "SE", "TAX", paste("M", mats_c2, sep="_"))
dat_fin_old <- dat[, c("M_120", "M_12")]
colnames(dat_fin_old) <-  c("short.run", "long.run")

dat_fin <-  dat[, c("M_12", "M_120")]
colnames(dat_fin) <-  c("short.run", "long.run")

## compare with tsDyn
library(tsDyn)
data(zeroyld)
head(zeroyld)
all.equal(zeroyld, as.data.frame(dat_fin), check.attributes = TRUE)

## hansen code
dat_2 <- dat[,7:62]
rs <- rbind(as.matrix(seq(0,18,1)),21,24,30,as.matrix(seq(36,(36+7*12),12)))


short <- 12
long <- 120
which(rs==short)
which(rs==long)
short_i <- which(rs==short)
long_i <- which(rs==long)
dat_new <- dat_2[, c(long_i,short_i)]   

dat_new

## prepare
dat_fin_more <- dat[, c("Year", "Month", "M_12", "M_120")]
colnames(dat_fin_more) <- c("Year", "Month", "short.run", "long.run")
dat_fin_more

head(dat_fin_more)
tail(dat_fin_more)


## save
setwd("/home/matifou/Dropbox/Documents/tsDyn/tsDyn/")
zeroyld_meta <- dat_fin_more
# use_data(zerolyd_meta)
# write.csv(dat_fin_more, "/home/matifou/Dropbox/Documents/tsDyn/tsDyn/data/zeroyld_meta.csv")
write.csv(dat_fin_more, "/home/matifou/Dropbox/Documents/tsDyn/tsDyn/data/zeroyldMeta.csv")
