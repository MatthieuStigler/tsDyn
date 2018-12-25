library(tidyverse)
library(usethis)


## local
setwd("/home/matifou/Dropbox/Documents/tsDyn/tsDyn")

## this file is downloaded from the website. Right headers, not data used by Hansen
path_orig <-  "data-raw/zeroyld1"

## this one is from hansen. 'Right' data, wrong headers
path_z <-  "data-raw/zeroyld.dat"

dat <- read_table(path_z, col_names = FALSE) 
dat_orig <-  read_table(path_orig, skip = 1, col_names = FALSE)

## maturities
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
dat_c <-  dat %>% 
  mutate(Year = as.integer(Year),
         Month = as.integer(Month),
         Date = as.Date(paste(Year, Month, "01", sep="-"))) %>% 
  select(Date, Year, Month, everything())


dat_fin_old <- dat_c[, c("M_120", "M_12")]
colnames(dat_fin_old) <-  c("short.run", "long.run")

dat_fin <-  dat_c[, c("M_12", "M_120")]
colnames(dat_fin) <-  c("short.run", "long.run")

#############################
## compare with tsDyn (old)
#############################

# library(tsDyn)
# data(zeroyld)
# head(zeroyld)
# all.equal(zeroyld, as.data.frame(dat_fin), check.attributes = TRUE)

#############################
## hansen code (old)
#############################
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

## full one

#############################
## prepare for export
#############################

## data 2 + time cols
dat_fin_more <- dat_c[, c("Date", "Year", "Month", "M_12", "M_120")]
colnames(dat_fin_more) <- c("Date", "Year", "Month", "short.run", "long.run")
dat_fin_more

zeroyldMeta <- as.data.frame(dat_fin_more)

## data all
zeroyldFull <-  dat_c %>% 
  rename_all(str_to_title) %>% 
  as.data.frame()

sapply(zeroyldFull, class)

#############################

## export
#############################


use_data(zeroyldMeta, overwrite = TRUE)
# write.csv(dat_fin_more, "/home/matifou/Dropbox/Documents/tsDyn/tsDyn/data/zeroyldMeta.csv")
# write.csv2(zeroyldMeta, "data/zeroyldMeta.csv",
#            row.names = FALSE)


## test?
# test_csv <- read.csv("data/zeroyldMeta.csv")
# test_csv2 <- read.csv2("data/zeroyldMeta.csv")
# sapply(test_csv2, class)