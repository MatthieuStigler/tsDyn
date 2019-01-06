library(tidyverse)

IP <- read_table(file="data-raw/ipdat.txt", col_names = FALSE) %>% 
  rename(ind_prod=X1)

IP_c <- IP %>% 
  mutate(ind_prod = log(ind_prod),
         lag = dplyr::lag(ind_prod, 12),
         diff = 100 * (ind_prod - lag),
         n_row = 1:n()) %>% 
  filter(n_row>=157+12)


IP_c
qplot(x=n_row, y = diff, geom = "line", data =IP_c)

#transform as in Hansen 1999
dat<-log(IP$ind_prod)
dat2<-diff(dat, 12)*100 # dat=(dat(13:length(dat(:,1)))-dat(1:length(dat(:,1))-12))*100
dat3<-dat2[157:length(dat2)] #dat=dat(157:length(dat(:,1)));
dat4<-ts(dat3, start=c(1960,1), freq=12)

IP_c
head(dat4)

all.equal(IP_c$diff, as.numeric(dat4))

### 
library(usethis)
IIPUs <- dat4
use_data(IIPUs, overwrite = TRUE)
