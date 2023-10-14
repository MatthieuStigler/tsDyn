library(tsDyn)
suppressMessages(library(dplyr))
library(purrr)
library(tidyr)
library(tibble)
suppressWarnings(RNGversion("3.5.3"))

data(IIPUs)

### grid ###
grid <-  crossing(include = c( "const", "trend","none", "both"),
                  lag = 1:2L, 
                  thDelay = 0:1L, 
                  test = c("1vs", "2vs3")) %>% 
  filter(thDelay<lag)


### run ###
res <- grid %>% 
  mutate(dat = pmap(list(include, lag, thDelay, test),
                    ~setarTest(IIPUs[1:100], 
                               include= ..1,
                               m = ..2,
                               thDelay = ..3,
                               test = ..4,
                               nboot = 2,
                               seed = 4323)))


### show F tests ###
res %>% 
  mutate(Ftests = map(dat, ~enframe(.$Ftests))) %>% 
  select(-dat) %>% 
  unnest(Ftests)%>% 
  as.data.frame() %>% 
  head(10) %>% 
  print(digits=3) 

res %>% 
  mutate(Ftests = map(dat, ~enframe(.$SSRs))) %>% 
  select(-dat) %>% 
  unnest(Ftests)%>% 
  as.data.frame() %>% 
  head(10) %>% 
  print(digits=3) 

