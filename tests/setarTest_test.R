library(tsDyn)
suppressMessages(library(tidyverse))
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
                    ~setarTest(IIPUs, 
                               include= ..1,
                               m = ..2,
                               thDelay = ..3,
                               test = ..4,
                               nboot = 2)))


### show F tests ###
res %>% 
  mutate(Ftests = map(dat, ~enframe(.$Ftests))) %>% 
  select(-dat) %>% 
  unnest(Ftests)

res %>% 
  mutate(Ftests = map(dat, ~enframe(.$SSRs))) %>% 
  select(-dat) %>% 
  unnest(Ftests)

