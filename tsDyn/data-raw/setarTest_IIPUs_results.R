library(tsDyn)
data(IIPUs)
B <-  1000
test_1 <- setarTest(IIPUs, m=16, thDelay=5, nboot=B)
test_2 <- setarTest(IIPUs, m=16, thDelay=5, nboot=B, test = "2vs3")

test_12 <-  list(test_1 = test_1, test_2 = test_2)

library(usethis)
setarTest_IIPUs_results <- test_12
use_data(setarTest_IIPUs_results, overwrite = TRUE)
