This document explores and implements procedures outlined in Chapter 24 - "Multiple outcomes or time-points within a study" of Borenstein et al's (2009) book, "Introduction to Meta-Analysis".



```r
library(xtable)
```




# Function for variance of effect sizes from two outcomes


```r
variance_two_effects <- function(vy1, vy2, r) {
    # vy1 variance of effect size 1: iv to dv1 vy2 variance of effect size 2:
    # iv to dv2 r correlation between dv1 and dv2
    1/4 * (vy1 + vy2 + 2 * r * sqrt(vy1) * sqrt(vy2))
}

variance_two_effects(vy1 = 0.05, vy2 = 0.05, r = 0.5)
```

```
## [1] 0.0375
```




The above analysis shows the variance of an averaged effect where the variance of each individual effect is 0.05 and the correlation between the two dependent variables is 0.50. The variance is reduced but not by as much as would be the case if the variables were independent 

# Implement example in Table 24.3
The following analyses implement an analysis of the data provided in the first part of Table 24.3, and produce the output provided in the latter part of Table 24.3.



```r
# create data from Table 24.3
table_24_3 <- matrix(c(0.3, 0.05, 0.1, 0.05, 0.5, 0.2, 0.02, 0.1, 
    0.02, 0.6, 0.4, 0.05, 0.2, 0.05, 0.6, 0.2, 0.01, 0.1, 0.01, 0.4, 0.4, 0.06, 
    0.3, 0.06, 0.8), byrow = TRUE, ncol = 5)
table_24_3 <- data.frame(table_24_3)
names(table_24_3) <- c("y1", "vy1", "y2", "vy2", "r")

# calculate pooled variance
table_24_3$vybar <- apply(table_24_3, 1, function(X) variance_two_effects(X["vy1"], 
    X["vy2"], X["r"]))

# create average effect size per study
table_24_3$ybar <- (table_24_3$y1 + table_24_3$y2)/2
table_24_3$weight <- 1/table_24_3$vybar
table_24_3$es_weight <- table_24_3$weight * table_24_3$ybar

# average effect size
sum(table_24_3$es_weight)/sum(table_24_3$weight)
```

```
## [1] 0.1819
```

```r

# variance of average effect size
1/sum(table_24_3$weight)
```

```
## [1] 0.003629
```

```r

# standard error of average effect size
sqrt(1/sum(table_24_3$weight))
```

```
## [1] 0.06024
```






```r
# reproduce entirety of table 24.3
print(xtable(round(table_24_3, 3)), type = "html")
```

<!-- html table generated in R 2.15.1 by xtable 1.7-0 package -->
<!-- Fri Aug 17 17:27:38 2012 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> y1 </TH> <TH> vy1 </TH> <TH> y2 </TH> <TH> vy2 </TH> <TH> r </TH> <TH> vybar </TH> <TH> ybar </TH> <TH> weight </TH> <TH> es_weight </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD align="right"> 0.30 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.10 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.50 </TD> <TD align="right"> 0.04 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 26.67 </TD> <TD align="right"> 5.33 </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.02 </TD> <TD align="right"> 0.10 </TD> <TD align="right"> 0.02 </TD> <TD align="right"> 0.60 </TD> <TD align="right"> 0.02 </TD> <TD align="right"> 0.15 </TD> <TD align="right"> 62.50 </TD> <TD align="right"> 9.38 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD align="right"> 0.40 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.60 </TD> <TD align="right"> 0.04 </TD> <TD align="right"> 0.30 </TD> <TD align="right"> 25.00 </TD> <TD align="right"> 7.50 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD align="right"> 0.20 </TD> <TD align="right"> 0.01 </TD> <TD align="right"> 0.10 </TD> <TD align="right"> 0.01 </TD> <TD align="right"> 0.40 </TD> <TD align="right"> 0.01 </TD> <TD align="right"> 0.15 </TD> <TD align="right"> 142.86 </TD> <TD align="right"> 21.43 </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD align="right"> 0.40 </TD> <TD align="right"> 0.06 </TD> <TD align="right"> 0.30 </TD> <TD align="right"> 0.06 </TD> <TD align="right"> 0.80 </TD> <TD align="right"> 0.05 </TD> <TD align="right"> 0.35 </TD> <TD align="right"> 18.52 </TD> <TD align="right"> 6.48 </TD> </TR>
   </TABLE>



In summary, the procedure is as follows

* Calculate mean and variance of mean effect size for each study
* Then complete analysis as per normal for a fixed-effect or random-effects meta-analysis.

# Multiple outcomes
The following function extends the calculation of variances of multiple outcomes to three or more outcomes.



```r
matrix2indexed <- function(x) {
    cbind(expand.grid(i = seq(nrow(x)), j = seq(nrow(x))), r = as.vector(x))
}

variance_multiple_effects <- function(V, r) {
    # V vector of variances for effect sizes: IV-DV1, IV-DV2, IV-DV3 r matrix
    # of correlations between DVs return: variance of average effect size
    
    as.vector(r)
    x <- matrix(1:16, ncol = 4)
    r_long <- matrix2indexed(r)
    # remove diagonal values
    r_long <- r_long[r_long$i != r_long$j, ]
    # add corresponding variances
    r_long$vi <- V[r_long$i]
    r_long$vj <- V[r_long$j]
    
    part1 <- 1/length(V)^2
    part2 <- sum(V)
    part3 <- sum(apply(r_long, 1, function(X) X["r"] * sqrt(X["vj"]) * sqrt(X["vi"])))
    vybar <- part1 * (part2 + part3)
    vybar
}
```




## Replicate 2 outcome function


```r
r_matrix <- diag(2)
r_matrix[lower.tri(r_matrix)] <- c(0.5)
r_matrix[upper.tri(r_matrix)] <- r_matrix[lower.tri(r_matrix)]
c(multiple = variance_multiple_effects(c(0.05, 0.05), r_matrix), 
    two = variance_two_effects(0.05, 0.05, 0.5))
```

```
## multiple      two 
##   0.0375   0.0375 
```




The following shows that the multiple outcome function is returning the same value as the two outcome function.


## Multiple outcome example


```r
r_matrix <- diag(3)
r_matrix[lower.tri(r_matrix)] <- c(0.1, 0.2, 0.3)
r_matrix[upper.tri(r_matrix)] <- r_matrix[lower.tri(r_matrix)]
variance_multiple_effects(c(0.1, 0.2, 0.3), r_matrix)
```

```
## [1] 0.09384
```




Above is an example of getting the variance of a combined effect based on the variance of three effect sizes.
