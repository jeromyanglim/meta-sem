# Part 1 Prevost et al 2007 Analyses


The following implements some analysis from the following article

> Prevost, A.T., Mason, D., Griffin, S., Kinmonth, A.L., Sutton, S. & Spiegelhalter, D. (2007). Allowing for correlations between correlations in random-effects meta-analysis of correlation matrices.. Psychological methods, 12, 434. [PDF](http://people.cehd.tamu.edu/~vwillson/Videos/E642%20meta/Covarying%20rs%20in%20MA%20of%20corrs%20Prevost%20PsyMeth%2008.pdf)

Let's import the data from Table 1 in the article.



```r
library(MASS)  # for simulating multivariate normal data
table1 <- read.csv("data/table1.csv")
```




# Function for generating correlation matrices
To make conversion between vectors of correlations and correlation matrices, I wrote the following function.



```r
vector2rmatrix <- function(x) {
    # x: vector of correlations representing either (a) upper off-diagonal
    # elements of matrix arranged row-wise, or equivalently (b) lower
    # off-diagnoal elements of matrix arrange column-wise.
    V <- length(x)
    p <- (1 + sqrt(1 + 8 * V))/2
    if (p%%1 != 0) 
        stop("input vector has invalid length")
    r_matrix <- diag(rep(1, p))
    r_matrix[lower.tri(r_matrix)] <- x
    r_matrix[upper.tri(r_matrix)] <- t(r_matrix)[upper.tri(r_matrix)]
    r_matrix
}
```





# Simulation to assess asymptotic correlations
In order to understand the properties of correlations of correlation estimates and in order to check the results of exact solutions, I developed the following simulation.

1. It generates a data frame based on a specified correlation matrix and sample size.
2. Sample correlations are obtained for the sample
3. The simulation is iterated many times, and the sample correlations are correlated. 

I assume that this correlation between sample estimates of correlations provides a approximate numeric estimate of the asymptotic correlation of a pair of correlations



```r
generate_sample <- function(n = 100, r = vector2rmatrix(c(0.1, 0.2, 
    0.3, 0.4, 0.5, 0.6)), j = 1, k = 2, l = 3, m = 4) {
    # n: sample size for simulation r: correlation matrix j,k,l,m: return
    # correlation between variables j with k and l with m return two
    # correlations
    x <- mvrnorm(n = n, mu = rep(0, nrow(r)), Sigma = r)
    cor_x <- cor(x)
    c(cor_x[j, k], cor_x[l, m])
}

simulate_correlations <- function(k = 2000, n = 100, r = vector2rmatrix(c(0.1, 
    0.2, 0.3, 0.4, 0.5, 0.6))) {
    # k: the number of simulations to run n: sample size for each run r: true
    # correlation matrix return
    results <- lapply(seq(k), function(X) generate_sample(n = n, r = r))
    plot(sapply(results, function(X) X[1]), sapply(results, function(X) X[2]), 
        pch = ".")
    cor(sapply(results, function(X) X[1]), sapply(results, function(X) X[2]))
}
```






```r
simulate_correlations(k = 10000, n = 10, r = vector2rmatrix(rep(0, 
    6)))
```

![plot of chunk simulated_estimates](figure/simulated_estimates1.png) 

```
## [1] 0.00696
```

```r
simulate_correlations(k = 10000, n = 10, r = vector2rmatrix(c(0.9, 
    0, 0, 0, 0, 0.9)))
```

![plot of chunk simulated_estimates](figure/simulated_estimates2.png) 

```
## [1] -0.007028
```

```r
simulate_correlations(k = 10000, n = 10, r = vector2rmatrix(rep(0.95, 
    6)))
```

![plot of chunk simulated_estimates](figure/simulated_estimates3.png) 

```
## [1] 0.4951
```

```r
simulate_correlations(k = 10000, n = 10, r = vector2rmatrix(rep(0.6, 
    6)))
```

![plot of chunk simulated_estimates](figure/simulated_estimates4.png) 

```
## [1] 0.2477
```

```r
simulate_correlations(k = 10000, n = 1000, r = vector2rmatrix(rep(0.6, 
    6)))
```

![plot of chunk simulated_estimates](figure/simulated_estimates5.png) 

```
## [1] 0.2774
```




This simulation suggests that

* If two correlaions are based on variables that are independent of each other, then the estimates will be independent.
* Increasing the dependence, increases the correlation of the estimates.
* The sample size does not appear to be substantially related (if at all) to the correlation between estimates. The covariance of the estimates decreases with $n$, but so does the variance of both correlation estimates. Given that $\textrm{cor}(XY)=\frac{\textrm{cov}(XY)}{\sqrt{\textrm{var}(X)\textrm{var}(Y)}}$, changes in sample size should basically cancel out. But I need to formally look into the details.

# Prevost et al asymptotic correlation between sample correlations


```r
prevost_corr <- function(r, n, j = 1, k = 2, l = 3, m = 4, with_n_minus_3 = FALSE) {
    # r_matrix is a correlation matrix of size 4
    sigma_jk_lm <- r[j, l] * r[k, m] + r[j, m] * r[k, l] - r[l, m] * (r[j, l] * 
        r[k, l] + r[j, m] * r[k, m]) - r[j, k] * (r[j, l] * r[j, m] + r[k, l] * 
        r[k, m]) + 1/2 * r[j, k] * r[l, m] * (r[j, l]^2 + r[j, m]^2 + r[k, l]^2 + 
        r[k, m]^2)
    
    if (with_n_minus_3) {
        # formula that appears in the paper
        result <- sigma_jk_lm/((n - 3) * (1 - r[j, k]^2) * (1 - r[l, m]^2))
    } else {
        # formula that I think is correct
        result <- sigma_jk_lm/((1 - r[j, k]^2) * (1 - r[l, m]^2))
    }
    result
}
```




The above code implements Equations 1 and 2 in Prevost et al. However, it appears to me that the $n-3$ should not be in the denominator of Equation 1 in Prevost et al. Thus, I have altered the code to not use it.


# Olkins and Finn formula
Because  the Prevost et al formula was disagreeing with my simulation, I read Olkins and Finn *Correlations Redux* to try another approach.



```r
olkins_finn_corr <- function(r, n, i = 1, j = 2, k = 3, l = 4) {
    var_r <- function(r, n) (1 - r^2)^2/n
    cov_r <- function(r, n, i = i, j = j, k = k, l = l) {
        (1/2 * r[i, j] * r[k, l] * (r[i, k]^2 + r[i, l]^2 + r[j, k]^2 + r[j, 
            l]^2) + r[i, k] * r[j, l] + r[i, l] * r[j, k] - (r[i, j] * r[i, 
            k] * r[i, l] + r[j, i] * r[j, k] * r[j, l] + r[k, i] * r[k, j] * 
            r[k, l] + r[l, i] * r[l, j] * r[l, k]))/n
    }
    v1 <- var_r(r[i, j], 10)
    v2 <- var_r(r[k, l], 10)
    c1 <- cov_r(r, 10, i, j, k, l)
    c1/sqrt(v1 * v2)
}
```







```r
test_matrix <- vector2rmatrix(rep(0.5, 6))
for (i in 1:4) {
    print(simulate_correlations(10000, 100, test_matrix))
}
```

```
## [1] 0.2327
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-61.png) 

```
## [1] 0.2236
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-62.png) 

```
## [1] 0.2181
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-63.png) 

```
## [1] 0.2298
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-64.png) 

```r

olkins_finn_corr(test_matrix, n = 100)
```

```
## [1] 0.2222
```

```r
prevost_corr(test_matrix, n = 100, with_n_minus_3 = FALSE)
```

```
## [1] 0.2222
```

```r
prevost_corr(test_matrix, n = 100, with_n_minus_3 = TRUE)
```

```
## [1] 0.002291
```




* My formula based on Olkins and Finn (1995) agrees with the Prevost formula when I remove the $n-3$ from the denominator.
* The mean correlations from the simulations are similar to the values returned by the exact formulas.


# Analysis of Table 1
Let's replicate the paragraph that assesses the asymptotic correlation for all pairs of correlations in a set of correlation matrices provided in the paper.




```r
extract_matrix <- function(x, row, var_names = c("pb", "a", "sn", 
    "pbc", "i", "fb")) {
    id <- table1[row, "study"]
    n <- table1[row, "ni"]
    r_vector <- table1[row, -c(1, 2)]
    r_vector <- unlist(r_vector)
    r_matrix <- vector2rmatrix(r_vector)
    rownames(r_matrix) <- colnames(r_matrix) <- var_names
    
    list(id = id, n = n, r_matrix = r_matrix, r_vector = r_vector)
}

study1 <- extract_matrix(table1, 1)
prevost_corr(study1$r_matrix, study1$n, j = 1, k = 2, l = 3, m = 4)
```

```
## [1] 0.06451
```




The above code extracts the correlation matrix and sample size from a row of the correlation data provided in the paper which I have stored in teh data frame `table1`.



```r
correlation_pairs <- function(p) {
    result <- expand.grid(j = seq(p), k = seq(p), l = seq(p), m = seq(p))
    result <- result[result$j < result$k & result$l < result$m, ]
    result <- result[(p * result$j + result$k) < (p * result$l + result$m), 
        ]
    
    result$common_variable <- result$k == result$l | result$j == result$l | 
        result$k == result$m
    
    result
}

p <- nrow(study1$r_matrix)
indices <- correlation_pairs(p)

cat("common count: ", sum(indices$common_variable), "\nno common variable count: ", 
    sum(!indices$common_variable), "\n")
```

```
## common count:  60 
## no common variable count:  45 
```

```r
head(indices)
```

```
##     j k l m common_variable
## 439 1 2 1 3            TRUE
## 475 1 2 2 3            TRUE
## 481 1 3 2 3            TRUE
## 487 1 4 2 3           FALSE
## 493 1 5 2 3           FALSE
## 499 1 6 2 3           FALSE
```




The above code is used to generate the a data frame that includes the indices for all unique pairs of correlations where the indices represent the variable numbers of the two variables in the first correlation and the two variables in the second correlation.

It also records whether the pair of correlations includes a common variable, which is relevant to the subsequent analysis.



```r
extract_asymptotic_r <- function(row_number) {
    study <- extract_matrix(table1, row_number)
    p <- nrow(study$r_matrix)
    indices <- correlation_pairs(p)
    indices$study <- row_number
    indices$asymptotic_correlation <- sapply(seq(nrow(indices)), function(X) prevost_corr(study$r_matrix, 
        study$n, indices[X, "j"], indices[X, "k"], indices[X, "l"], indices[X, 
            "m"]))
    indices
}

combined_indices <- lapply(1:8, function(X) extract_asymptotic_r(X))

combined_indices <- do.call(rbind, combined_indices)

aggregate(asymptotic_correlation ~ common_variable, combined_indices, 
    mean)
```

```
##   common_variable asymptotic_correlation
## 1           FALSE                 0.1405
## 2            TRUE                 0.3232
```

```r

par(mfrow = c(2, 1))
hist(combined_indices[combined_indices$common_variable == TRUE, "asymptotic_correlation"], 
    xlim = c(-0.1, 0.8), main = "Correlation of correlations have common variables")
hist(combined_indices[combined_indices$common_variable == FALSE, 
    "asymptotic_correlation"], xlim = c(-0.1, 0.8), main = "Correlation of correlations have distinct variables")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 


* The above analysis extracts the asymptotic correlations for each pair of correlations for each of the eight study correlation matrices. These are then combined into an overall data.frame.
* The results obtained replicate those reported in Prevost et al (2007) regarding the higher correlation of correlations for pairs of correlations that include a common variable. 






