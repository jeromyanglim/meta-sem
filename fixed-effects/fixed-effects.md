# Fixed-effects meta analysis

The following document implements the fixed effects meta-analysis procedure outlined by Borenstein, Hedges, Higgins, and Rothstein (2009) pages 65 to 97.


# Generate sample data


```r
set.seed(1234)

sample_r_data <- function(true_r, n) {
    # true_r: true population correlations n: vector of sample sizes for each
    # study
    rcor <- function(true_r, n) {
        x <- rnorm(n)
        y <- rnorm(n, true_r * x, sqrt(1 - true_r^2))
        cor(x, y)
    }
    
    # sample correlations
    Y <- sapply(n, function(X) rcor(true_r, n = X))
    
    # note that this is only approximate
    approx_variance <- function(r, n) (1 - r^2)^2/(n - 1)
    
    # Variance of correlation
    Vy <- approx_variance(Y, n)
    
    # Standard error of correlation
    SEy <- sqrt(Vy)
    
    data.frame(n, Y, Vy, SEy)
}
```





# Performing a fixed effect meta-analysis


```r
fixed_effect_meta_analysis <- function(Y, Vy, ci = 0.95) {
    # Y: vector of effect sizes (e.g., r correlations) Vy: vector ci: scalar,
    # confidence interval: 0 < ci < 1
    stopifnot(length(Y) == length(Vy))
    
    # assign weight to each study
    W <- 1/Vy
    
    # mean effect
    M <- sum(W * Y)/sum(W)
    
    # variance of summary effect
    Vm <- 1/sum(W)
    
    # standard error of summary effect
    SEm <- sqrt(Vm)
    
    # Confidence interval for summary effect
    SEm_multiple <- abs(qnorm((1 - ci)/2))
    LLm <- M - SEm * SEm_multiple
    ULm <- M + SEm * SEm_multiple
    
    # z-value for test of null hypothesis
    Z <- M/SEm
    
    # two tailed p-value
    p <- 2 * (1 - pnorm(abs(Z)))
    
    # return values
    list(studies = cbind(Y, Vy, W), M = M, Vm = Vm, SEm = SEm, ci = ci, ci_limits = c(LLm, 
        ULm), test = c(Z = Z, p = p))
}
```





# Putting it all together


```r
meta1 <- sample_r_data(true_r = 0.3, n = seq(50, 250, 50))
fixed_effect_meta_analysis(meta1$Y, meta1$Vy)
```

```
## $studies
##           Y       Vy      W
## [1,] 0.2019 0.018778  53.25
## [2,] 0.3972 0.007166 139.55
## [3,] 0.1253 0.006502 153.79
## [4,] 0.2323 0.004497 222.35
## [5,] 0.3844 0.002917 342.83
## 
## $M
## [1] 0.2949
## 
## $Vm
## [1] 0.001097
## 
## $SEm
## [1] 0.03312
## 
## $ci
## [1] 0.95
## 
## $ci_limits
## [1] 0.2300 0.3598
## 
## $test
##     Z     p 
## 8.905 0.000 
## 
```

```r

meta2 <- sample_r_data(true_r = 0.1, n = rep(100, 10))
fixed_effect_meta_analysis(meta2$Y, meta2$Vy)
```

```
## $studies
##               Y       Vy      W
##  [1,]  0.123140 0.009797 102.07
##  [2,] -0.008036 0.010100  99.01
##  [3,]  0.100947 0.009896 101.05
##  [4,]  0.105039 0.009879 101.22
##  [5,]  0.070516 0.010001  99.99
##  [6,]  0.138007 0.009720 102.88
##  [7,]  0.053284 0.010044  99.56
##  [8,]  0.015813 0.010096  99.05
##  [9,]  0.062392 0.010023  99.78
## [10,]  0.091153 0.009934 100.67
## 
## $M
## [1] 0.07573
## 
## $Vm
## [1] 0.0009947
## 
## $SEm
## [1] 0.03154
## 
## $ci
## [1] 0.95
## 
## $ci_limits
## [1] 0.01391 0.13754
## 
## $test
##       Z       p 
## 2.40105 0.01635 
## 
```






# Testing on examples from book


```r
# Table 14.2
fixed_effect_meta_analysis(Y = c(0.095, 0.277, 0.367, 0.664, 0.462, 
    0.185), Vy = c(0.033, 0.031, 0.05, 0.011, 0.043, 0.023))
```

```
## $studies
##          Y    Vy     W
## [1,] 0.095 0.033 30.30
## [2,] 0.277 0.031 32.26
## [3,] 0.367 0.050 20.00
## [4,] 0.664 0.011 90.91
## [5,] 0.462 0.043 23.26
## [6,] 0.185 0.023 43.48
## 
## $M
## [1] 0.4093
## 
## $Vm
## [1] 0.004163
## 
## $SEm
## [1] 0.06452
## 
## $ci
## [1] 0.95
## 
## $ci_limits
## [1] 0.2828 0.5357
## 
## $test
##         Z         p 
## 6.343e+00 2.255e-10 
## 
```




The results are basically are approximately those provided on page 90 to 91. 
I used rounded values of $V_M$ and $Y$ from Table 14.2, and I used a more precise value for 95% confidence interval calculation.

