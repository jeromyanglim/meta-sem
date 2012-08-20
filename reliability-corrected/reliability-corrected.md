The following implements some of the functions from Hunter and Schmidt (2004) page 95 onwards on "error of measurement and correction for attenuation.

# Correction based on true 


```r
correct_true_correlation <- function(rho_xy, rho_xx, rho_yy) {
    # calculate reliabilty based on population correlations and reliabilities
    stopifnot(rho_xy >= -1 & rho_xy <= 1)
    stopifnot(rho_xx >= 0 & rho_xx <= 1)
    stopifnot(rho_yy >= 0 & rho_yy <= 1)
    stopifnot(rho_yy >= 0 & rho_yy <= 1)
    stopifnot(rho_xy <= (sqrt(rho_xx) * sqrt(rho_yy)))
    
    rho_xy/(sqrt(rho_xx) * sqrt(rho_yy))
}

correct_true_correlation(0.5, 0.8, 0.8)
```

```
## [1] 0.625
```




# Uncorrect true correlation
Alternatively, the corrected true correlation may be known along with the reliability of observed x and observed y. From this, the uncorrected true correlation can be calculated.



```r
uncorrect_true_correlation <- function(rho_UT, rho_xx, rho_yy) {
    sqrt(rho_xx) * sqrt(rho_yy) * rho_UT
}

uncorrect_true_correlation(0.5, 0.8, 0.8)
```

```
## [1] 0.4
```

```r
correct_true_correlation(0.4, 0.8, 0.8)  # restore
```

```
## [1] 0.5
```




# Standard error of correlation


```r
standard_error_correlation <- function(r, n) {
    (1 - r^2)/sqrt(n - 1)
}

ci_correlation <- function(r, n, ci = 0.95) {
    z <- qnorm(1 - (1 - ci)/2)
    se_r <- standard_error_correlation(r, n)
    ci_diff <- z * se_r
    ci_r <- r + c(-1, 1) * ci_diff
    names(ci_r) <- c("lower_ci", "upper_ci")
    ci_r
}

ci_correlation_reliability_correct <- function(r, n, rho_xx, rho_yy, 
    ci = 0.95) {
    ci_r <- ci_correlation(r, n, ci = ci)
    ci_rho <- correct_true_correlation(ci_r[1], rho_xx, rho_yy)
    ci_rho[2] <- correct_true_correlation(ci_r[2], rho_xx, rho_yy)
    names(ci_rho) <- c("lower_ci", "upper_ci")
    ci_rho
}
```





# Page 97 Example


```r
rho_TU <- 0.6
rho_xx <- 0.45
rho_yy <- 0.55
N <- 100
r_xy <- 0.2

# calculate true correlation
(rho_xy <- uncorrect_true_correlation(rho_TU, rho_xx, rho_yy))
```

```
## [1] 0.2985
```

```r

# percentage reduction from corrected to uncorrected
round(((rho_TU - rho_xy)/rho_TU), 3)
```

```
## [1] 0.503
```

```r

# standard error of uncorrected correlation
standard_error_correlation(rho_xy, N)
```

```
## [1] 0.09155
```

```r

# standard error of uncorrected observed correlation
standard_error_correlation(r_xy, N)
```

```
## [1] 0.09648
```

```r

# ci uncorrected observed
round(ci_correlation(r_xy, N, 0.95), 2)
```

```
## lower_ci upper_ci 
##     0.01     0.39 
```

```r

# ci corrected
round(ci_correlation_reliability_correct(r_xy, N, rho_xx, rho_yy, 
    0.95), 2)
```

```
## lower_ci upper_ci 
##     0.02     0.78 
```



