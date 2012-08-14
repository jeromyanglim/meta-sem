# Random-effects meta analysis

The following document implements the random-effects meta-analysis procedure outlined by Borenstein, Hedges, Higgins, and Rothstein (2009) pages 72 to 75.



```r
set.seed(1234)
library(metafor)
```





# Function to generate data
The data generating mechanism is

$$\begin{align}
Y_i  & = \mu + \zeta_i + \epsilon_i \\
\zeta_i & \sim N(0, \tau^2) \\
\epsilon_i & \sim N(0, V_{Y_i})
\end{align}
$$

where 

* $Y_i$ is the observed effect for the $i$th study, 
* $\mu$ is the overall mean effect, 
* $\zeta_i$ is the true deviation of the effect of $i$th study from $\mu$, and 
* $\tau^2$ is the variance of the study deviation effects
* $\epsilon_i$ is the deviation of the $i$th study from the studies true effect.
* $V_{Y_i}$ is the variance of the effect due to sampling error.

Note that I have assumed that $\zeta_i$ and $\epsilon_i$ both follow a normal distribution. Other distributions would be possible.




```r
sample_r_data_random <- function(true_mu, true_tau, n) {
    # true_r: true population correlations n: vector of sample sizes for each
    # study
    rcor <- function(true_r, n) {
        x <- rnorm(n)
        y <- rnorm(n, true_r * x, sqrt(1 - true_r^2))
        cor(x, y)
    }
    
    # true_r
    true_r <- rnorm(n = length(n), mean = true_mu, sd = true_tau)
    
    # sample correlations
    Y <- sapply(seq(n), function(X) rcor(true_r[X], n = n[X]))
    
    # note that this is only approximate
    approx_variance <- function(r, n) (1 - r^2)^2/(n - 1)
    
    # Variance of correlation
    Vy <- approx_variance(Y, n)
    
    # Standard error of correlation
    SEy <- sqrt(Vy)
    
    data.frame(n, true_r, Y, Vy, SEy)
}

```





# Performing a random-effects meta-analysis


```r
random_effects_meta_analysis <- function(Y, Vy, ci = 0.95) {
    # Y: vector of effect sizes (e.g., r correlations) Vy: vector ci: scalar,
    # confidence interval: 0 < ci < 1
    stopifnot(length(Y) == length(Vy))
    
    get_tau_square <- function(Y, Vy) {
        W <- 1/Vy
        df <- length(Y) - 1
        Q <- sum(W * Y^2) - sum(W * Y)^2/sum(W)
        C <- sum(W) - sum(W^2)/sum(W)
        T_square <- (Q - df)/C
        list(T_square = T_square, W = W, Q = Q, C = C, df = df)
    }
    
    T_square <- get_tau_square(Y, Vy)
    
    
    Vy_star <- Vy + T_square$T_square
    W_star <- 1/Vy_star
    M_star <- sum(W_star * Y)/sum(W_star)
    
    # variance of summary effect
    Vm_star <- 1/sum(W_star)
    
    # standard error of summary effect
    SEm_star <- sqrt(Vm_star)
    
    # Confidence interval for summary effect
    SEm_multiple <- abs(qnorm((1 - ci)/2))
    LLm_star <- M_star - SEm_star * SEm_multiple
    ULm_star <- M_star + SEm_star * SEm_multiple
    
    # z-value for test of null hypothesis
    Z_star <- M_star/SEm_star
    
    # two tailed p-value
    p_star <- 2 * (1 - pnorm(abs(Z_star)))
    
    # return values
    list(studies = cbind(Y, Vy, W = T_square$W, SEy = sqrt(Vy), W_star), T_square = T_square, 
        M_star = M_star, Vm_star = Vm_star, SEm_star = SEm_star, ci = ci, ci_limits = c(LLm_star, 
            ULm_star), test = c(Z_star = Z_star, p_star = p_star))
}

meta1 <- sample_r_data_random(0.3, 0.1, rep(100, 5))
random_effects_meta_analysis(meta1$Y, meta1$Vy)
```

```
## $studies
##             Y       Vy      W     SEy W_star
## [1,]  0.13927 0.009713 102.96 0.09855  33.08
## [2,]  0.34567 0.007831 127.69 0.08849  35.28
## [3,]  0.35236 0.007749 129.06 0.08803  35.38
## [4,] -0.05902 0.010031  99.69 0.10015  32.74
## [5,]  0.30492 0.008310 120.34 0.09116  34.69
## 
## $T_square
## $T_square$T_square
## [1] 0.02052
## 
## $T_square$W
## [1] 102.96 127.69 129.06  99.69 120.34
## 
## $T_square$Q
## [1] 13.49
## 
## $T_square$C
## [1] 462.5
## 
## $T_square$df
## [1] 4
## 
## 
## $M_star
## [1] 0.2215
## 
## $Vm_star
## [1] 0.005842
## 
## $SEm_star
## [1] 0.07643
## 
## $ci
## [1] 0.95
## 
## $ci_limits
## [1] 0.07169 0.37131
## 
## $test
##   Z_star   p_star 
## 2.897912 0.003757 
## 
```





# Testing on examples from book


```r
# Table 14.2
random_effects_meta_analysis(Y = c(0.095, 0.277, 0.367, 0.664, 0.462, 
    0.185), Vy = c(0.033, 0.031, 0.05, 0.011, 0.043, 0.023), ci = 1 - 2 * (1 - 
    pnorm(1.96)))
```

```
## $studies
##          Y    Vy     W    SEy W_star
## [1,] 0.095 0.033 30.30 0.1817  14.43
## [2,] 0.277 0.031 32.26 0.1761  14.86
## [3,] 0.367 0.050 20.00 0.2236  11.59
## [4,] 0.664 0.011 90.91 0.1049  21.15
## [5,] 0.462 0.043 23.26 0.2074  12.61
## [6,] 0.185 0.023 43.48 0.1517  16.87
## 
## $T_square
## $T_square$T_square
## [1] 0.03628
## 
## $T_square$W
## [1] 30.30 32.26 20.00 90.91 23.26 43.48
## 
## $T_square$Q
## [1] 11.74
## 
## $T_square$C
## [1] 185.9
## 
## $T_square$df
## [1] 5
## 
## 
## $M_star
## [1] 0.3577
## 
## $Vm_star
## [1] 0.01093
## 
## $SEm_star
## [1] 0.1045
## 
## $ci
## [1] 0.95
## 
## $ci_limits
## [1] 0.1528 0.5626
## 
## $test
##    Z_star    p_star 
## 3.4216216 0.0006225 
## 
```

```r

# Table 14.9
random_effects_meta_analysis(Y = c(0.5493, 0.6931, 0.4236, 0.2027, 
    0.8673, 0.4847), Vy = c(0.027, 0.0115, 0.0455, 0.0025, 0.0175, 0.0213), 
    ci = 1 - 2 * (1 - pnorm(1.96)))
```

```
## $studies
##           Y     Vy      W    SEy W_star
## [1,] 0.5493 0.0270  37.04 0.1643  9.170
## [2,] 0.6931 0.0115  86.96 0.1072 10.689
## [3,] 0.4236 0.0455  21.98 0.2133  7.840
## [4,] 0.2027 0.0025 400.00 0.0500 11.827
## [5,] 0.8673 0.0175  57.14 0.1323 10.045
## [6,] 0.4847 0.0213  46.95 0.1459  9.676
## 
## $T_square
## $T_square$T_square
## [1] 0.08205
## 
## $T_square$W
## [1]  37.04  86.96  21.98 400.00  57.14  46.95
## 
## $T_square$Q
## [1] 36.26
## 
## $T_square$C
## [1] 381
## 
## $T_square$df
## [1] 5
## 
## 
## $M_star
## [1] 0.5328
## 
## $Vm_star
## [1] 0.01688
## 
## $SEm_star
## [1] 0.1299
## 
## $ci
## [1] 0.95
## 
## $ci_limits
## [1] 0.2781 0.7874
## 
## $test
##    Z_star    p_star 
## 4.101e+00 4.114e-05 
## 
```



