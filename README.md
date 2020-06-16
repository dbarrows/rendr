
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rendr

<!-- badges: start -->

[![R build
status](https://github.com/dbarrows/rendr/workflows/R-CMD-check/badge.svg)](https://github.com/dbarrows/rendr/actions)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

An R package for simulating reaction and reaction-diffusion systems.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
devtools::install_github('dbarrows/rendr')
```

## Reaction system solvers

Networks are created using the
[`bondr`](https://dexter.barrows.io/bondr) package. See [Creating
Networks](https://dexter.barrows.io/bondr/#creating-networks) for more
details.

``` r
library(rendr)
```

``` r
(network <- bondr::network_examples())
#> [2m#  Reaction network: 3 reactions x 4 species[22m
#> [34m     Reactants    Products     Rate
#> [39m[2mR1[22m       S + E [90m->[39m SE        0.00166
#> [2mR2[22m          SE [90m->[39m S + E       1e-04
#> [2mR3[22m          SE [90m->[39m E + P         0.1
```

### RRE

A deterministic solver that uses the Reaction Rate Equation (RRE) in
conjunction with an ode solver.

``` r
(sys <- rsys(network = network,
             state = c(301, 120, 0, 0),
             T = 30))
#> [90m$network[39m
#> [2m#  Reaction network: 3 reactions x 4 species[22m
#> [34m     Reactants    Products     Rate
#> [39m[2mR1[22m       S + E [90m->[39m SE        0.00166
#> [2mR2[22m          SE [90m->[39m S + E       1e-04
#> [2mR3[22m          SE [90m->[39m E + P         0.1
#> 
#> [90m$state[39m
#>   S   E  SE   P 
#> 301 120   0   0 
#> 
#> [90m$T[39m
#> [1] 30
```

``` r
(sol <- rre(sys))
#> [90mNetwork[39m
#> [2m#  Reaction network: 3 reactions x 4 species[22m
#> [34m     Reactants    Products     Rate
#> [39m[2mR1[22m       S + E [90m->[39m SE        0.00166
#> [2mR2[22m          SE [90m->[39m S + E       1e-04
#> [2mR3[22m          SE [90m->[39m E + P         0.1
#> 
#> [90mSolution[39m
#> [90m# A tibble: 100 x 5[39m
#>     Time     S     E    SE      P
#>    [3m[90m<dbl>[39m[23m [3m[90m<dbl>[39m[23m [3m[90m<dbl>[39m[23m [3m[90m<dbl>[39m[23m  [3m[90m<dbl>[39m[23m
#> [90m 1[39m 0      301  120     0    0    
#> [90m 2[39m 0.303  285. 104.   16.2  0.255
#> [90m 3[39m 0.606  271.  90.9  29.1  0.947
#> [90m 4[39m 0.909  260.  80.6  39.4  1.99 
#> [90m 5[39m 1.21   250.  72.1  47.9  3.32 
#> [90m 6[39m 1.52   241.  65.2  54.8  4.88 
#> [90m 7[39m 1.82   234.  59.5  60.5  6.62 
#> [90m 8[39m 2.12   227.  54.8  65.2  8.53 
#> [90m 9[39m 2.42   221.  50.9  69.1 10.6  
#> [90m10[39m 2.73   216.  47.6  72.4 12.7  
#> [90m# … with 90 more rows[39m
```

A function is provided for easy visualisation of solutions.

``` r
ggplot2::theme_set(wplot::theme_wc())
plot(sol)
```

<img src="man/figures/README-unnamed-chunk-6-1.svg" width="100%" />

### SSA

Generate a single realisation of the Chemical Master Equation (CME)
solution via the Stochastic Solution Algorithm (SSA).

``` r
ssa(sys, all.out = TRUE) %>% plot()
```

<img src="man/figures/README-ssa-1.svg" width="100%" />

## Reaction-diffusion system solvers

Reaction-diffusion systems created via the `rdsys` class, which are
constructed similarly to `rsys`s.

``` r
(sys <- rdsys(
    network = network('
             0 <-> U,  4e3, 2
             0  -> V,  1.2e4
        2U + V  -> 3U, 12.5e-8
    '),
    volume = volume(
        dims = c(40, 1, 1),
        h = 1/40,
        seed = c(25, 75)
    ),
    D = c(1e-3, 1e-1),
    T = 3.5
))
#> [90m$network[39m
#> [2m#  Reaction network: 4 reactions x 2 species[22m
#> [34m     Reactants    Products      Rate
#> [39m[2mR1[22m           0 [90m->[39m U             4000
#> [2mR2[22m           U [90m->[39m 0                2
#> [2mR3[22m           0 [90m->[39m V            12000
#> [2mR4[22m      2U + V [90m->[39m 3U        1.25e-07
#> 
#> [90m$volume[39m
#> [2m#[22m[34m dims: [39m40 x 1 x 1
#> [2m#[22m[34m h: [39m0.025
#> [2m#[22m[34m states: [39m[2m
#> # [22m[90m# A tibble: 40 x 5[39m
#>       x     y     z     U     V
#>   [3m[90m<int>[39m[23m [3m[90m<int>[39m[23m [3m[90m<int>[39m[23m [3m[90m<dbl>[39m[23m [3m[90m<dbl>[39m[23m
#> [90m1[39m     1     1     1    25    75
#> [90m2[39m     2     1     1    25    75
#> [90m3[39m     3     1     1    25    75
#> [90m4[39m     4     1     1    25    75
#> [90m5[39m     5     1     1    25    75
#> [90m# … with 35 more rows[39m
#> 
#> [90m$D[39m
#>     U     V 
#> 0.001 0.100 
#> 
#> [90m$T[39m
#> [1] 3.5
```

### ISSA

Generate a single realisation of the Reaction-diffusion Master Equation
(RDME) solution via the Inhomogeneous Stochastic Solution Algorithm
(ISSA).

``` r
issa(sys) %T>%
    print() %>%
    plot()
#> Starting ISSA simulation with parameters:
#>  - Reactions:   4
#>  - Species:     2
#>  - Dimensions:  40x1x1
#>  - h:           0.025
#>  - time:        [0, 3.5]
#> ....................................................................................................
#> [90mNetwork[39m
#> [2m#  Reaction network: 4 reactions x 2 species[22m
#> [34m     Reactants    Products      Rate
#> [39m[2mR1[22m           0 [90m->[39m U             4000
#> [2mR2[22m           U [90m->[39m 0                2
#> [2mR3[22m           0 [90m->[39m V            12000
#> [2mR4[22m      2U + V [90m->[39m 3U        1.25e-07
#> 
#> [90mSolution[39m
#> [2m# 100 time points x 2 species[22m
```

<img src="man/figures/README-issa-1.svg" width="100%" />

### NSM

Generate a single realisation of the (RDME) solution via the Next
Subvolume Method (NSM). Usage of the `nsm` function is the same as with
the [ISSA solver](#issa). The NSM algorithm is usually faster for
systems with a large number of subvolume (voxels) relative to the
reaction network size. You may have to try both the ISSA and NSM solver
to see which is faster for a given reaction-diffusion system.

``` r
system.time(issa(sys, verbose = FALSE))
#>    user  system elapsed 
#> 249.551   0.201 250.134
system.time(nsm(sys, verbose = FALSE))
#>    user  system elapsed 
#> 139.534   0.124 139.795
```

``` r
sys_small <- sys
sys_small$volume <- volume(
        dims = c(2, 1, 1),
        h = 1/40,
        seed = c(25, 75)
    )
system.time(issa(sys_small, verbose = FALSE))
#>    user  system elapsed 
#>   0.629   0.003   0.632
system.time(nsm(sys_small, verbose = FALSE))
#>    user  system elapsed 
#>   1.042   0.003   1.046
```
