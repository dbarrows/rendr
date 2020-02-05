
<!-- README.md is generated from README.Rmd. Please edit that file -->

# reactor

<!-- badges: start -->

<!-- badges: end -->

An R package for simulating reaction and reaction-diffusion systems.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("dbarrows/reactor")
```

## Reaction system solvers

Networks are created using the
[`bondr`](https://github.com/dbarrows/bondr) package. See [Creating
Networks](https://github.com/dbarrows/bondr#creating-networks) for more
details.

``` r
library(reactor)

(mm_network <- parse_network(bondr::mm_string))
#> # Reaction network: 3 reactions x 4 species
#>     Reactants    Products     Rate
#> 1       S + E -> SE        1.66e-3
#> 2          SE -> S + E        1e-4
#> 3          SE -> E + P        1e-1
```

### RRE

A deterministic solver that uses the Reaction Rate Equation (RRE) in
conjunction with an ode solver.

``` r
(model <- rmodel(network = mm_network,
                 state = c(301, 120, 0, 0),
                 tspan = c(0, 30)))
#> $network
#> # Reaction network: 3 reactions x 4 species
#>     Reactants    Products     Rate
#> 1       S + E -> SE        1.66e-3
#> 2          SE -> S + E        1e-4
#> 3          SE -> E + P        1e-1
#> 
#> $state
#>   S   E  SE   P 
#> 301 120   0   0 
#> 
#> $tspan
#> [1]  0 30
```

``` r
(sol <- rre(model))
#> # A tibble: 100 x 5
#>     Time     S     E    SE      P
#>    <dbl> <dbl> <dbl> <dbl>  <dbl>
#>  1 0      301  120     0    0    
#>  2 0.303  285. 104.   16.2  0.255
#>  3 0.606  271.  90.9  29.1  0.947
#>  4 0.909  260.  80.6  39.4  1.99 
#>  5 1.21   250.  72.1  47.9  3.32 
#>  6 1.52   241.  65.2  54.8  4.88 
#>  7 1.82   234.  59.5  60.5  6.62 
#>  8 2.12   227.  54.8  65.2  8.53 
#>  9 2.42   221.  50.9  69.1 10.6  
#> 10 2.73   216.  47.6  72.4 12.7  
#> # … with 90 more rows
```

A function is provided for easy visualisation of solutions.

``` r
rsol_plot(sol)
```

<img src="man/figures/README-unnamed-chunk-5-1.svg" width="100%" />

### SSA

Generate a single realisation of the Chemical Master Equation (CME)
solution via the Stochastic Solution Algorithm (SSA).

``` r
ssa(model) %>%
    rsol_plot()
```

<img src="man/figures/README-unnamed-chunk-6-1.svg" width="100%" />

## Reaction-diffusion system solvers

Reaction-diffusion models created via the `rdmodel` class, which are
constructed similarly to `rmodel`s.

``` r
(model <- rdmodel(
    network = parse_network("
             0 <-> U,  4e3, 2
             0  -> V,  1.2e4
        2U + V  -> 3U, 12.5e-8
    "),
    volume = volume(
        dims = c(40, 1, 1),
        h = 1/40,
        seed = c(25, 75)
    ),
    D = c(1e-3, 1e-1),
    tspan = c(0, 3.5)
))
#> $network
#> # Reaction network: 4 reactions x 2 species
#>     Reactants    Products     Rate
#> 1           0 -> U             4e3
#> 2           U -> 0               2
#> 3           0 -> V           1.2e4
#> 4      2U + V -> 3U        12.5e-8
#> 
#> $volume
#> # dims: 40 x 1 x 1
#> # h: 0.025
#> # states: 
#> # # A tibble: 40 x 5
#>       x     y     z     U     V
#>   <int> <int> <int> <dbl> <dbl>
#> 1     1     1     1    25    75
#> 2     2     1     1    25    75
#> 3     3     1     1    25    75
#> 4     4     1     1    25    75
#> 5     5     1     1    25    75
#> # … with 35 more rows
#> 
#> $D
#>     U     V 
#> 0.001 0.100 
#> 
#> $tspan
#> [1] 0.0 3.5
```

### ISSA

Generate a single realisation of the Reaction-diffusion Master Equation
(RDME) solution via the Inhomogeneous Stochastic Solution Algorithm
(ISSA).

``` r
issa(model) %>%
    rdsol_plot()
#> Starting ISSA simulation with parameters:
#>  - Reactions:   4
#>  - Species:     2
#>  - Dimensions:  40x1x1
#>  - h:           0.025
#>  - time: [0, 3.5]
#> ....................................................................................................
```

<img src="man/figures/README-unnamed-chunk-8-1.svg" width="100%" />
