# Compute cohort index from age and period index

Compute cohort index from age and period index

## Usage

``` r
coh(agegroup, period, noa, periods_per_agegroup)
```

## Arguments

- agegroup:

  age group index

- period:

  period index

- noa:

  number of age groups in total

- periods_per_agegroup:

  periods per age group

## Value

cohort index

## Examples

``` r
# last agegroup in first period equals first cohort
coh(10, 1, 10, 5)  
#> [1] 1

# first agegroup in last period equals last cohort 
coh(1, 8, 10, 5) 
#> [1] 53
```
