# assam: Approximate and scalable species archetype models

<!-- badges: start -->

<!--[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) -->

<!-- badges: end -->

`assam` is an R package associated with the manuscript "Approximate and scalable finite mixture models for clustering species environmental responses in ecology" by [Hui](https://francishui.netlify.app/), [Menendez](https://findanexpert.unimelb.edu.au/profile/1043920-patricia-menendez-galvan), [Foster](https://scholar.google.com.au/citations?user=WPRiGeYAAAAJ&hl=en), and [Woolley](https://people.csiro.au/W/S/Skip-Woolley), which is currently in preparation/review. Such model-based clustering of species according to their environmental responses, into a small number of so-called "archetypal responses" is also known in ecology as species archetype models, after [Dunstan et al., (2011)](https://doi.org/10.1016/j.ecolmodel.2010.11.030), [Dunstan et al., (2013)](https://link.springer.com/article/10.1007/s13253-013-0146-x) and [Hui et al., (2013)](https://doi.org/10.1890/12-1322.1).

# Installation

Currently, `assam` is available and can be installed from Github with the help of devtools package using:

```         
pak::pkg_install("fhui28/assam")
```

# Getting started

Users are recommended to either:

1.  Check out the examples in the help files for the `assam` and `passam` functions in the package.

2.  Examine the `runsims.R` scripts inside the folders `ms/setting1` and `ms/setting2` and `ms/setting_small`, which are template `R` scripts that can be adapted to run simulation settings 1, 2, and the smaller simulation setting (where asSAMs was also compared with a Bayesian method to fit SAMs using [cmdstanr](https://mc-stan.org/cmdstanr/)), respectively, in the associated manuscript.

    - The `runsims.R` scripts inside the folders `ms/setting1` and `ms/setting2` can be straightforwardly modified to run simulation settings 3 and 4 as well, which are also described in the web appendices of the manuscript.
    - The `ms` folder additionally contains a `GBRseabed_synthetic.RData` data file which mimics the format of the real [Great Barrier Reef Seabed Biodiversity project](https://www.gbif.org/dataset/3b7f6307-2f03-4f25-9f64-6aa4c3a3aea1) community ecology data, sourced from the Australian Institute of Marine Science (AIMS), that was analyzed as part of the manuscript.

# If you find any bugs and issues...

If you find something that looks like a bug/issue, please use Github issues and post it up there. As much as possible, please include in the issue:

1.  A description of the bug/issue;
2.  Paste-able code along with some comments that reproduces the problem e.g., using the [reprex](https://cran.r-project.org/web/packages/reprex/index.html) package. If you also have an idea of how to fix the problem, then that is also much appreciated.
3.  Required data files etc...

Alternatively, please contact the corresponding author at [fhui28\@gmail.com](mailto:fhui28@gmail.com)
