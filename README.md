<h1 align="center">
  &nbsp;Autocorrelation-informed home range estimation:<br> a review and practical guide</h1>
<div align="center">

&nbsp;&nbsp;&nbsp;
<a href="mailto:i.simoes-silva@hzdr.de"><img border="0" alt="Email" src="https://assets.dryicons.com/uploads/icon/svg/8007/c804652c-fae4-43d7-b539-187d6a408254.svg" width="35" height="35"></a>&nbsp;&nbsp;&nbsp;
<a href="https://twitter.com/ecoisilva"><img border="0" alt="Twitter" src="https://assets.dryicons.com/uploads/icon/svg/8385/c23f7ffc-ca8d-4246-8978-ce9f6d5bcc99.svg" width="35" height="35"></a>&nbsp;&nbsp;&nbsp;

</div>

### About the manuscript:

Preprint is available on [EcoEvoRxiv](https://ecoevorxiv.org/23wq7/). Click [here](https://ecoevorxiv.org/23wq7/download) to download the full text.

> Home range estimation is a key output from tracking datasets, but the inherent properties of animal movement can lead traditional methods to under- or overestimated their size. **Autocorrelated Kernel Density Estimation (AKDE)** methods were designed to be statistically efficient while explicitly dealing with the complexities and biases of modern movement data, such as *autocorrelation*, *small sample sizes*, and *missing or irregularly sampled data*.

> This repository is a companion piece to our manuscript *"Autocorrelation-informed home range estimation: a review and practical guide"*, and provides:
1. [R tutorial](https://ecoisilva.github.io/AKDE_minireview/code/AKDE_R-tutorial.html) (or as a [.pdf](files/SuppFile2_R-tutorial.pdf) file).
1. Simulation [data](data/data_sims.csv) and [code](code/AKDE_sims.R).

# R tutorial:

> The **AKDE** family of home range estimators will be run using **R software** and the `ctmm` package (Calabrese *et al.*, [(2016)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12559). If you are not familiar with `R`, make sure you follow these steps:

1. Install **R** from <https://www.r-project.org>.
2. Install **RStudio Desktop** from [here](https://rstudio.com/products/rstudio/download/#download) for a graphical interface for R.
3. Install the required `R` packages with the following code:

```r
install.packages("ctmm")
```
We provide a guide to **home range estimation** using the following workflow:

-  **Step 1.** -- Formatting and loading an animal tracking dataset;
-  **Step 2.** -- Checking for the *range residency* assumption;
-  **Step 3.** -- Selecting the best-fit movement model through *model selection*;
-  **Step 4.** -- Feeding a movement model into the *home range estimator*;
-  **Step 5.** -- Evaluating additional *biases*, applying *mitigation measures*.

# Simulation data and code:

> To quantify the level of improvement offered by each mitigation measure and to explore the tradeoff between accuracy and computational cost, we performed a detailed simulation study. For more details, check our manuscript [here](https://www.biorxiv.org/).

![Methods comparison - error and computational cost](files/methods-comparison.png)


# Useful links:

- [Coding Club](https://ourcodingclub.github.io/tutorials.html): collection of coding tutorials with examples in R.
- [`learn` R package](https://rstudio.github.io/learnr/): interactive tutorials for R.

For further documentation on using the `ctmm` package:
- https://cran.r-project.org/web/packages/ctmm/ctmm.pdf
- https://ctmm-initiative.github.io/ctmm/
