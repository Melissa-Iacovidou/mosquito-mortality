# Omitting age-dependent mosquito mortality in malaria models underestimates the effectiveness of insecticide-treated nets

### Authors: M.A. Iacovidou<sup>1,2</sup>, P. Barreaux<sup>3</sup>, S.E.F. Spencer<sup>2,4</sup>, M.B. Thomas<sup>5</sup>, E.E. Gorsich<sup>2,6</sup>, K.S. Rock<sup>1,2</sup>

**1** Mathematics Institute, University of Warwick, UK \
**2** The Zeeman Institute for Systems Biology and Infectious Disease Epidemiology Research, University of Warwick, UK \
**3** Liverpool School of Tropical Medicine, UK \
**4** Department of Statistics, University of Warwick, UK \
**5** Department of Biology, University of York, UK \
**6** School of Life Sciences, University of Warwick, UK


This repository contains the code for the paper "Omitting age-dependent mosquito mortality in malaria models underestimates the effectiveness of insecticide-treated nets".

The code is written in **Julia** and requires the following packages:

<table>
<td>
  <p> XLSX.jl </p>
  <p> Plots.jl </p>
  <p> DataFrames.jl </p>
</td>
<td>
  <p> Survival.jl </p>
  <p> Optim.jl </p>
  <p> Distributions.jl </p>
</td>
<td>
  <p> Measurements.jl </p>
  <p> LinearAlgebra.jl </p>
  <p> Random.jl </p>
</td>
<td>
  <p> Statistics.jl </p>
  <p> LaTeXStrings.jl </p>
  <p> StatsPlots.jl </p>
</td>
<td>
  <p> QuadGK.jl </p>
  <p> Cuba.jl </p>
</td>
</table>

### Files
**Modules**
- _ParameterFitting.jl_ - Fits the data to all the functions using maximum likelihood estimation. Calculates Wald confidence intervals and sets up distributions for Monte Carlo simulations. Exports all the functions and fitted variables needed for the plots in _SurvivalMortalityPlots.jl_ and _VectorialCapacityPlots.jl_. 
- _VCCalculations.jl_ - Exports all the functions used in the four steps for the calculation of the vectorial capacity, in addition to the calculation of the percentage decrease between the two treatments. These are used in plotting the various figures in _VectorialCapacityPlots.jl_.

**Plotting**
- _InitialDataPlots.jl_ - Plots the figures that can be found in the Supporting information and the initial data figure with all the replicates added together.
- _SurvivalMortalityPlots.jl_ - Plots the various survival and mortality figures found in the paper using the module _ParameterFitting.jl_.
- _VectorialCapacityPlots.jl_ - Plots the figures for the calculations in steps 1-4 regarding the vectorial capacity.

In addition, a Mathematica file (_ErlangEIP.nb_) is included where the values for the EIP distribution variables are obtained, which are used in the calculations.

### Notes
⚠️ For the code to run one must download the data from the Supporting information and edit the relevant lines in _ParameterFitting.jl_ & _InitialDataPlots.jl_.

⚠️ The plots are further edited using **Inkscape**, hence there may be some discrepancies in the visual appearance of the figures.
