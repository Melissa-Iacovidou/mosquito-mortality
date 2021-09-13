# Omitting age-dependent mosquito mortality in malaria models underestimates the effectiveness of insecticide-treated nets

### Authors: M.A. Iacovidou<sup>1,2</sup>, P. Barreaux<sup>3</sup>, M.B. Thomas<sup>4</sup>, E.E. Gorsich<sup>2,5</sup>, K.S. Rock<sup>1,2</sup>

**1** Mathematics Institute, University of Warwick, UK \
**2** The Zeeman Institute for Systems Biology and Infectious Disease Epidemiology Research, University of Warwick, UK \
**3** Liverpool School of Tropical Medicine, UK \
**4** Department of Biology, University of York, UK \
**5** School of Life Sciences, University of Warwick, UK


This repository contains the code for the paper "Omitting age-dependent mosquito mortality in malaria models underestimates the effectiveness of insecticide-treated nets".

The code is written in **Julia** and requires the following packages:

<table>
<td>
<p> XLSX.jl </p>
<p> Plots.jl </p>
  </td>
<td>
<p> DataFrames.jl </p>
<p> LsqFit.jl </p>
  </td>
<td>
<p> Measurements.jl </p>
<p> Statistics.jl </p>
  </td>
<td>
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
- _ParameterFitting.jl_ - Fits the data to all the functions. Exports all the functions and fitted variables (with and without errors) needed for the plots in _SurvivalMortalityPlots.jl_ and _VectorialCapacityPlots.jl_. 
- _VCCalculations.jl_ - Exports all the functions used in the four steps for the calculation of the vectorial capacity, in addition to the calculation of the percentage decrease between the two treatments. These are used in plotting the various figures in _VectorialCapacityPlots.jl_.

**Plots**
- _InitialDataPlots.jl_ - Plots the figures that can be found in the Supporting information and the initial data figure with all the replicates added together.
- _SurvivalMortalityPlots.jl_ - Plots the various survival and mortality figures found in the paper using the module _ParameterFitting.jl_.
- _VectorialCapacityPlots.jl_ - Plots the figures for the calculations in steps 1-4 regarding the vectorial capacity.

In addition, a Mathematica file (_ErlangEIP.nb_) is included where the values for the EIP distribution variables are obtained, which are used in the calculations.

### Notes
⚠️ For the code to run one must download the data from the Supporting information and edit the relevant lines in _ParameterFitting.jl_ & _InitialDataPlots.jl_.

⚠️ The plots are further edited using **Inkscape**, hence there may be some discrepancies in the visual appearance of the figures.
