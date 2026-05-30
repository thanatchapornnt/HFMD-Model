# HFMD-Model

## Overview
A climate-forced SEIRS compartmental model to simulate and forecast 
the seasonal dynamics of Hand, Foot, and Mouth Disease (HFMD) 
in Bangkok, Thailand using monthly data from 2011 to 2022.

## Data Sources
- Disease data: DDC Report 506, Department of Disease Control
- Climate data: TerraClimate (temperature, relative humidity, rainfall)

## Methods
- Deterministic SEIRS model solved with deSolve (R)
- Climate forcing: β(t) = β₀ × exp(α₁T + α₂RH + α₃Rain)
- Parameter estimation via Maximum Likelihood Estimation (bbmle::mle2)
- Stochastic simulation: 100 runs with 95% confidence envelope

## Tools
R · deSolve · ggplot2 · dplyr · bbmle · patchwork

## Author
Thanatchaporn Nakthang  
Department of Physics, Faculty of Science, Naresuan University
