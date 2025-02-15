# Setting up the model
Download CmdStanR and follow the help guide:
https://mc-stan.org/cmdstanr/articles/cmdstanr.html

Stan files are located in the stan folder

# Model (MARSS in a Bayesian framework)
## Why state-space?
SSMs differentiate process and observation error, which allows estimation of unobserved states while simultaneously accounting for sampling uncertainty. In the context of population dynamics, process error refers to true changes in abundance due to birth, deaths, immigration, emigrations, while observation error refers to randomness in surveys sampling abundance. This separation of process and observation error may reduce spurious interpretations of extinction risk, which should be based only on process variability

An in-depth description about the statistical underpinnings of the multivariate (multi-population) autoregressive state-space model:
https://nwfsc-timeseries.github.io/MARSS-Manual/chap-marss.html

## Why Bayesian? Scarce ecological data, some parameters will be difficult to estimate in a likelihood framework.
The priors in the Stan examples are based on pinnipeds, so they will likely need to change for your taxa.

## Framework
The MARSS models were originally designed to be ran at the species level, simulatneously estimating all conspecific populations (multivariate). This can be modified to higher taxonomic levels if you so wish. Model run times will exponentially increase with # of time series in the MARSS.

## Why not just fit a smoother like a GAM over each time series?
MARSS essentially does this, but the cool thing is squeezing as much information from the data as possible. Information meaning process and error parameters. We can also explore more hypotheses about different observation params by each survey method, or different process params by subspecies, or covariance of process changes based on spatial distance... which are all beyond the scope of this class.

# Script 
The outputs will give you model diagnostics, abundance predictions, and error parameters. I plot the fits to the data, and the total abundance. 
