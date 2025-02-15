# marbiotrends

# Setting up the model
Download CmdStanR and follow the help guide:
https://mc-stan.org/cmdstanr/articles/cmdstanr.html

Stan files are located in the stan folder

# Model (Multi-variate autoregressive state-space model in a Bayesian framework)
An in-depth description about the statistical underpinnings of the model:
https://nwfsc-timeseries.github.io/MARSS-Manual/chap-marss.html
Why Bayesian? Scarce ecological data, some parameters will be difficult to estimate in a likelihood framework.
The priors in the Stan examples are based on pinnipeds, so they will likely need to change for your taxa.
The MARSS models were originally designed to be ran at the species level, simulatneously estimating all conspecific populations (multivariate). This can be modified to higher taxonomic levels if you so wish. Model run times will exponentially increase with # of time series in the MARSS. 

# Script 
The outputs will give you model diagnostics, abundance predictions, and error parameters.
Finally, I plot the fits to the data, and the total abundance. 
