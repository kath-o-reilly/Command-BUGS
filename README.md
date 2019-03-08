# Command-BUGS

This site accompanies the paper in the COMMAND Special Issue on using the BUGS languague in infectious disease modelling.

The associated R code is in the repo, and the datasets used are in the folder.

## Step 1: which BUGS should I use?

See Table 1 in the paper that summarises the content provided within each version of BUGS; WinBUGS, OpenBUGS, JAGS. If you are familiar with using R and RStudio, I would recommend either OpenBUGS or JAGS. From a personal perspective, I find the error messages easier to interpret in JAGS.

For the code in the repo, all have been check in JAGS.

## Step 2: how do I start using BUGS?

Whatever version you decide, installation of the software is required. The following sites are where you download the software, and follow the instructions from there;

WinBUGS

OpenBUGS

JAGS http://mcmc-jags.sourceforge.net/

## Step 3: Check that the software works

The folder 'software checks' contains the same model and dataset but for each version of BUGS.

The model is the simple linear model described in the paper, where the data has been simulated with alpha=5, beta=-1.5 and tau=0.01. The dataset has 100 observations. 

## Step 4: Trialling other examples in the paper
