# Reproduces the results in the paper "Exploring the influence of patient variability on propofol target-controlled infusion performance" by Ylva Wahlquist
# Date: 2023-11-02

# Setup the TCI as a QP problem. Solve for a reference patient (35 year old, male, 70 kg, 170 cm) with venous sampling and no opiates added. The reference patient is the same as in the Eleveld model, with given inter-individual variability in the PK parameters (also exist for the PD parameters, but not used here).
# Use the Mahalanobis distance and sample from the α-quantile of the multivariate normal distribution with mean 0 and covariance matrix Q (eye(1/σ_i)) to get N patient models that are within the α-quantile of the reference patient.

using Pkg
cd(@__DIR__)
Pkg.activate("..")

using Distances, Distributions, OSQP, Plots, Random, SparseArrays
using CSV, DataFrames

##
include("expm.jl")
include("tcifunctions.jl")

##############################


Random.seed!(123456)

N = 61      # Nbr of sampling instances
h = 1 / 6   # Sampling period [min]
n_pat_sample = 1000 # nbr of patients to sample and simulate

## Define the reference patient and find u by solving the TCI problem
referencecovariates, pkpdparams, n, Phi, Gamma = getrefpatient(h) # Get the reference patient covariates and PKPD parameters
t, u, BIS_ref, N, r, CeMax, BIS_baseline = TCIreference(pkpdparams, n, Phi, Gamma, N, h) # Solve the TCI problem for the reference patient

plotBISu(t, u, BIS_ref) # Plot the reference patient u and simulated BIS


## Sample patients similar to the reference patient of quantile α with variability in PK parameters
PKvariability = true
PDvariability = false

n_pat_sample = 1000
α = 0.5 # 50th quantile of the multivariate normal distribution to sample from
BIS_all = samplesimulatepatients(h, referencecovariates, BIS_baseline, n_pat_sample, α, PKvariability, PDvariability; seed=123456)

plotsampledpatients(n_pat_sample, BIS_all, t, BIS_ref)

# Save to file
saveresultstofile(α, t, u, BIS_all, BIS_ref) # Save reference u, reference y and sampled patients BIS to file
savehistdatatofile(BIS_all, α, n_pat_sample) # Save histogram data of final BIS (at 10 min) to file

num60, num40, totaloutside = countBISoutside(BIS_all)
percentoutside = totaloutside / n_pat_sample * 100

## Now with α = 0.01

α = 0.01 # 1st quantile of the multivariate normal distribution to sample from
BIS_all = samplesimulatepatients(h, referencecovariates, BIS_baseline, n_pat_sample, α, PKvariability, PDvariability; seed=123456)

plotsampledpatients(n_pat_sample, BIS_all, t, BIS_ref)

# Save to file
saveresultstofile(α, t, u, BIS_all, BIS_ref) # Save reference u, reference y and sampled patients BIS to file
savehistdatatofile(BIS_all, α, n_pat_sample) # Save histogram data of final BIS (at 10 min) to file

num60, num40, totaloutside = countBISoutside(BIS_all)
percentoutside = totaloutside / n_pat_sample * 100


## Sample patients similar to the reference patient of quantile α with variability in PK AND PD parameters
PKvariability = true
PDvariability = true

n_pat_sample = 1000
α = 0.5 # 50th quantile of the multivariate normal distribution to sample from
BIS_all = samplesimulatepatients(h, referencecovariates, BIS_baseline, n_pat_sample, α, PKvariability, PDvariability; seed=123456)

plotsampledpatients(n_pat_sample, BIS_all, t, BIS_ref)

# Save to file
saveresultstofile(α, t, u, BIS_all, BIS_ref) # Save reference u, reference y and sampled patients BIS to file
savehistdatatofile(BIS_all, α, n_pat_sample) # Save histogram data of final BIS (at 10 min) to file

num60, num40, totaloutside = countBISoutside(BIS_all)
percentoutside = totaloutside / n_pat_sample * 100

