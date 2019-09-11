### PMCMC for Bayesian changepoint inference
The following code comprises methodology to perform Bayesian inference of multiple structural breaks in time series.
The code primarily implements routines provided in [Nick Whiteley, Christophe Andrieu, and Arnaud Doucet. "Bayesian computational methods for inference in multiple change-points models." (2011)] building on routines developed in [Paul Fearnhead, and Zhen Liu. "On‐line inference for multiple changepoint problems." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 69.4 (2007): 589-605.].

The two scripts provided provide examples of Bayesian inference in changepoint problems with 'example_MCMC.m' implementing an expensive but exact Gibbs sampler iteratively sampling from the full conditional over changepoints and then the full conditional of the parameter defining the a priori distribution over changepoint locations (in this case we assume that changepoints occur independently within the time series with a certain constant probability).
The second script 'example_particle_Gibbs.m' performs exactly the same inference though with particle approximations introduced as per the two aforementioned papers. In particular the conditional stratified optimal resampling routine is implemented, building on the highly efficient stratified optimal resampling of [Fearnhead 2007].



