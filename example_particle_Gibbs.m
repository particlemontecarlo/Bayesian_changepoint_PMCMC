clear
rng(1)
addpath(genpath('csmc_cp'));
addpath(genpath('tests'));
addpath(genpath('fearnhead_cp'));
addpath(genpath('tests'));



%% test CSMC

T = 200;
tau_star = [0,50,100,150];

s2 = 1;
means = [0,5,10,5];
obs_arr = zeros(1,T);
% generate data
for i = 1:T
    segment_current = sum(i>tau_star);
    obs_arr(i) = randn(1)*sqrt(s2) + means(segment_current);
end

figure(1)
plot(obs_arr)

% set parameters
params.pGeo = 0.02;
params.Y = obs_arr;
params.sigma02 = 20;
params.sigma2 = s2;
[~,T] = size(params.Y);
tau0 = [0:10:(T-1)];
tau = tau0;



%% run test for small number of particles using iterated conditional SMC

N = 10;
M = 1e3;

tau_collect = [tau];
like_ests = zeros(1,M);
pGeo_collect = zeros(1,M);
a_prior = 1;
b_prior = 1;
for m=1:M
    disp(m)
    % cSMC to sample kernel leaving the full conditional over changepoints
    % invariant
    [SS_all,log_W_all,log_W_bar_all] = forwardFilteringCSMC(N,tau,params);
    [ tau ] = bwdsSampling(params,SS_all,log_W_all);
    [like_est] = get_likelihood_est(log_W_bar_all);
    
    % update changepoint parameters conditional on changepoints
    n_success = length(tau)-1;
    params.pGeo = betarnd(a_prior + n_success,b_prior + T - n_success);
    
    % save results
    like_ests(m) = like_est;
    tau_collect = [tau_collect,tau];
    pGeo_collect(m) = params.pGeo;
end





figure(2);
subplot(2,1,1)
histogram(tau_collect(tau_collect>0),0:T)
title('Posterior distribution of changepoints particle Gibbs')

subplot(2,1,2)
histogram(pGeo_collect)
title('Posterior distribution of geometric parameter for changepoint prior')



