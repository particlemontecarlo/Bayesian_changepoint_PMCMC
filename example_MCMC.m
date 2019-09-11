%% test algorithm with simple changepoint problem - varying the means
% compare between Fearnhead and special case of particle filter
clear
rng(1)
addpath(genpath('smc_cp'));
addpath(genpath('fearnhead_cp'));
addpath(genpath('tests'));




T = 200;
taustar = [0,50,100,150];

s2 = 1;
means = [0,5,10,5];
obs_arr = zeros(1,T);
% generate data
for i = 1:T
    segment_current = sum(i>taustar);
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


%% compare SMC and Fearnhead with N<T
M = 1e4;
tau_collect = [];
pGeo_collect = zeros(1,M);
a_prior = 1;
b_prior = 1;
for m=1:M
    [SS_all_Fearnhead,log_W_all_Fearnhead,~] = forwardFilteringFearnhead(params,T);
    [ tau ] = bwdsSampling(params,SS_all_Fearnhead,log_W_all_Fearnhead);
    
    n_success = length(tau)-1;
    params.pGeo = betarnd(a_prior + n_success,b_prior + T - n_success);
    tau_collect = [tau_collect tau];
    pGeo_collect(m) = params.pGeo;
end

figure(2)
subplot(2,1,1)
histogram(tau_collect(tau_collect>0),0:T)
title('Posterior distribution of changepoints for Fearnhead')

subplot(2,1,2)
histogram(pGeo_collect(pGeo_collect~=0))
title('Posterior distribution of p for Fearnhead')







