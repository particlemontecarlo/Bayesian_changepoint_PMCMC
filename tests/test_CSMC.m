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

%% run test for large number of particles comparing with exact algorithm
N = T;

[SS_all,log_W_all,log_W_bar_all] = forwardFilteringCSMC(N,tau,params);
[SS_all_Fearnhead,log_W_all_Fearnhead,~] = forwardFilteringFearnhead(params,T);

% test support is the same
assert(all(all(SS_all_Fearnhead==SS_all)))

% test that the differences in recursions are the same up to machine
% precision
diff_array = zeros(T,T);
log_W_all_1 = log_W_all;
log_W_all_2 = log_W_all_Fearnhead;
log_W_all_1(log_W_all_1==-Inf) = 0;
log_W_all_2(log_W_all_2==-Inf) = 0;
diff_measure = abs(log_W_all_1-log_W_all_2)<1e-8;
if all(all(diff_measure))
    disp('Pass test comparing filtering recursions for large number of particles')
else
    error('Difference in SMC and Fearnhead for large number particles')
end


%% run test for small number of particles using iterated conditional SMC

N = T/10;
M = 1e3;
[SS_all_Fearnhead,log_W_all_Fearnhead,~] = forwardFilteringFearnhead(params,T);
tau_collect_Fearnhead = [];
for m=1:M
    [ tau_Fearnhead ] = bwdsSampling(params,SS_all_Fearnhead,log_W_all_Fearnhead);
    tau_collect_Fearnhead = [tau_collect_Fearnhead tau_Fearnhead];
end


tau_collect = [tau];
like_ests = zeros(1,M);
for m=1:M
    m
    [SS_all,log_W_all,log_W_bar_all] = forwardFilteringCSMC(N,tau,params);
    [ tau ] = bwdsSampling(params,SS_all,log_W_all);
    [like_est] = get_likelihood_est(log_W_bar_all);
    like_ests(m) = like_est;
    tau_collect = [tau_collect,tau];
end


figure(2);
subplot(2,1,1)
histogram(tau_collect_Fearnhead(tau_collect_Fearnhead>0),0:T)
title('Posterior distribution of changepoints for Fearnhead')

subplot(2,1,2)
[~,tau_length] = size(tau_collect);
tau_burnin = tau_collect(floor(tau_length/2):tau_length);
histogram(tau_burnin(tau_burnin>0),0:T)
title('Posterior distribution of changepoints iterated cSMC')







