%% test algorithm with simple changepoint problem - varying the means
% compare between Fearnhead and special case of particle filter
clear
rng(1)
addpath(genpath('smc_algorithm'));
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
params.sigma02 = 10;
params.sigma2 = s2;
[~,T] = size(params.Y);
N = T;


% run forward filtering
[SS_all_Fearnhead,log_W_all_Fearnhead] = forwardFilteringFearnhead(params,T);
[SS_all,log_W_all] = changepoint_SMC(N,params);

% test support is the same
assert(all(all(SS_all_Fearnhead==SS_all)))

% test that the differences in recursions are the same up to machine
% precision
diff_array = zeros(T,T);
log_W_all(log_W_all==-Inf) = 0;
log_W_all_Fearnhead(log_W_all_Fearnhead==-Inf) = 0;
diff_measure = abs(log_W_all-log_W_all_Fearnhead)<1e-8;
assert(all(all(diff_measure)))








