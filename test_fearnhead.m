%%%% test algorithm with simple changepoint problem - varying the means
clear
rng(1)
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

figure(3)
plot(obs_arr)


%%
params.pGeo = 0.02;
params.Y = obs_arr;
params.sigma02 = 10;
params.sigma2 = s2;
[~,T] = size(params.Y);

% run forward filtering
[SS_all,log_W_all] = forwardFilteringFearnhead(params,T);

M = 1e3;
tau_collect = [];
for m=1:M
    [ tau ] = bwdsSampling(params,SS_all,log_W_all);
    tau_collect = [tau_collect tau];
end

figure(4)
histogram(tau_collect(tau_collect>0),0:T)


%% is it possible to compare likelihood estimates for the two methods???






