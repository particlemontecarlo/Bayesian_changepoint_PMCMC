%% test algorithm with simple changepoint problem - varying the means
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


%% sample tau given pGeo
params.pGeo = 0.02;
params.Y = obs_arr;
params.sigma02 = 10;
params.sigma2 = s2;
[~,T] = size(params.Y);


N = T;

R = 1000;
a_prior = 1;
b_prior = 1;

tau_collect = [];
pGeo_collect = zeros(1,R);
[SS_all,log_W_all] = changepoint_SMC(N,params);

for r=1:R
    disp(r)
    [ tau ] = bwdsSampling(params,SS_all,log_W_all);
    tau_collect = [tau_collect tau];
    %n_success = length(tau)-1;
    %params.pGeo = betarnd(a_prior + n_success,b_prior + T - n_success);
    %pGeo_collect(r) = params.pGeo;
end

figure(2)
histogram(tau_collect(tau_collect>0),0:T)


%% run with particle number less than the number of observations
N = 10;

R = 1000;
a_prior = 1;
b_prior = 1;

tau_collect = [];
pGeo_collect = zeros(1,R);
[SS_all,log_W_all] = changepoint_SMC(N,params);

for r=1:R
    disp(r)
    [ tau ] = bwdsSampling(params,SS_all,log_W_all);
    tau_collect = [tau_collect tau];
    %n_success = length(tau)-1;
    %params.pGeo = betarnd(a_prior + n_success,b_prior + T - n_success);
    %pGeo_collect(r) = params.pGeo;
end

figure(5)
histogram(tau_collect(tau_collect>0),0:T)

