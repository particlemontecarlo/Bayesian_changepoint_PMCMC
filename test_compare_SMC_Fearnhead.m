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

%% compare Fearnhead and SMC with large number of particles
N = T;

% run forward filtering
[SS_all_Fearnhead,log_W_all_Fearnhead] = forwardFilteringFearnhead(params,T);
[SS_all,log_W_all] = changepoint_SMC(N,params);

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



%% compare SMC and Fearnhead with N<T
N = T/2;

% run forward filtering
[SS_all_Fearnhead,log_W_all_Fearnhead] = forwardFilteringFearnhead(params,T);
[SS_all,log_W_all] = changepoint_SMC(N,params);


M = 1e3;
tau_collect_Fearnhead = [];
tau_collect_SMC = [];

for m=1:M
    [ tau ] = bwdsSampling(params,SS_all_Fearnhead,log_W_all_Fearnhead);
    tau_collect_Fearnhead = [tau_collect_Fearnhead tau];
    [ tau ] = bwdsSampling(params,SS_all,log_W_all);
    tau_collect_SMC = [tau_collect_SMC tau];
end

figure(2)
subplot(2,1,1)
histogram(tau_collect_Fearnhead(tau_collect_Fearnhead>0),0:T)

subplot(2,1,2)
histogram(tau_collect_SMC(tau_collect_SMC>0),0:T)



