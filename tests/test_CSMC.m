clear
rng(1)
addpath(genpath('csmc_cp'));
addpath(genpath('tests'));
addpath(genpath('fearnhead_cp'));
addpath(genpath('tests'));



%% test conditional stratified resampling
% come up with realistic test scenario
% set up weights and particles
N = 200;
SS_nm1 = 1:N+1;
[~,card_SS] = size(SS_nm1);
W_nm1 = abs(randn(1,N+1));
W_nm1 = W_nm1/sum(W_nm1);
logW = log(W_nm1);
[C,logC] = solveC(logW);

% get the unsafe particles
unsafe_mask = logW+logC<=0;
I_nm1 = SS_nm1(unsafe_mask);

% make sure tau is in unsafe set
tau_kappa = I_nm1(1);

L_nm1 = card_SS-sum(unsafe_mask);
W_for_I_nm1 = W_nm1(unsafe_mask);

% resample, making sure tau kappa is preserved
[ O_nm1 ] = conditionalStratifiedResampling( tau_kappa,I_nm1,W_for_I_nm1,N,L_nm1 );

assert(sum(O_nm1&(tau_kappa==I_nm1))==1)
disp('Passed test of tau persisting after resampling')
S_nm1 = zeros(1,N);
S_nm1(~unsafe_mask) = 1;
S_nm1(unsafe_mask) = O_nm1;



%% test CSMC

T = 200;
N = 50;
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
M = 1e2;
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
[~,tau_length] = size(tau_collect);
tau_burnin = tau_collect(floor(tau_length/2):tau_length);
histogram(tau_burnin(tau_burnin>0),T)







