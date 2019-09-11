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







