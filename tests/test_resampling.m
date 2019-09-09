clear
rng(1)
%%
% resampling should ensure the conditioning particle survives while
% approximating the optimal resampling
% should have
% 1. random support
% 2. weights of particlesrng(1)
addpath(genpath('smc_cp'));
addpath(genpath('tests'));



%%
% test solution to resampling scheme
% one test in a simple case and one test in a random case
N = 1;
W = [0.25,0.75]; % solution should be 2
logW = log(W);
[Cdev,logCdev] = solveC(logW);
Ntest = sum(min(1,W*Cdev));
assert(abs(Ntest-N)<=1e-6)



N = 1e2;
W = (rand(1,N+1));
W = W/sum(W);
W = W/sum(W);
logW = log(W);
tic;
[C,logC] = solveC(logW);
toc;
Ntest = sum(min(1,W*C));
assert(abs(Ntest-N)<=1e-8)

%%
logWu = [-685.2660 -641.6245 -570.9773 -354.5661 -210.1526  -85.8920   -0.3898   -1.2358   -3.4343];
Wu  = exp(logWu);
W = Wu/sum(Wu);
logW = log(W);
N = 10;
[C,logC] = solveC(logW);
logW+logC<=0



%% test stratified resampling
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
L_nm1 = card_SS-sum(unsafe_mask);
W_for_I_nm1 = W_nm1(unsafe_mask);

[ O_nm1 ] = stratifiedResampling( I_nm1,W_for_I_nm1,N,L_nm1 );
S_nm1 = zeros(1,N);
S_nm1(~unsafe_mask) = 1;
S_nm1(unsafe_mask) = O_nm1;



%% test forward filtering
SS_nm1 = [0,1,4,5];
log_W_nm1 = log([0.1,0.1,0.4,0.4]);
n = 7;
N = 3;
params.pGeo = 0.1;
params.Y = randn(1,n);
params.sigma02 = 1;
params.sigma2 = 1;


[SS_updated,Wbar_updated] = forwardFilteringSMCRecursion(n,N,SS_nm1,log_W_nm1,params)

%%
SS_nm1 = [0];
log_W_nm1 = log([1]);
n = 2;
N = 3;
params.pGeo = 0.1;
params.Y = randn(1,n);
params.sigma02 = 1;
params.sigma2 = 1;
[SS_updated,Wbar_updated] = forwardFilteringSMCRecursion(n,N,SS_nm1,log_W_nm1,params)


