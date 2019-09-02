%%%%Stratified resampling for changepoints
function logWnew = conditionalFiltering(logW,S,Suppn,N,n,tau_star)

% % % description of variables
% w     : T x min(N+1,n)    weights of random support
% wbar  : 1 x min(N+1,n)    weights of random support afte resampling
% SS    : 1 x min(N+1,n)    random support vector
% S     : 1 x min(N+1,n)    indicator array
% In    : 1 x min(N+1,n)    array used in resampling


% we have that SS captures where a possible changepoint is with SSn at most 
% roughly of length n and that Sn is an indicator array of size n, 
% suggesting whether the particle is alive or dead
% for each particle indexed by elements in SSn we have a weight


% initialise arrays
logW = zeros(T,N+1);
logWbar = zeros(T,N+1);
SS = zeros(T,N+1);
S = zeros(T,N+1);


% n=1
S(1,1) = 1;


% iterate forward
for n=2:n

    % find the constant used in resampling
    if n-1<N
        logC = Inf;
    else
        [~,logC] = solveC(logW);
    end

    % find the support points that are to haev weight 0 by resampling
    SS_nm1 = SS(n-1,:);
    logW_nm1 = logW(n-1,:);
    killed_particle_mask = logW_nm1<-logC;
    Lnm1 = sum(killed_particle_mask);
    if Lnm1>0
        % resample
        for i=1:(N-Lnm1)
            % provide the location and weights of the particles
            [ O_nm1 ] = stratifiedResampling( I_nm1,logW_nm1 )

        end
    end

    % propagate particles deterministically
    % as an input we should have 
    % 1) Data at time n
    % 2) Weights of particles at time n
    % 3) Value of constant
    % 
    %
    % Should return 
    % 1) new weights of particles after propagation
    data_n = 0;
    [logWbar,predData,predStatic] = propagateParticles(data_n,SS_nm1,logW_nm1,logCnm1,n);
    
    % update the support
    SSn = SS
    
    for i=1:card_SS
        xnm1 = SS_nm1(i);
        summand_arr(i) = fn(n-1,xnm1) * Wnm1(xnm1);
    end
    
    
    % update temporary variables
    
end


end


function [logWbar,predData,predStatic] = propagateParticles(data_n,SS_nm1,logW_nm1,logCnm1,n)
%%% this function is used to determinstically update the particles

% get the number of active particles
[card_SS,~] = size(SS_nm1);
warning('size of SS must be defined')


% get value of last Wbar
summand_arr = zeros(1,card_SS);
for i=1:card_SS
    xnm1 = SS_nm1(i);
    summand_arr(i) = fn(n-1,xnm1) * Wnm1(xnm1);
end
Wbar_n_nm1 = gn(n-1)*sum(summand_arr);


% get values of other Wbars
Wbar_arr = zeros(1,card_SS);
for i=1:card_SS
    xnm1 = SS_nm1(i);
    Wbar_arr(xnm1) = gn(xnm1)*fn(xnm1,xnm1)*Wnm1(xnm1)...
        /min(1,Wnm1(xnm1)*exp(logCnm1));
end

end



