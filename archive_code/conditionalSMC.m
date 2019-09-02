function [ tau,loglik ] = conditionalSMC( X,z,M,theta,tau,GDPPrior,n_quadrature_points )
%CONDITIONALSMC Performs conditional SMC sampling of changepoints,
%conditioned on an existing changepoint configuration
%   Method follows that of Whiteley (2010) and Fearnhead (2007). Method
%   assumes there exists a changepoint at t=0, for a series indexed between
%   1:t

if tau(1)~=0
    error('Invalid conditioning changepoint')
end


[~,T] = size(X);

%Allocate matrices for particle support and particle weights
[S] = zeros(T,M);
logW  = repmat(-Inf,T,M);
%Define start support
S(1,1) = 0;

%Allocate weight 1 to W1(0)
logW(1,1) = 0;
N=1;
loglik = zeros(1,T);


[ logp,predData,predStatic ] = predictiveStart( X(:,1),z(1),GDPPrior,n_quadrature_points );
loglik(1) = logp;

%Iterate through time, evaluating particle weights
%h = waitbar(0);
for t=2:T
    %waitbar(t/T,h,'Conditional SMC progress: ')
    
    if t<=M
        logWtm1 = logW(t-1,1:N);
        St = [S(t-1,1:N) t-1];
        S(t,1:N+1) = St;
        N = N+1;
    end
    
    
    %Evaluate new weights of particles
    [logWt,predData,predStatic] = propagateParticles(t,N,logWtm1,X(:,t),z(t),...
        predData,predStatic,theta,GDPPrior,n_quadrature_points);
    
    %Normalise new weights
    logWn = logWt-max(logWt);
    logW(t,1:N) = logWn - log(sum(exp(logWn)));
    loglik(t) = log(sum(exp(logWn))) + max(logWt);
    
    %Resample if N>M
    if t>= M && t<T
        recentTau = tau(tau<=t);
        kappa = recentTau(end);
        %If kappa == t, conditional particle will necessarily be created, 
        %otherwise will need to be preserved in resampling
        if kappa==t
            kappaIndex = [];
            kappaIst = true;
        else
            kappaIndex = find(S(t,:)==kappa,1);
            if isempty(kappaIndex)
                error('Cannot find conditioned particle in support')
            end
            kappaIst = false;
        end


        logweights = conditionalStratifiedResampling(logW(t,1:N),M-1,kappaIndex,kappaIst);

        pSurvived = logweights~=-Inf;
        N = sum(pSurvived)+1;
        
        S(t+1,1:N) = [S(t,pSurvived) t];
        
        A = size(S);
        if (A(2)>M)
            error('Number of particles exceeds M')
        end
        
        logWtm1 = logweights(pSurvived);
        predData = predData(pSurvived);

    end
    
    
end



loglik = cumsum(loglik);

%Backwards sampling
taui = S(T,discretesample(exp(logW(T,:)),1));

tau = taui;
while true
    if taui<=1
        if taui~=0
            tau = [0 tau];
        end
        break
    end
    
    logWtau = logW(taui,:);
    Stau = S(taui,:);
    
    logWsupport = logWtau(logWtau>-Inf);
    Ssupport = Stau(logWtau>-Inf);
    nsupport = length(Ssupport);
    
    logpsupport = zeros(1,nsupport);
    for i = 1:nsupport
        logpsupport(i) = logWsupport(i) + log(transition(taui,Ssupport(i),theta.pGeo));
    end
    
    logpn = logpsupport-max(logpsupport);
    taui = Ssupport(discretesample(exp(logpn)/sum(exp(logpn)),1));
    tau = [taui tau];

end



end


%Function evaluating weights based on posterior predictive and transition
%parobabilities
function [logWt,predDyn,predStatic] = propagateParticles(t,N,logWtm1,Xt,zt,...
    predDyn,predStatic,theta,GDPPrior,n_quadrature_points)



pGeo = theta.pGeo;
logWt = zeros(1,N);

pre_calc = pre_calculation(predStatic,Xt,zt);
for j = 1:(N-1)
    [logp,predDyn(j)] = predictive(Xt,zt,predDyn(j),predStatic,pre_calc );
    logWt(j) = logp + log(1-pGeo) + logWtm1(j);
end

% plot(exp(predDyn(1).logfProd(:,1)))

tmp = zeros(1,N-1);
for i = 1:(N-1)
    tmp(i) = log(pGeo) + logWtm1(i);
end

tmp_m = max(tmp);
[logp,predDyn(N),predStatic] = predictiveStart(Xt,zt,GDPPrior,n_quadrature_points);
logWt(N) = logp + log(sum(exp(tmp-tmp_m))) + tmp_m;

end

% this can be used as the additional density evaluations need only be
% evaluated for one particle
function [pre_calc] = pre_calculation(predStatic,Xt,zt)

rhoGrid = predStatic.rhoGrid;
rhoVar = predStatic.rhoVar;

logf1 = logpXcondRhoZ(rhoGrid,Xt',zt,rhoVar);
pre_calc = logf1;
end


%Recursively evaluates the posterior predictive
function [logp,predDynj] = predictive(~,~,predDynj,predStatic,pre_calc)


logfProdOld = predDynj.logfProd;

% rhoGrid = predStatic.rhoGrid;
h = predStatic.h;
logrhoInt = predStatic.logrhoInt;
% rhoVar = predStatic.rhoVar;
logdenom = predDynj.logdenom;


%logf1 = logpXcondRhoZ(rhoGrid,Xt',zt,rhoVar);
logf1 = pre_calc;
logfProdNew = logf1 + logfProdOld;
a1 = logfProdNew + logrhoInt - logdenom;
max_a1 = max(a1,[],1);
diff_mat = bsxfun(@minus,a1,max_a1);
logpAll = log(h/3) + log(sum(exp(diff_mat ))) + max_a1 ;

predDynj.logfProd = logfProdNew;
predDynj.logdenom = logpAll + logdenom;

if any(isinf(logpAll))
    warning('Underflow in prediction')
    logpAll = logpAll(~isinf(logpAll));
end


logp = sum(logpAll);

end


function [ logp,predDynj,predStatic ] = predictiveStart( Xt,zt,GDPPrior,n_quadrature_points )
%POSTPREDICTIVESTART 

thresh = 1e-10;
a = -1+thresh;
b = 1-thresh;


h = (b-a)/(n_quadrature_points+1);
SimpsonFactor =  [1 repmat([4 2],1,n_quadrature_points/2 - 1) 4 1]';
rhoGrid = linspace(a,b,n_quadrature_points+1)';

if GDPPrior
    rhoPrior = rhoPriorGDP(rhoGrid);
else
    rhoPrior = rhoPriorNormal(rhoGrid);
end
logrhoInt = log(rhoPrior)+log(SimpsonFactor);
rhoVar = (1-rhoGrid.^2);


logf1 = logpXcondRhoZ(rhoGrid,Xt',zt,rhoVar);
pnew = h/3 * sum(exp(logf1+logrhoInt),1);
logfProd = logf1;
logdenom = log(pnew);


logp = sum(logdenom);

predDynj.logfProd = logfProd;
predDynj.logdenom = logdenom;

predStatic.rhoGrid = rhoGrid;
predStatic.logrhoInt = logrhoInt;
predStatic.rhoVar = rhoVar;
predStatic.h = h;

end

%%%Transition probabilities
function p = transition(tau,Ct,pGeo)

if Ct ==(tau-1)
    p = pGeo;%(G(tau-Ct-1,pGeo) - G(tau-2-Ct,pGeo)) / ( 1 - G(tau-2-Ct,pGeo) );
else
    p = (1-pGeo);%(1 - G(tau-Ct-1,pGeo)) / ( 1 - G(tau-2-Ct,pGeo) );
end

end


function logp = logpXcondRhoZ(Rho,Xti,zt,rhoVar)

mu = Rho*zt;
logp = - 0.5*log(2*pi*rhoVar) - 0.5*((Xti-mu).^2)./rhoVar ;

end



%Generalised double Pareto density
function rhoPrior = rhoPriorGDP(Rho)

rhoPrior = GDPdensity(Rho./(1-Rho.^2).^0.5,3,1) .* abs((1 - Rho.^2).^(-3/2));
end

%Generalised double Pareto density
function rhoPrior = rhoPriorNormal(Rho)

rhoPrior = normpdf(Rho./(1-Rho.^2).^0.5) .* abs((1 - Rho.^2).^(-3/2));
end


%Generalised Pareto density for first component
function rhoPrior = rhoPriorGP(Rho)

rhoPrior = GPdensity(Rho./(1-Rho.^2).^0.5,3,1) .* abs((1 - Rho.^2).^(-3/2));
end




function [ logp ] = G( s,pGeo )
%G Distribution function of distance between two successive changepoints

p = 1 - (1-pGeo)^s;
end





