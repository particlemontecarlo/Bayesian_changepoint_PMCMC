function [tau] = changepoint_SMC_old(N,params)


[~,T] = size(params.Y);

n = 2;
SS_nm1 = [0];
W_nm1 = [1];


% run forward filtering
SS_all = zeros(T,T);
SS_all(1,1) = SS_nm1;
W_all = zeros(T,T);
W_all(1,1) = W_nm1;
for n=2:T
    [SS_updated,Wbar_updated] = forwardFiltering(n,N,SS_nm1,W_nm1,params);
    W_nm1 = Wbar_updated/sum(Wbar_updated);
    W_nm1 = Wbar_updated/sum(Wbar_updated);
    SS_nm1 = SS_updated;
    SS_all(n,1:min(n,N+1)) = SS_updated;
    W_all(n,1:min(n,N+1)) = W_nm1;
end

% backward sampling
taui = SS_all(T,find(  mnrnd(1,W_all(T,:))   ));
tau = taui;
while taui>1

    Wtau = W_all(taui,:);
    Stau = SS_all(taui,:);
    
    logWsupport = log(Wtau(Wtau>0));
    Ssupport = Stau(Wtau>0);
    nsupport = length(Ssupport);
    
    logpsupport = zeros(1,nsupport);
    for i = 1:nsupport
        logpsupport(i) = logWsupport(i) + log(fn(taui,Ssupport(i),params));
    end
    logpn = logpsupport-max(logpsupport);
    taui = Ssupport(  find(mnrnd(1,exp(logpn)/sum(exp(logpn))))  );
    tau = [taui tau];
end
if taui~=0
    tau = [0 tau];
end


end
