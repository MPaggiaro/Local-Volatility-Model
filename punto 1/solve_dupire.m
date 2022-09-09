function [ k, C ] = solve_dupire( T, K, V, Expiry, N, M, K_min, K_max, Scheme)
% Solve Dupire equation for the asset X up to Expiry
%
% T.. LV expiries, K.. LV nodes, V.. LV matrix
% Expiry.. expiry of the call options 
% N,M.. grid dimension in time and strike
% K_min,K_max.. min and max strikes in the discretized grid 
% scheme.. finite difference scheme 'f','b','cn'
%
% C(i) is the price of a call option on X at strike k(i)=exp(h(i))

% create time and strike grid
t = linspace(0,Expiry,N);
h = linspace(log(K_min),log(K_max),M);
k = exp(h);

dt = t(2)-t(1);
I = eye(M);
C = max(0,1-exp(h))';
%C_0 = C;

if strcmp(Scheme,'f')
    % forward scheme 
    for i=1:length(t)
    	A = build_A(T,K,V,t(i),h);
        C = (I-dt*A)*C;
    end
elseif strcmp(Scheme,'b')
    % backward scheme 
    for i=1:length(t)-1
    	A = build_A(T,K,V,t(i+1),h);
        C = (I+dt*A)\C;
    end    
elseif strcmp(Scheme,'cn')
    % Crank-Nicolson scheme 
    for i=1:length(t)-1
    	A = build_A(T,K,V,t(i),h);
    	Ap = build_A(T,K,V,t(i+1),h);
        C = (I+0.5*dt*Ap)\( (I-0.5*dt*A)*C );
    end
end

end

