function [ k, C ] = solve_dupire_mod( T, K_norm, V, N, M, K_min, K_max, Scheme)

%the modified function solves numerically the Dupire's PDE 
%using the requested theta-method.

h = linspace(log(K_min),log(K_max),M);
C=zeros(M,length(T));
k=exp(h);

%Anyway, we need some boundary conditions:
%the prices at time 0 are respectively equal to their strikes;
%the price at any time for the strike k_max (last row) is always 0. So

I=eye(M);

for j=1:length(T)
    
    if j==1
        C(:,j) = max(0,1-k)';
        t = linspace(0,T(j),N);
    else
        C(:,j) = C(:,j-1);
        t = linspace(T(j-1),T(j),N);
    end
    
    dt=t(2)-t(1);
    
    if strcmp(Scheme,'f')
        % forward scheme
        for i=1:length(t)
            A = build_A(T,K_norm,V,t(i),h);
            C(:,j)=(I-dt*A)*C(:,j);
        end
    elseif strcmp(Scheme,'b')
        % backward scheme
        for i=1:length(t)-1
            A = build_A(T,K_norm,V,t(i+1),h);
            C(:,j)=(I+dt*A)\C(:,j);
        end
    elseif strcmp(Scheme,'cn')
        % Crank-Nicolson scheme
        for i=1:length(t)-1
            A = build_A(T,K_norm,V,t(i),h);
            Ap = build_A(T,K_norm,V,t(i+1),h);
            C(:,j)=(I+0.5*dt*Ap)\((I-0.5*dt*A)*C(:,j));
        end
    end
end