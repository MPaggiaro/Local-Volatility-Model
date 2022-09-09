function [vol] = model_volatility_mod(T,K_norm,V,N,M,K_min,K_max,Scheme)
% Compute model implied volatility of a local volatility model
% T.. LV expiries, K_norm.. LV strikes, V.. LV matrix
% N,M,K_min,K_max,Scheme.. settings for the Dupire solver
%
% vol(i,j).. model implied volatility at time T(i) and strike K_norm(i,j) 

[rows, cols] = size(K_norm);
vol = zeros(rows,cols);
%the whole vector T is now passed to the function
[ k, C ] = solve_dupire_mod( T, K_norm, V, N, M, K_min, K_max, Scheme);
%hence, C is a matrix
for i=1:length(T)
    price = interp1(k,C(:,i),K_norm(:,i));
    vol(:,i) = blsimpv(1,K_norm(:,i),0,T(i),price);
end

end
