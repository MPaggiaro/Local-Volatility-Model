function [eta] = localvol(T,K,V,qt,qk)

% returns the local volatility function eta(qt,qk) given the values V of the function at the nodes (T,K)
% the interpolation is piecewise constant in time and spline with flat extrapolation in strike

idx = find(T>=qt,1);
k = K(:,idx);
v = V(:,idx);

eta = interp_flat_extrap(k,v,qk,'spline');

%figure('Name',['Local Volatility funtion at time ' num2str(T(idx),'%f')]);
%plot(qk,eta,'s',k,v,':.');
end

