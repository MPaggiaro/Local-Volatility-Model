function [fwd_vol] = fwd_from_spot_vol(T,spot_vol)

% T.. expiries of the market options (row vector of size N)
% spot_vol.. (M x N) matrix of spot volatilities
%
% returns the (M x N) matrix fwd_vol of forward volatilities

[rows, ~] = size(spot_vol);

fwd_vol = spot_vol;

var = spot_vol.^2 .* repmat(T,rows,1);
fwd_var = var(:,2:end) - var(:,1:end-1);
fwd_time = T(2:end) - T(1:end-1);

fwd_vol(:,2:end) = sqrt(fwd_var ./ repmat(fwd_time,rows,1));

end

