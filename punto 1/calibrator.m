function [V, ModelVol, MaxErr] = calibrator(T,K_norm,MktVol,threshold,MaxIter,N,M,K_min,K_max,Scheme)
% Calibrate a local volatility matrix for the normalized asset X
%
% T.. expiries of the market options 
% K_norm.. normalized strikes of the market options 
% MktVol.. implied volatilities of the market options 
% threshold,MaxIter.. calibration threshold and max nb of iterations 
% N,M.. grid dimension in time and strike (Dupire)
% K_min,K_max.. min and max strikes in the discretized grid (Dupire) 
% scheme.. finite difference scheme 'f','b','cn' (Dupire)
%
% V is the calibrated local volatility matrix
% ModelVol is the model implied volatilities at the normalized strikes
% MaxErr is the vector of calibration error at each iteration

nb_iter = 1;
MaxErr(nb_iter) = 100;

% initial guess for the LV matrix
V = MktVol;

% forward market volatility
mkt_fwd_vol = fwd_from_spot_vol(T,MktVol); %slide 107

% calibration progress
h = waitbar(0,'Calibration...');
    
while ( MaxErr(nb_iter) > threshold && nb_iter < MaxIter )
        
    % update iteration number
    nb_iter = nb_iter+1;
    
    % compute model implied volatilities
    ModelVol = model_volatility(T,K_norm,V,N,M,K_min,K_max,Scheme);
    
    % compute max_err
    MaxErr(nb_iter) = max(max(abs(ModelVol-MktVol)));
    waitbar(threshold/MaxErr(nb_iter));
    
    % compute model fwd vol 
    model_fwd_vol = fwd_from_spot_vol(T,ModelVol);
    
    % compute new LV parameters
    V = V ./ model_fwd_vol .* mkt_fwd_vol;

end

close(h);
MaxErr = MaxErr(2:nb_iter);

end

