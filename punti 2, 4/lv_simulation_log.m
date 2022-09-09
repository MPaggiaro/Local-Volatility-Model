function[S] = lv_simulation_log(T,Fwd,V,K,N,M,expiry)
		% Monte Carlo simulation for the asset S under the LV model
        % T.. LV expiries, F.. forwards at T, V.. LV matrix, K.. LV strikes
		% M.. MC timesteps, N.. MC simulations, expiry.. MC time horizon
		
		logX = zeros(N,1);
        
        % normalized market strikes
        [rows, ~] = size(K);
        K_norm = K ./ repmat(Fwd, rows, 1);
        
        W = randn(N,M,length(expiry));

        for k = 1:length(expiry)
            if k==1
                dt = expiry(k)/M;
                t = (0:M)'*dt;
            else
                dt = (expiry(k)-expiry(k-1))/M;
                t = expiry(k-1) + (1:M)'*dt;
            end
            
            for j = 1:M
                sigma = localvol(T,K_norm,V,t(j),exp(logX(:)));
                logX(:) = logX(:) - 0.5 * dt * sigma.^2 + sigma .* W(:,j,k) * sqrt(dt);
            end

            fwd = interp1(T,Fwd,expiry(k));
            S(k,:) = fwd * exp(logX(:));
        end
	end

