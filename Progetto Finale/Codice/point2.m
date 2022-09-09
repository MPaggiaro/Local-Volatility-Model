%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  SCRIPT PER IL TITOLO ECORP  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% market expiries
 T = [0.140  0.216 	 0.312 	 0.466 	 0.734 	 0.984 	 1.060 	 1.482 	 1.981 	 2.058 	 2.978 ];
% forwards at market expiries
Fwd = [2694.93	2696.09	2699.20	2698.08	2704.40	2709.25	2712.81	2720.39	2730.30	2732.90	2763.60];

% discount factors at market expiry
D =[ 0.9953772 	 0.9940096 	 0.9919588 	 0.9887443 	 0.9827208 	 0.9770636 	 0.9752305 	 0.9649902 	 0.9529284 	 0.9509734 	 0.9266855 
];

% market strikes
K = [
    2425	2350	2225	2100	1925	1775	1750	1600	1475	1450	1275
    2625	2600	2575	2525	2475	2425	2425	2350	2300	2275	2200
    2675	2675	2650	2650	2625	2625	2625	2600	2575	2575	2575
    2700	2700	2700	2700	2700	2700	2725	2725	2725	2725	2775
    2725	2725	2725	2750	2775	2800	2800	2825	2875	2875	2950
    2750	2775	2800	2825	2875	2925	2925	3000	3075	3075	3250
    2825	2850	2900	2975	3075	3150	3175	3300	3450	3475	3600
    ];

% LV matrix
MktVol = [
    0.1209	0.1255	0.1367	0.1460	0.1594	0.1706	0.1723	0.1815	0.1882	0.1895	0.1975
    0.0906	0.0944	0.0988	0.1084	0.1200	0.1284	0.1295	0.1395	0.1466	0.1481	0.1576
    0.0843	0.0855	0.0918	0.0983	0.1096	0.1164	0.1178	0.1269	0.1346	0.1352	0.1448
    0.0803	0.0825	0.0869	0.0943	0.1052	0.1123	0.1121	0.1210	0.1285	0.1292	0.1386
    0.0795	0.0810	0.0850	0.0906	0.1007	0.1069	0.1084	0.1164	0.1226	0.1234	0.1336
    0.0758	0.0779	0.0796	0.0866	0.0955	0.1007	0.1022	0.1092	0.1154	0.1163	0.1255
    0.0768	0.0738	0.0764	0.0810	0.0876	0.0923	0.0925	0.0983	0.1041	0.1042	0.1177
 ];


% spot price
S0=2687.20;

% normalized market strikes
[rows, cols] = size(K);
Knorm = K ./ repmat(Fwd, rows, 1);

%%
%%%%%%% 2.1 MODEL CALIBRATION %%%%%%

% We use the optimized solver implemented in point1.

% Dupire solver setting
N = 10;
M = 200;
Kmin = 0.1;
Kmax = 3;

% calibration settings
tol = 0.001;
MaxIter = 100;

% calibration
[V, ModelVol, MaxErr] = calibrator_mod(T,Knorm,MktVol,tol,MaxIter, N,M,Kmin, Kmax, 'cn');

% plot calibration error at each iteration
figure (1)
grid on
plot(MaxErr);
title('Trend of the approximation error')

%%
%%%%%% 2.2 PRICING OF TWO CALL OPTIONS %%%%%%%

% option data
expiry = 0.5;
strike = [0.9*S0,1.1*S0];

%%% PRICING CALLS via DUPIRE EQUATION %%%
[k,C] = solve_dupire(T,Knorm,V,expiry,N,M,Kmin,Kmax,'cn');

% discount factor and forward at the expiry date:
discount_factor = interp1(T,D,expiry);
fwd = interp1(T,Fwd,expiry);

norm_strike = strike/fwd;
c=interp1(k,C,norm_strike);

% Price trasformation:
call_prices_dup = fwd * discount_factor * c;

% Implied volatility:
impl_vol_dup = blsimpv(1,norm_strike,0,expiry,c);

%%% PRICING CALLS via MC %%%
N = 1000000;
M = 100;

% Running simulation:
S = lv_simulation_log(T,Fwd,V,K,N,M,expiry);

% Pricing call options:
call_prices_MC(1) = discount_factor*mean(max(S(1,:) - strike(1),0));
call_prices_MC(2) = discount_factor*mean(max(S(1,:) - strike(2),0));

% Implied volatility:
impl_vol_MC(1) = blsimpv(fwd,strike(1),0,expiry,call_prices_MC(1)/discount_factor);
impl_vol_MC(2) = blsimpv(fwd,strike(2),0,expiry,call_prices_MC(2)/discount_factor);

%%
%%%%%%%% 2.3 Monte Carlo simulation %%%%%%%%

expiry = [2 2.5];
strike = [0.9, 1.1];

discount_factor = interp1(T,D,expiry(2));

N=100000; %MC simulations
M=100; %timesteps

% MC simulation
S = lv_simulation_log(T,Fwd,V,K,N,M,expiry);

% option price
option_price = [discount_factor*mean(max(S(2,:) - strike(1)*S(1,:),0));
     discount_factor*mean(max(S(2,:) - strike(2)*S(1,:),0))]

fwd(1) = interp1(T,Fwd,expiry(1));
fwd(2) = interp1(T,Fwd,expiry(2));
impl_vol_fwd = [blsimpv(fwd(2),strike(1)*fwd(1),0,expiry(2)-expiry(1),option_price(1)/discount_factor);
                      blsimpv(fwd(2),strike(2)*fwd(1),0,expiry(2)-expiry(1),option_price(2)/discount_factor)]


%%
%%%%%%%%%% 2.4 SPOT & FORWARD SMILE %%%%%%%%%%%%%%%

% Using data obtained in 2.1, we can plot the two skews:
figure(2)
plot(K(:,1),V(:,1),'r',K(:,1),ModelVol(:,1),'b','LineWidth',1.5);
title('Volatility smiles');
legend('Forward volatility','Spot volatility');
xlabel('Strike K')
grid on

% We observe that the slope of the forward volatility
% is bigger than the slope of the spot volatility,
% as we expected from LVM.

