clc;
clearvars;

% market expiries
 T = [0.140  0.216 	 0.312 	 0.466 	 0.734 	 0.984 	 1.060 	 1.482 	 1.981 	 2.058 	 2.978 ];

% forwards at market expiries
Fwd = [2694.93	2696.09	2699.20	2698.08	2704.40	2709.25	2712.81	2720.39	2730.30	2732.90	2763.60];

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
% normalized market strikes
[rows, cols] = size(K);
K_norm = K ./ repmat(Fwd, rows, 1);

% repmat crea una matrice del tipo
%[A B C D;
%A B C D;
% A B C D]

% Dupire solver settings
N = 10;
M = 200;
K_min = 0.1;
K_max = 3;
Scheme = 'cn'; %crank-nicholson

% calibration settings
Threshold = 0.0010;
MaxIter = 100;

tstart = tic;
[V, ModelVol, MaxErr] = calibrator(T,K_norm,MktVol,Threshold,MaxIter,N,M,K_min,K_max,Scheme);
telapsed(1)=toc(tstart);
tstart=tic;
[V_mod, ModelVol_mod, MaxErr_mod] = calibrator_mod(T,K_norm,MktVol,Threshold,MaxIter,N,M,K_min,K_max,Scheme);
telapsed(2)=toc(tstart);

% plot calibration error at each iteration
figure(1)
plot(MaxErr);
figure(2)
plot(MaxErr_mod);

% plot local volatility function vs market implied volatility
figure(3)
plot(K(:,1),V(:,1),K(:,1),MktVol(:,1)); %V in blue, Mkt in red
figure(4)
plot(K(:,1),V_mod(:,1),K(:,1),MktVol(:,1)); %V in blue, Mkt in red

telapsed