%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  EXERCISE 4: EURUSD  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% spot
Y0 = 1.18950;

% market expiries; coincide with expiries of the LV matrix
T = [ 0.038  0.085 	 0.170 	 0.247 	 0.499 	 0.751 	 1.000 	 2.000 	 3.003 	 4.003 ];

% forwards at market expiries
Fwd = [1.190422	1.191420 1.193375	1.195351	1.202167	1.209294	1.216577	1.245994	1.269555	1.292629];

% Domestic dicount factors:
Dd = [ 0.9991560 	 0.9984030 	 0.9969780 	 0.9955680 	 0.9905900 	 0.9851000 	 0.9799950 	 0.9580750 	 0.9367260 	 0.9143290 ];

% Foreign discount factors:
Df = [0.999930461	 1.0000145 	 1.0002258 	 1.0004651 	 1.0011388 	 1.0014927 	 1.0023030 	 1.0035777 	 0.9997690 	 0.9936008 ];

% market deltas
Delta = [ 0.1 0.25 0.5 0.75 0.9];

% LV matrix

% NB: the market volatility from the database has got some irregularities.
MktVol = [
0.0779	0.0740	0.0666	0.0671	0.0677	0.0691	0.0702	0.0737	0.0759	0.0776
0.0745	0.0646	0.0649	0.0660	0.0668	0.0685	0.0689	0.0729	0.0751	0.0768
0.0740	0.0671	0.0643	0.0646	0.0666	0.0678	0.0688	0.0727	0.0749	0.0767
0.0720	0.0687	0.0651	0.0637	0.0673	0.0687	0.0697	0.0736	0.0756	0.0772
0.0779	0.0690	0.0670	0.0670	0.0688	0.0701	0.0712	0.0748	0.0768	0.0785
];


% NB: for more regular results, try this data:
% MktVol = [
%     0.06188	0.06087	0.06312	0.06862	0.06925	0.07863	0.08138	0.08237	0.087	0.09188;
%     0.0605	0.0585	0.06	0.06475	0.06475	0.072	0.07375	0.07425	0.07825	0.0825;
%     0.06	0.0575	0.0585	0.063	0.063	0.0685	0.07	0.071	0.0745	0.078;
%     0.0615	0.0595	0.061	0.06575	0.06625	0.071	0.07275	0.07475	0.07775	0.0805;
%     0.06362	0.06263	0.06487	0.07037	0.07175	0.07688	0.07963	0.08313	0.086	0.08813];


%%
%%%%%%% 4.1 %%%%%%%%%%%%%%%%%

% market strikes
K = zeros(length(Delta),length(T));

% find K such that BS_Delta(K,Fwt,T,MktVol) = Delta

for i = 1:length(Delta)
    for j = 1:length(T)
        K(i,j) = fzero(@(Strike) blsdelta(Fwd(j),Strike,0,T(j),MktVol(i,j))-(1-Delta(i)), Fwd(j));
    end
end

% normalized market strikes
[rows, cols] = size(K);
K_norm = K ./ repmat(Fwd, rows, 1);

% Dupire solver settings
N = 20;
M = 300;
K_min = 0.5;
K_max = 2.5;
Scheme = 'cn';

% calibration settings
Threshold = 0.0010;
MaxIter = 100;

[V, ModelVol, MaxErr] = calibrator_mod(T,K_norm,MktVol,Threshold,MaxIter,N,M,K_min,K_max,Scheme);

% plot calibration error at each iteration
figure (1)
plot(MaxErr);
grid on
title('Trend of the approximation error')

% plot local volatility function vs market implied volatility
figure (2)
plot(K(:,1),MktVol(:,1),'r','LineWidth',1.5);
hold on
%% 
%%%%%%%%%%% 4.2 %%%%%%%%%%%%%%%%%

Y1 = 1.01*Y0;

Fwd_new = Y1*Df./Dd;

% Having updated Y0 and forward, the procedure is the same as before:
K_new = zeros(length(Delta),length(T));

% find K such that BS_Delta(K,Fwt,T,MktVol) = Delta

for i = 1:length(Delta)
    for j = 1:length(T)
        K_new(i,j) = fzero(@(Strike) blsdelta(Fwd_new(j),Strike,0,T(j),MktVol(i,j))-(1-Delta(i)), Fwd_new(j));
    end
end

% normalized market strikes
K_new_norm = K ./ repmat(Fwd, rows, 1);

% We see that the difference K_norm - K_new_norm
% is zero, as hypotesis.

[V_new, ModelVol_new, MaxErr_new] = calibrator_mod(T,K_new_norm,MktVol,Threshold,MaxIter,N,M,K_min,K_max,Scheme);

plot(K_new(:,1),ModelVol_new(:,1),'b','LineWidth',1.5);
grid on
legend ('Volatility with spot Y(0)','Volatility with spot Y''(0)')
xlabel('Strike')
ylabel('Volatility')

% We see the behaviour of a sticky-delta market.