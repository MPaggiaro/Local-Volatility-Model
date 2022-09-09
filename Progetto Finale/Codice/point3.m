%%%%%%%%%%%%%%%%%%% SCRIPT PER IL TITOLO FAIL %%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% market expiries
T =[0.0082 0.0192];

% forwards at market expiries
Fwd = [4.50708	4.50708];

% discount factors at market expiry
D =[ 1.0000000 	 1.0000000 ];

% market strikes
K = [
    4.41060	4.41060
    4.45479	4.45479
    4.50809	4.50809
    4.56930	4.56930
    4.63267	4.63267
    ];

% LV matrix
MktVol = [
    0.1903	0.1203
    0.1975	0.1275
    0.2050	0.1350
    0.2175	0.1475
    0.4323	0.1623
    ];

% spot price
S0=4.50708;

% normalized market strikes
[rows, cols] = size(K);
Knorm = K ./ repmat(Fwd, rows, 1);

%%%%%%% 2.1 CALIBRATION DEL MODELLO %%%%%%
% Dupire solver setting
N = 1000;
M = 200;
Kmin = 0.7;
Kmax = 1.3;

% calibration settings
tol = 0.001;
MaxIter = 100;

% calibration
%%%%%% ERRORE !!!!%%%%% [v, ModelVol, MaxErr] = calibrator(T,Knorm,MktVol,tol,MaxIter, N,M,Kmin, Kmax, 'cn');
C=0.*K;
C(:,1)=blsprice(S0,K(:,1),0,T(1),MktVol(:,1));
C(:,2)=blsprice(S0,K(:,2),0,T(2),MktVol(:,2));

%Grafico andamento della volatilità in funzione dello strike
nodi=linspace(Kmin,Kmax,600);
figure
[smile]=interp_flat_extrap(Knorm(:,1),MktVol(:,1),nodi,'spline');
plot(nodi,smile, 'b', 'LineWidth', 4);
title('Market Volatilities - normalized strikes')
grid on
hold on
[smile]=interp_flat_extrap(Knorm(:,2),MktVol(:,2),nodi,'spline');
plot(nodi,smile, 'r', 'LineWidth', 4);
legend('Maturity 18 dicembre','Maturity 22 dicembre');

%Grafico K-Call
figure
plot(K(:,1),C(:,1), 'b', 'LineWidth', 4);
title('Call prices - strikes')
hold on
grid on
plot(K(:,2),C(:,2), 'r', 'LineWidth', 4);
legend('Maturity 18 dicembre','Maturity 22 dicembre');
