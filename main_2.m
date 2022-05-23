addpath(genpath(pwd))

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRES
Mx = 3; %My auto (ratio)
J_b = 5; % On suppose le réseau d'antennes carré de côté
% max 11 GPU atm
c= 3e8;
f = 100e6;
lambda = c/f;
dist = lambda/2.1; % J_b antennes espacées de dist (normalisé)
N = 3;
SNR = 0.1;
iMEM = 1; % nombre itérations MEM 
RATIO = 10e14; % ratio erreur Kalman
mode = 'CPU';
% A_ev matrice d'évolution définie plus loin (need Pix ligne 69)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% choix / gen image
% image = imread('bh_pic.jpg');
% adapted = adaptImg(img, Mx);

vadapted = reshape(linspace(0,1,9),9,1);
adapted = reshape(vadapted, 3, 3);

%%% Calcul paramètres
My = size(adapted, 2);
z = antennes(J_b, dist);
J = size(z, 1);
Pix = numel(adapted);
I = matI(Mx, My);

A_ev = eye(Pix);

%%% Kalman direct en fait mdr
X = KalmainV3(vadapted,J,z,I,N,f,c,lambda,Mx,My,iMEM,SNR,A_ev,RATIO,mode);
dispKalman(X,N,Mx,My);