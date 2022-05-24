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
N = 25;
SNR = 10;
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

% A_ev = eye(Pix);
A_ev = zeros(Pix);
A_ev(7,1) = 1; A_ev(4,2)=1; A_ev(1,3)=1; A_ev(8,4)=1;A_ev(5,5)=1;A_ev(2,6)=1;A_ev(9,7)=1;A_ev(6,8)=1;A_ev(3,9)=1;

%%% Kalman direct en fait mdr
[X,tX] = KalmainV4(vadapted,J,z,I,N,f,c,lambda,Mx,My,iMEM,SNR,A_ev,RATIO,mode);
% dispKalman_t(X,tX,N,Mx,My);

z_err = z+ randn(size(z))*0.01;
H = matF(J,Pix,z,lambda,I);
H_err = matF(J,Pix,z_err,lambda,I);