addpath('./fonctions_e1','Kalman')
addpath('./fonctions_rec')
addpath('./img_t')
addpath('./fonctions_clean')

close all
% clear all

img = imread('points.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRES
Mx = 32; %My auto (ratio)
J_b = 6; % On suppose le réseau d'antennes carré de côté
% max 11 GPU atm
c= 3e8;
f = 100e6;
lambda = c/f;
dist = lambda/2.1; % J_b antennes espacées de dist (normalisé)
N = 4;
SNR = 0.1;
iMEM = 1000; % nombre itérations MEM 
RATIO = 10e14; % ratio erreur Kalman
mode = 'GPU';
% A_ev matrice d'évolution définie plus loin (need Pix ligne 69)

% OPTIONS DE LANCEMENT
inter_fig = 0; % figures autres que originale / traitée
mode_R = 0; % 0 FFT - 1 S - 2 direct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%A FAIRE
% vectoriser les données à traiter sous forme "truc x 1"
% permet de faciliter les calculs, notamment
% pour le
%            y = H*x
% size(H) = J^2*Q
% H matrice FFT -> réponse du réseau à la réception des données
% y = R
% R(xx) = H.X
%
% Xdirty = pinv(H)*R -> a tester
% pinv -> pseudo inversion, plus efficace avec des matrices
% mal conditionnée, + fonctionne avec des matrices pas ²
%
% A : J x Q
% S : Q x N
% X : J x N
% b : J x N (pour avoir une consistance dans les calculs)
% R = (A.S + b)^H
% <=>
% R = (X.X^H)/N + (b.b^H)/N
% R =   R(xx)   +  R(noise)
% On peut tester les deux pour voir que c'est équivalent
%  pas  exactement égaux tho -> processus aleatoires
%
% Puissance du signal
% Ps -> puissance par antenne
% dans R -> diagonale -> autocorrélations
% -> puissances
% Ps = trace(R)/J

% adaptation image
adapted = adaptImg(img, Mx);
vadapted = vectorize(adapted);

%calcul des paramètres
My = size(adapted, 2);
z = antennes(J_b, dist);
J = size(z, 1);
Pix = numel(adapted);
I = matI(Mx, My);

A_ev = eye(Pix);

%figure();imagesc(adapted,'CDataMapping','scaled');colorbar;axis image

%calcul des matrices / def matrices vides
switch mode_R
    case 0
        [F,R,noise_cov] = matR_FT(J,z,Pix,vadapted,lambda,I,SNR,N);
    otherwise
        [F,R,noise_cov] = matR_FT(J,z,Pix,vadapted,lambda,I,SNR,N); %default to 0
end

A = matA(z,I,Pix,f,c);
Y = zeros(Mx*My, 1);
%B = ifft(F);

% reconstruction 
Y_bf = beamforming(A,R,Mx,My);
Y_mvdr = MVDR(A,R,Mx,My);
Y_aar = AAR(A,R,Mx,My);

% ecart(vadapted,Y_bf,Mx,My);
% ecart(vadapted,Y_mvdr,Mx,My);
% ecart(vadapted,Y_aar,Mx,My);

tic
Y_em = MEM(Y_bf,iMEM);
toc

Y_em=normarr(Y_em);
Y_bf = normarr(Y_bf);
Y_mvdr = normarr(Y_mvdr);
Y_aar = normarr(Y_aar);

y = y_to_im(Y_em,Mx,My);
y_d = y_to_im(Y_bf, Mx, My);
figure();imagesc(flip(abs(y_d),1));colorbar; axis image;title("BF")
figure();imagesc(flip(abs(y),1));colorbar; axis image;title("MEM")
psnr(abs(Y_bf),vadapted)
psnr(abs(Y_em),vadapted)
ssim(abs(Y_bf),vadapted)
ssim(abs(Y_em),vadapted)

% Kalman
tic
%X = KalmainV3(vadapted,J,z,I,N,f,c,lambda,Mx,My,iMEM,SNR,A_ev,RATIO,mode);
toc
%dispKalman(X,N,Mx,My);
%Pix
if inter_fig
    figure(),image(img);
    figure(),pcolor(imag(A)),colorbar
    figure(),pcolor(abs(R)),colorbar
end