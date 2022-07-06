addpath(genpath(pwd));
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRES
Mx = 10; %My auto (ratio)
J_b = 6; % On suppose le réseau d'antennes carré de côté
% max 11 GPU atm
c= 3e8;
f = 1e9;
lambda = c/f;
dist = 10*lambda/2.1; % J_b antennes espacées de dist 
N = 6; % itérations Kalman
W = 6; % nombre de réalisations
SNR = 10;
iMEM = 5; % nombre itérations MEM 
RATIO = 10e14; % ratio erreur Kalman
mode = 'CPU';
Nreal = 1000;
% A_ev matrice d'évolution définie plus loin (need D ligne 35)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% choix / gen image
img = imread('bh_square.jpg');
adapted = adaptImg(img, Mx); 
vadapted = vectorize(adapted);
clear img

%%% Calcul paramètres
My = size(adapted, 2);
z = antennes(J_b, dist);
J = size(z, 1);
D = numel(adapted);
I = matI(Mx, My);

A_ev = eye(D);
A = matA(z,I,D,f,c);

erreurs = [0.001 0.003 0.005 0.01 0.03 0.05 0.1 0.3];
npos = numel(erreurs);

% calcul invariants kalman
H = matF(J,D,z,lambda,I);
[nY,nR,nQ,tX] = dataGen_im(H,J,N,D,vadapted,A_ev,Nreal,SNR,RATIO);

% X_0 = true_image ;
X_0 = MEM(MVDR(A,reshape(nY(:,1),J,J),Mx,My),iMEM);
X_0 = normarr(X_0);
tX(:,1) = X_0;

P_0 = X_0*X_0' - mean(X_0,'all');

% def listes vides
wK = zeros(W,N,J^2,D);
wK_err = zeros(W,N,J^2,D);
wX = zeros(W,N,D);
wX_err = zeros(W,N,D);
e_x = zeros(W,1);
e_K = zeros(W,1);
e_x_e = zeros(W,1);
de_x = zeros(W,1);
de_K = zeros(W,1);
de_x_e = zeros(W,1);

for err = 1:npos
    for realisation = 1:W
        % 1. Gen données aléatoires
        poids = erreurs(err);
        
        z_err = z+ randn(size(z))*poids;
        I_err = I;
        % I_err(4,:) = I_err(4,:) + randn(size(I_err(4,:)))*poids;
%         I_err = I + randn(size(I))*poids;
        H_err = matF(J,Pix,z_err,lambda,I_err);

        C = H ./ H_err;
        
        % 2. Tourner Kalman_XPU
        [X,K] = Kalman_CPU_V3(A_ev,H,X_0,P_0,nY,nR,nQ,D,N);
        [X_e,K_e] = Kalman_CPU_V3(A_ev,H_err,X_0,P_0,nY,nR,nQ,D,N);

        % 3. Calcul erreurs
        errX = zeros(N,1);
        errK = zeros(N,1);
        errX_e = zeros(N,1);
        for k = 1:N
            errX(k) = fro(abs(X(:,k)/max(abs(X(:,k)))) - tX(:,k));
            errX_e(k) = fro(abs(X_e(:,k)/max(abs(X_e(:,k)))) - tX(:,k));
            errK(k) = fro(K_err(:,:,k)-K(:,:,k));
        end
    end
end