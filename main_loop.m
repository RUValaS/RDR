addpath(genpath(pwd));
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRES
Mx = 2; %My auto (ratio) -- nombre de pixels
J_b = 6; % On suppose le réseau d'antennes carré de côté
% max 11 GPU atm
c= 3e8;
f = 1e9;
lambda = c/f;
dist = 10*lambda/2.1; % J_b antennes espacées de dist 
N = 7; % itérations Kalman
W = 1; % nombre de réalisations
SNR = 10;
iMEM = 2; % nombre itérations MEM 
RATIO = 10e14; % ratio erreur Kalman
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
dz = dz_c(J,z);

A_ev = eye(D);
A = matA(z,I,D,f,c);

% erreurs = [0.001 0.003 0.005 0.01 0.03 0.05 0.1 0.3];
% erreurs = [0.001 0.003 0.005 0.008 0.01 0.03];
% erreurs = [0.001 0.3 0.05];
erreurs = [0.001];
npos = numel(erreurs);

% calcul invariants kalman
H = matF(J,D,z,lambda,I);

[nY,nR,nQ,tX] = dataGen_im(H,J,N,D,vadapted,A_ev,Nreal,SNR,RATIO);
% X_0 = true_image ;
% X_mvdr = MEM(MVDR(A,reshape(nY(:,1),J,J),Mx,My),iMEM);
X_mvdr = MVDR(A,reshape(nY(:,1),J,J),Mx,My);
X_mvdr = normarr(abs(X_mvdr));
% X_0 = vadapted;
X_0 = X_mvdr;
% tX(:,1) = vadapted;

P_0 = X_0*X_0' - mean(X_0,'all');

% Nouveau calcul X_0 et P_0 exactes
% NX_0 = repmat(vadapted,1,Nreal);
% vX_0 = abs(0.001*randn(D,Nreal) + NX_0);
% errX_0 = vX_0 - NX_0;
% 
% X_0 = vX_0(:,randi(Nreal));
% P_0 = (errX_0 * errX_0')/Nreal;

clear NX_0 vX_0 errX_0

% def listes vides
wK = zeros(D,J^2,N,W);
whos wK
wK_err = zeros(D,J^2,N,W);
wX = zeros(D,N,W);
wX_err = zeros(D,N,W);
e_x = zeros(W,1);
e_K = zeros(W,1);
e_x_e = zeros(W,1);
e_th = zeros(W,1);
e_xx = zeros(W,1);

errX = zeros(N,1);
errK = zeros(N,1);
errX_e = zeros(N,1);
errth = zeros(N,1);
errXX = zeros(N,1);

optot = npos*W;
tic
for err = 1:npos
    poids = erreurs(err);
    for realisation = 1:W
        elapsedTime = toc;
        opDone = W*(err-1) + realisation;
        opRemain = optot - opDone + 1;
        tMoyOp = elapsedTime/opDone;
        ETA = tMoyOp*opRemain;
        fprintf('Erreur : %u/%u -- Réalisation : %u/%u -- Time : %.4f -- ETA : %.4f\n',err,npos,realisation,W,elapsedTime,ETA);
        % 1. Gen données aléatoires
        z_err = z;
%         z_err = z+ randn(size(z))*poids;
%         I_err = I;
        I_err(4,:) = I_err(4,:) + randn(size(I_err(4,:)))*poids;
%         I_err = I + randn(size(I))*poids;
        H_err = matF(J,D,z_err,lambda,I_err);

        % Calcul Hcorr
        % 1. poser pb
        
%         y = H_err*X_mvdr;
        y = nY(:,1);
        
        Iv = reshape(I_err,[],1);

        prb = @(x)prb_f(x,y,Iv,dz,J,D,lambda,X_mvdr);


        % 2. init
%         dbiasPrio = zeros(1,2*D)*err/100;
%         xparfait = I_err-I;
%         dbiasPrioInit = reshape(xparfait,1,2*D) + dbiasPrio;
%         x0 = [dbiasPrioInit reshape(X_mvdr,1,[])];
%         lb = [-0.1*ones(1,2*D) zeros(1,D)];
%         ub = [0.1*ones(1,2*D) ones(1,D)];

        x0 = zeros(1,2*D);
        lb = [-0.1*ones(1,2*D)];
        ub = [0.1*ones(1,2*D)];


        % 3. lsqnonlin
        options = optimoptions(@lsqnonlin,'Display','Iter','FiniteDifferenceType','central');
%         options = optimoptions(@lsqnonlin,'MaxIterations',1000,'StepTolerance',1e-18,'Display','Iter','FunctionTolerance',1e-20,'MaxFunctionEvaluations',1e6,'FiniteDifferenceType','central');
        % ,'Algorithm','levenberg-marquardt'
        x = lsqnonlin(prb,x0,lb,ub,options);

        % 4. nouveau X_0, calcul Hcorr
%         X_0 = x(2*D+1:3*D);
        Icorr = - reshape(x(1:2*D),D,2) + I_err;
        Hcorr = matF(J,D,z_err,lambda,Icorr);
        
        amel(I_err,Icorr,I)
%         C = H ./ H_err;
        
        % 2. Tourner Kalman_XPU
        [X,K,Xp] = Kalman_CPU_V6(A_ev,H,X_0,P_0,nY,nR,nQ,D,N);
        [X_e,K_e,Xp_e] = Kalman_CPU_V6(A_ev,H_err,X_0,P_0,nY,nR,nQ,D,N);

        % 3. Calcul erreurs
        for k = 1:N
            errX(k) = fro(abs(X(:,k)/max(abs(X(:,k)))) - tX(:,k));
            errX_e(k) = fro(abs(X_e(:,k)/max(abs(X_e(:,k)))) - tX(:,k));
            errXX(k) = fro(abs(X(:,k)/max(abs(X(:,k)))) - abs(X_e(:,k)/max(abs(X_e(:,k)))));
            errK(k) = fro(K_e(:,:,k)-K(:,:,k));
            if k>=2
                errth(k) = fro((K_e(:,:,k) - K(:,:,k))*nY(:,k) + (eye(D) - K_e(:,:,k)*H_err)*Xp_e(:,k) - (eye(D) - K(:,:,k)*H)*Xp(:,k));
            end
        end % end calculs erreurs
       
        % 4. append in vect errs
        wK(:,:,:,realisation) = K;
        wK_err(:,:,:,realisation) = K_e;
        wX(:,:,realisation) = X;
        wX_err(:,:,realisation) = X_e;

        e_x(realisation) = mean(errX(5:end));
        e_K(realisation) = mean(errK(5:end));
        e_x_e(realisation) = mean(errX_e(5:end));
        e_th(realisation) = mean(errth(5:end));
        e_xx(realisation) = mean(errXX(5:end));
    end % end real

    de_x = std(errX);
    de_K = std(errK);
    de_x_e = std(e_x_e);
    de_xx = std(e_xx);

    % Sauvegarde données
    nom = strcat(strrep(sprintf('%.3f',poids),'.','_'),'.mat');
%     save(nom,"wK","wK_err","wX","wX_err","e_x","e_K","e_x_e","de_x","de_K","de_x_e","e_th","de_xx","e_xx"); % ajouter Xp Xp_e
end % end err

function s = prb_f(x,y,Iv,dz,J,D,lambda,xMVDR)
%         prb = @(x)([real(y); imag(y)] - [cos(-2*(pi/lambda)*(dz(:,1)*(Iv(1:D) - x(1:D)) + dz(:,2)*(Iv(D+1:D*2) - x(D+1:2*D))));
%             sin(-2*(pi/lambda)*(dz(:,1)*(Iv(1:D) - x(1:D)) + dz(:,2)*(Iv(D+1:D*2) - x(D+1:2*D))))] * x(2*D+1:3*D));
Hmodif = zeros(J^2,D);
for q=1:D
    Hmodif(: ,q) = exp(-1i*2*pi/lambda*(dz(:,1)*(Iv(q) - x(q))+dz(:,2)*(Iv(D+q) - x(D+q))));
end % .* aussi non ?
ymodif = Hmodif*xMVDR;
s = [real(y);imag(y)] - [real(ymodif);imag(ymodif)];
    
end