function [X,out_2] = Kalman_CPU_V3(A,H,X_0,P_0,nY,nR,nQ,Pix,n,C)
%KALMAN_CPU Estime via un filtre de Kalman la sortie X-out depuis X_in
%   variables:
%   - H : observation (J^2 x Pix) -> matrice de transformation FFT non
%           uniforme
%   - A : matrice d'évolution (Pix x Pix)
%   - X_in : image initiale (Pix x 1)
%   - R : matrice de covariance du bruit d'observation
%   - Q : matrice de covariance de l'erreur sur le modèle d'évolution
%   - I : nombre itérations

X = zeros(Pix,n);
P = zeros(Pix,Pix,n);

Ks = zeros(Pix,size(H,1),n);
% e_x = zeros(n,1);

%%% INITIALISATION
X(:,1) = X_0;
P(:,:,1) = P_0;
% figure();
for k=2:n
%     fprintf('itération : %u \n',k);
    %%% PREDICTION
    Xp = A*X(:,k-1) + transpose(mean(nQ(:,:,k-1)));
    Pp = A*P(:,:,k-1)*(A.') + nQ(:,:,k);
    yp = H*Xp + mean(nR(:,k));
    %%% KALMAN GAIN
    K = Pp*(H.')*pinv(H*Pp*(H.') + nR(:,k));
    %%% ESTIMATION
    X(:,k) = Xp + K*(nY(:,k)-yp);
    P(:,:,k) = Pp - K*H*Pp;
%     subplot(n-1,1,k-1);
%     imshow(abs(K));colorbar;axis image;colormap("parula");title(mean(K, 'all'))
    Ks(:,:,k) = K; % loggin K
%     e_x(k) = fro(K*( (ones(size(C)) - C) .* H )*Xp);
end
out_2 = Ks;
% out_2 = e_x;
end