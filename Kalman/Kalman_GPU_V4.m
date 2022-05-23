function X = Kalman_GPU_V4(A,H,X_0,P_0,nY,nR,nQ,Pix,n)
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

%%% INITIALISATION
X(:,1) = X_0;
P(:,:,1) = P_0;

for k=2:n
    fprintf('itération : %u \n',k);
    %%% PREDICTION
    Xp = A*X(:,k-1);
    Pp = A*P(:,:,k-1)*(A.') + nQ(:,:,k);
    %%% KALMAN GAIN
    S = gpuArray(H*Pp*(H.') + nR(:,k));
    T = gpuArray(Pp*(H.'));
    K = gather(T*pinv(S));   
    clear S
    clear T
    %%% ESTIMATION
    X(:,k) = Xp + K*(nY(:,k)-H*Xp);
    P(:,:,k) = Pp - K*H*Pp;
end
end

