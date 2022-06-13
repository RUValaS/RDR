function [X, tX] = KalmainV4(true_image,J,z,I,N,f,c,lambda,Mx,My,iMEM,SNR_visibilites,A,RATIO,mode)
%KALMAIN Estime via un filtre de Kalman la sortie X-out depuis X_in
%normalise en plus les images
%utilise true_image en X_0
%renvoie les images

%%%%%%%%%%%%% def / calcul paramètres
Pix = Mx*My;
H = matF(J,Pix,z,lambda,I);
Nreal = 1000;

% génération n images
[nY,nR,nQ,tX] = dataGen_im(H,J,N,Pix,true_image,A,Nreal,SNR_visibilites,RATIO);

% génération image initiale
% ajout bruit image
poids = 0.1;
X_0 = true_image ;
tX(:,1) = X_0;

P_0 = X_0*X_0';
% P_0 = P_0 + wishrnd(eye(Pix),Pix+2)*poids;
% P_0 = eye(Pix);

%ajout erreur sur A (!!!APRES!!! dataGen => sinon c'est pas une erreur)
% H_err = abs(H .*( poids*exp(1j*randn(size(H))) ) );
z_err = z+ randn(size(z))*poids;
H_err = matF(J,Pix,z_err,lambda,I);


%%%%%%%%%%%% Kalman
X=-1;
if mode=='CPU'
    X = Kalman_CPU_V2(A,H_err,X_0,P_0,nY,nR,nQ,Pix,N);
end
if mode=='GPU'
    X = Kalman_GPU_V5(A,H_err,X_0,P_0,nY,nR,nQ,Pix,N);
end

Psr = zeros(N,1);
for k=1:N
        Psr(k) = psnr(abs(X(:,k)),tX(:,k));
%     Psr(k) = ssim(abs(X(:,k)/max(abs(X(:,k)))),tX(:,k));
end
figure();plot(Psr);title('PSNR = f(it)')
mean(Psr(6:end))

%{
z_err = z+ randn(size(z))*0.01;
H = matF(J,Pix,z,lambda,I);
H_err = matF(J,Pix,z_err,lambda,I);
%}

end

