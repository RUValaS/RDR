function [X_err, tX] = KalmainV4(true_image,J,z,I,N,f,c,lambda,Mx,My,iMEM,SNR_visibilites,A,RATIO,mode)
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
poids = 0.001;
X_0 = true_image ;
tX(:,1) = X_0;

P_0 = X_0*X_0' - mean(X_0,'all');
% P_0 = P_0 + wishrnd(eye(Pix),Pix+2)*poids;
% P_0 = eye(Pix);

%ajout erreur sur A (!!!APRES!!! dataGen => sinon c'est pas une erreur)
% H_err = abs(H .*( poids*exp(1j*randn(size(H))) ) );
% z_err = z+ randn(size(z))*poids;
I_err = I;
% I_err(4,:) = I_err(4,:) + randn(size(I_err(4,:)))*poids;
I_err = I + randn(size(I))*poids;
H_err = matF(J,Pix,z,lambda,I_err);

C = H ./ H_err;


%%%%%%%%%%%% Kalman
X_err=-1;
if mode=='CPU'
    [X_err,K_err] = Kalman_CPU_V3(A,H_err,X_0,P_0,nY,nR,nQ,Pix,N,C);
    [~,K] = Kalman_CPU_V3(A,H,X_0,P_0,nY,nR,nQ,Pix,N);
end
if mode=='GPU'
    X_err = Kalman_GPU_V5(A,H_err,X_0,P_0,nY,nR,nQ,Pix,N);
end

errX = zeros(N,1);
errK = zeros(N,1);
errP = zeros(N,1);
for k=1:N
%         Psr(k) = psnr(abs(X(:,k)),tX(:,k));
    errX(k) = fro(abs(X_err(:,k)/max(abs(X_err(:,k))))-tX(:,k));
    errK(k) = fro(K_err(:,:,k)-K(:,:,k));
    errP(k) = fro(abs(X_err(4,k)/max(abs(X_err(:,k))))-tX(4,k));
end
% figure();plot(errX);title('erreur X')
% fprintf('||e_X|| = %e ---- ||errK|| = %e ---- ||errP|| = %e \n',mean(errX(4:end)),mean(errK(4:end)),mean(errP(4:end)));
fprintf('||e_X|| = %e ---- ||errK|| = %e \n',mean(errX(4:end)),mean(errK(4:end)));
% fprintf('||e_X|| = %e ---- ||e_X||th = %e \n',mean(Psr(4:end)),mean(e_x(4:end)));

%{
figure();plot(e_x);title("Ev e_x")
hold on
plot(Psr);
hold off
%}

%{
z_err = z+ randn(size(z))*0.01;
H = matF(J,Pix,z,lambda,I);
H_err = matF(J,Pix,z_err,lambda,I);
%}

end

