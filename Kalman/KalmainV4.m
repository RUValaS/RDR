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
P_0 = nQ(:,:,1);
X_0 = true_image;

%%%%%%%%%%%% Kalman
X=-1;
if mode=='CPU'
    X = Kalman_CPU_V2(A,H,X_0,P_0,nY,nR,nQ,Pix,N);
end
if mode=='GPU'
    X = Kalman_GPU_V4(A,H,X_0,P_0,nY,nR,nQ,Pix,N);
end

end

