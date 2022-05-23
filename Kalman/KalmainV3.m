function X = KalmainV3(true_image,J,z,I,N,f,c,lambda,Mx,My,iMEM,SNR_visibilites,A,RATIO,mode)
%KALMAIN Estime via un filtre de Kalman la sortie X-out depuis X_in
%normalise en plus les images

%%%%%%%%%%%%% def / calcul paramètres
Pix = Mx*My;
Arr = matA(z,I,Pix,f,c);
H = matF(J,Pix,z,lambda,I);
Nreal = 1000;

% génération n images
[nY,nR,nQ] = dataGen_wo_noise(H,J,N,Pix,true_image,A,Nreal,SNR_visibilites,RATIO);

% génération image initiale
[R,~] = matR_FT_V2(J,z,Pix,true_image,lambda,I,SNR_visibilites,Nreal);
Y_mvdr = MVDR(Arr,R,Mx,My);
X_0 = MEM(Y_mvdr,iMEM);
X_0 = normarr(X_0);
P_0 = nQ(:,:,1);

%%%%%%%%%%%% Kalman
X=-1;
if mode=='CPU'
    X = Kalman_CPU_V2(A,H,X_0,P_0,nY,nR,nQ,Pix,N);
end
if mode=='GPU'
    X = Kalman_GPU_V4(A,H,X_0,P_0,nY,nR,nQ,Pix,N);
end

end

