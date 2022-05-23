function [nY,nR,nQ] = dataGen(H,J,N,Pix,true_image,A,Nreal,SNR_visibilites,RATIO)
%DATAGEN genere les data pour kalman


%init matrices
nY = zeros(J^2,N); % observations => R
nX = zeros(Pix,N); % images EXACTES
nR = zeros(J^2,J^2,N);
nQ = zeros(Pix,Pix,N);
% calculs nX et nY

for k=1:N
    if k == 1
        nX(:,1) = true_image;
    else
        nX(:,k) = A*nX(:,k-1);
    end
    b = sqrt(mean(nX(:,k).^2,'all')/RATIO)*randn(size(nX(:,k))); % erreur evolution
    xb = abs(nX(:,k)+b);
    b = nX(:,k)-xb;
    nQ(:,:,k) = (b*b')/Pix;

    Y_wo_noise = H*(nX(:,k)+b);
    Y_wo_noise = reshape(Y_wo_noise,J,J);
    P_e = trace(Y_wo_noise);
    b2 = (sqrt(P_e/2)*randn(J,Nreal) + 1i*sqrt(P_e/2)*randn(J,Nreal))/(sqrt(Nreal)*sqrt(SNR_visibilites));
    BCov = (b2*b2')/J;

    nY(:,k) = reshape((Y_wo_noise + BCov),[],1);
    BCov = BCov(:);
    nR(:,:,k) = (BCov*BCov')/J^2;
end

end

