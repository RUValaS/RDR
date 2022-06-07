function  dispKalman_t(X,tX,N,Mx,My)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figure();
for k=1:N
    subplot(2,N,k);imshow(y_to_im(abs(tX(:,k)),Mx,My));colorbar;axis image;colormap("parula")
    subplot(2,N,k+N);imshow(y_to_im(abs(X(:,k)/max(abs(X(:,k)))),Mx,My));colorbar;axis image;colormap("parula");title(ssim(abs(X(:,k)/max(abs(X(:,k)))),tX(:,k)))
end
end

