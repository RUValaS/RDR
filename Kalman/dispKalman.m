function  dispKalman(X,N,Mx,My)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figure();
for k=1:N
    subplot(1,N,k);imshow(y_to_im(abs(X(:,k)),Mx,My));colorbar;axis image;colormap("parula")
end
end

