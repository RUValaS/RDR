function ecart(vadapted,Y,Mx,My)
%DIFF calcule et affiche l'écart entre vadapted et Y
% différence avec Y
Y_norm = Y/abs(max(Y));
Y_diff = abs(vadapted - Y_norm);

uyR = y_to_im(Y_norm, Mx, My);
figure(),imagesc(flip(abs(uyR),1)),colorbar,axis image

uyR = y_to_im(Y_diff, Mx, My);
%figure(),pcolor(abs(uyR)),colorbar,axis equal

%psnr(abs(Y_norm),vadapted)
%ssim(abs(Y_norm),vadapted)
end

