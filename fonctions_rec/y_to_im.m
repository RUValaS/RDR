function uyR = y_to_im(Y,Mx,My)
%Y_TO_IM met le vecteur Y sous forme d'image reconstruite (les méthodes
%utilisées induisant une rotation de l'image, celle-ci est corrigée)

uyR = zeros(Mx, My);

uy = reshape((Y),My,Mx);
tuy = uy';

for c=1:My
    uyR(:, My-c+1) = tuy(:, c);
end
end

