function y = MVDR(A,R,Mx,My)
%MVDR reconstruction de l'image par méthode de maximum de vraisemblance
Rm = pinv(R);
%w = Rm*A/((A')*Rm*A);
%y = (w')*R*w;

y = zeros(Mx*My, 1);

for dir=1:(Mx*My)
    Ad = A(:, dir);
    w = Rm*Ad/((Ad')*Rm*Ad);
    y(dir) = (w')*R*w;
end
end

