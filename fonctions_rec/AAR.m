function y = AAR(A,R,Mx,My)
%AAR Reconstruction par méthode d'adaptative angle reconstruction
Rm = pinv(R);
Rm2 = Rm^2;

y = zeros(Mx*My, 1);

for dir=1:(Mx*My)
    Ad = A(:, dir);
    
    y(dir) = (Ad')*Rm*Ad/((Ad')*Rm2*Ad);
end
end

