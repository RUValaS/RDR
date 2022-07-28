function F = matF(J,Q,z,lambda,I)
%MATF Summary of this function goes here
%   Detailed explanation goes here

F = zeros(J^2,Q);

dz = dz_c(J, z);

for q=1:Q
    F(: ,q) = exp(-1i*2*pi/lambda*(dz(:,1)*I(q,1)+dz(:,2)*I(q,2)));
end % .* aussi non ?
end

