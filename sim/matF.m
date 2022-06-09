function F = matF(J,Q,z,lambda,I)
%MATF Summary of this function goes here
%   Detailed explanation goes here
dz = zeros(J,J,2);
F = zeros(J^2,Q);

for j = 1 : J
    for i = 1:J
        dz(i,j,1:2) = z(j,:)-z(i,:);
    end
end
dz = reshape(dz,J^2,2);

for q=1:Q
    F(: ,q) = exp(-1i*2*pi*(dz(:,1)*I(q,1)+dz(:,2)*I(q,2)));
end % .* aussi non ?
end

