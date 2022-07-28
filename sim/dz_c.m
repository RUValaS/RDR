function dz_o = dz_c(J,z)
dz_o = zeros(J,J,2);
for j = 1 : J
    for i = 1:J
        dz_o(i,j,1:2) = z(j,:)-z(i,:);
    end
end
dz_o = reshape(dz_o,J^2,2);
end