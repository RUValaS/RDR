function z = antennes(J_b, dist)
z = zeros(J_b, 2); % num, pos
for i=1:J_b
    for j=1:J_b
        z((i-1)*J_b+j, :) = [(i-1)*dist; (j-1)*dist];
    end
end
