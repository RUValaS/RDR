function A = matA(z, I_q, Q, f, c)

%%%%%%%%%%%%%%%%%%%%%
% Q sources
% J antennes
% z positions des antennes Jx2 (n0)x(coords)
% I_q vecteur cosinus de dir
% s_q(n,k) <=> signal of qth source @ time sample n, freq f_k
%%%%%%%%%%%%%%%%%%%%%

J = size(z, 1);
A = zeros(J,Q);

for i=1:J
    for q=1:Q
        A(i,q) = exp(-1j*2*pi*(f/c)*(z(i,:))*transpose(I_q(q, :)));
    end
end