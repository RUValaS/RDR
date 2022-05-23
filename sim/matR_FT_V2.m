function [R,noise_cov] = matR_FT_V2(J,z,Q,vadapted,lambda,I,SNR,N)
%MATR_FT MATR_FT sans retour de F

dz = zeros(J,J,2);
F = zeros(J^2,Q);

for j = 1 : J
    for i = 1:J
        dz(i,j,1:2) = z(j,:)-z(i,:);
    end
end
dz = reshape(dz,J^2,2);

for q=1:Q
      F(: ,q) = exp(-1i*2*pi/lambda*(dz(:,1)*I(q,1)+dz(:,2)*I(q,2)));
end

R_wo_noise = F*vadapted;
R_wo_noise = (reshape(R_wo_noise,J,J));
Ps = trace(R_wo_noise);  % puissance moyenne sur chaque antenne
Pb = Ps/SNR;
noise = sqrt(Pb/2)*(randn([J,N]) +1i*randn([J,N]))/sqrt(N);
noise_cov = noise*noise'/J;
R = R_wo_noise + noise_cov;

end

