function B = matB(z,Q,J,lambda,I)
%MATB Summary of this function goes here
%   Detailed explanation goes here

for dir=1:Q
   for n=1:J
      for k=1:J
         B(dir) = exp(-1i*(2*pi/lambda)*(z(n,:)-z(k,:))*(I(dir,:)')); 
      end
   end
end

end

