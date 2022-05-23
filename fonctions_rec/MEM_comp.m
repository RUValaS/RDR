function x = MEM_comp(x_d,q,true_im)

    function y = f(x)
       y = -x.*log(x); 
    end

vSSIM = zeros(q,1);
vPSNR = zeros(q,1);

s = zeros(size(x_d));
x = x_d;
%close all
%figure(),
for i=1:q
   x = x/max(x,[],'all');
   s = s + f(x);
   x = x_d - s;
   vSSIM(i) = ssim(normarr(abs(x)),true_im);
   vPSNR(i) = psnr(normarr(abs(x)),true_im);
end

figure();plot(vSSIM);title("SSIM")
figure();plot(vPSNR);title("PSNR")

end

