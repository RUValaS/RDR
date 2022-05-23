function x = MEM(x_d,q)
    function y = f(x)
        y = -x.*log(x);
    end

s = zeros(size(x_d));
x = x_d;
for i=1:q
    x = x/max(x,[],'all');
    s = s + f(x);
    x = x_d - s;
end
end

