function I = matI(Mx, My)

Q = Mx*My;
I = zeros(Q, 2);

X = linspace(-1,1,Mx);
Y = linspace(-1,1,My);

for i=1:Mx
    for j=1:My
        I((i-1)*My + j, :) = [X(i) Y(j)];
    end
end

