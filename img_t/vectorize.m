function vimg = vectorize(img)

nb = numel(img);
vimg = zeros(nb,1);
nx = size(img, 1);
ny = size(img, 2);

for i=1:nx
    for j=1:ny
        vimg((i-1)*ny+j)=img(i,j);
    end
end
