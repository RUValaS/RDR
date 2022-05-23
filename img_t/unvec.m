function img = unvec(vimg, nx, ny)

if numel(vimg)~=nx*ny
    ME = MException('myComponent:inputError', 'mauvaises dimensions');
    throw(ME)
end

img = zeros(nx, ny);
for i=1:nx
    for j=1:ny
        img(i,j) = vimg((i-1)*ny+j);
    end
end
