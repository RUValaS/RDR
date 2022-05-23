function img = adaptImg(s, mx)

%jpg img -> s= m*n*3 array
%to mx-1 * my-1 greyscaled & normalized array

%resize
sr = imresize(s, [mx NaN]); %s_resized
% reshape(img, Nx.Ny, 1)
% -> vectorise l'image

%Greyscaling and [0;1]
max =0;
img = zeros(size(sr, 1), size(sr, 2));
for i=1:size(sr, 1)
    for j=1:size(sr, 2)
        m = mean(sr(i,j,:));
        if m>max
            max=m;
        end
        img(i,j) = m;
    end
end
img=img/max;