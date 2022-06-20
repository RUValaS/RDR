function Ad = rnA(A,F,p)
%RNA déforme la matrice d'entrée A

Ad=A; % default si mode incohérent

pente = 0.1;
GausSize = 3; % taille matrice gaussienne (eg 3 = 3x3)
GaussCentr = 0.9; % poids centre 

% génération champ de distortion


N =size(A); % size in pixels of image
       % frequency-filter width

[X,Y] = ndgrid(1:N(1),1:N(2));
i = min(X-1,N(1)-X+1);
j = min(Y-1,N(2)-Y+1);
H = exp(-.5*(i.^2+j.^2)/F^2);
Z = real(ifft2(H.*fft2(randn(N))));

surf(Z);

Ad = p*Z + A;

end

