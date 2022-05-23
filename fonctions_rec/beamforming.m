% BEAMFORMING
%
% on cherche le filtre w tq
% y = w^H * R * w
% est maximisé
%
% la manière la plus simple est avec w = A 
% au facteur près pour normaliser
%
% dans le cas du BF -> y = a^H * R * a
% or R = x*x^H/N 
% où N peut être négligé car constant
% donc y = (a^H * x)(x^H * a)
% où chaque facteur représente la corrélation entre signaux
% pour une direction donnée
%
% -> méthode la plus basique
% voir ensuite les méthodes MVDR (capon) et AAR
%
% défauts : lobes secondaires
% cf PPT pour les définitions des fitres
% -> a implémenter pour faire les différences
% cf Bible pour les défauts de chaque méthode

function y = beamforming(A, R,Mx,My)
y = zeros(Mx*My, 1);
for dir=1:(Mx*My)
    Ad = A(:, dir);
    y(dir) = (Ad')*R*Ad;
end

