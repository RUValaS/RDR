function res = amel(avant,after,ref)
%AMEL Summary of this function goes here
%   Detailed explanation goes here
res = fro(avant-ref) - fro(after-ref);
end

