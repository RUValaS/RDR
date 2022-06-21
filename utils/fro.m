function fr = fro(A)
%FRO Calcul norme frobenius
fr = trace(A*(A'))/numel(A);
end

