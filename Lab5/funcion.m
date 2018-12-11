function [f1,f2,l] = funcion(W,U)
S = size(W);
l = ones(S);
for i = 1:S(1)
    for j = 1:S(1)
        k = sub2ind([S(1), S(1)], i, j);
        f1(i,j) = W(i,j+S(1));
        f2(i,j) = U(k);
    end
end
end
