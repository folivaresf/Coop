function [f1,f2,l] = funcion(W,U)
S = size(W);
l = ones(S);
for i = 1:S(1)
    for j = 1:S(1)
        f1(i,j) = W(i,j+S(1));
        f2(i,j) = U(i,j);
    end
end
end
