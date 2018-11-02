function [A] = unvec(v)
    A = [];
    s = size(v);
    n = sqrt(s(1));
    for i=1:n
        A = [A, v((i-1)*n + 1: i*n)];
    end
end