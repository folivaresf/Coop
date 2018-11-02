function [v] = vec(A)
    v = [];
    s = size(A);
    for i=1:s(2)
        v = [v; A(:,i)];
    end
end