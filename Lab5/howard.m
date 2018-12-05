% P6
function [V, u] = howard(V0, u0)
    W(1) = V0;
    v(1) = u0;
    dif = 1;
    k = 0;
    
    tol = 1e-2;
    while(dif > tol)
       W(k+1) = A(v(k))\b;
       f = @(v) A(v)*W(k+1);
       v(k+1) = fminsearch(f, 1);
       dif = norm(W(k+1) - W(k));
       k = k+1;
    end
    V = W(k);
    u = v(k);
end