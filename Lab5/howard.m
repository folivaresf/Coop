% P6
function [W, U] = howard(V0, u0, X, V)
    W = V0;
    U = u0;
    dif = 1;
    S1 = size(X);
    [f1,f2,l]= funcion([X V],uo);
    Q = Qfuncion(f1,f2);
    A = Q - lambda*eye(S1(1)*S1(2));
    S2 = size(Q);
    b = ones(S2(2));
    tol = 1e-2;
    while(dif > tol) 
       [f1,f2,l]= funcion([X V],U);
       Q = Qfuncion(f1,f2);
       A = Q - lambda*eye(S1(1)*S1(2));
       W1 = A\b;
       U = argmin(A*W1);
       dif = norm(W1 - W);
       W = W1;
    end
end