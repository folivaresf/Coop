% P6
function [W, U] = howard(V0, u0, X, V)
    lambda = 100;

    W = V0;
    U = u0;
    dif = 1;
    S1 = size(X);
    [f1,f2,l]= funcion([X V], u0);
    Q = Qfuncion(f1,f2);
    A = Q - lambda*eye(S1(1)*S1(2));
    S2 = size(Q);
    b = ones(S2(2), 1);
    tol = 1e-2;
    while(dif > tol) 
%       [f1,f2,l]= funcion([X V],U);
%       Q = Qfuncion(f1,f2);
%       A = Q - lambda*eye(S1(1)*S1(2));
        A = funcion_A(U, X, V, lambda);
        W1 = A\b;
%       U = argmin(A*W1);
        YES = @(v) funcion_A(v, X, V, lambda)*W1;
        for i = 1:S2(1)
            for j = 1:S2(2)
                k = sub2ind([S2(1), S2(2)], i, j);
%                 U(k) = @(v) subsref(YES(v), struct('type', '()', 'subs', {{k}}));
                if W1(k)<= 0
                    U(k) = 0;
                elseif i <= (S1(1)+1)
                    if W1(k + S1(1) + 1) - W1(k) < 0
                       U(k) = 1;
                    else
                        U(k) = 0;
                    end
                elseif (i > (S1(1)+1)) && (i <= S1(2)-S1(1)+1)
                    if W1(k+S(1)+1) - W1(k-S1(1)-1) > 0
                       if W1(k-S1(1)-1) - W(k) < 0
                          % FALTA MATRACA 
                       end
                    end
                end
                end
            end
        end
        dif = norm(W1 - W);
        W = W1;
    end
end