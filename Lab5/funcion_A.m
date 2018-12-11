%funcion A
function A = funcion_A(v, X, V, lambda)
    [f1,f2,l]= funcion([X V], v);
    Q = Qfuncion(f1,f2);
    S1 = size(X);
    S2 = size(Q);
    A = Q - lambda*eye(S1(1)*S1(2));
end