%Matriz A(u)
function [A] = A(u)
    n = 100;
    nx = n;
    nv = n;
    hx = 0.01;
    hv = 0.01;
    
    N = (nx+1)*(nv+1);
    A = zeros(N, N);
    [f1, f2] = funcion(u);
    
    for i = 1:N
        for j = 1:N
            if i==j
                A(i, j) = lambda - abs(f1/hx) - abs(f2/hv);
            end
            if j==i-1
                A(i, j) = (-f1/hx)*(f1<0);
            end
            if j==i+1
                A(i, j) = (f1/hx)*(f1>=0);
            end
            if j==i-(nx+1)
                A(i, j) = (-f2/hv)*(f2<0);
            end
            if j==i+(nx+1)
                A(i, j) = (f2/hv)*(f2>=0);
            end
        end
    end
end