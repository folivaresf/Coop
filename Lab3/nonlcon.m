function [c, ceq] = nonlcon(z)
    x0 = 1;
    v0 = 0;
    k = 1;
    a = 2;
    
    % Cambiar valor de N
%     N = 2000;
    % N = 1000;
    % N = 500;
    % N = 300;
    % N = 200;
    N = 100;
    
    h = z(1)/N;
    x(1) = x0;
    v(1) = v0;
    
    for i = 1:1:N
        x(i+1) = x(i) + h*v(i);
        v(i+1) = v(i) + h*(-k*x(i) - a*v(i) + z(i+1));
    end
    
    c = [];
    ceq = [x(N), v(N)];
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    