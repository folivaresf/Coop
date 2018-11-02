clear all
% Parametros
k = 1;
a = 2;
%Cambiar valor de N
% N = 2000;
% N = 1000;
% N = 500;
% N = 300;
% N = 200;
N = 100;

x0 = [1, zeros(1,N)];
% x0 = [1, -ones(1,N)];
% x0 = [1, ones(1,N)];
% x0 = [-1, zeros(1,N)];

lb = [-Inf, -ones(1,N)];
ub = [Inf, ones(1,N)];
fun = @(z) z(1);

[z, fval] = fmincon(fun, x0, [], [], [], [],...
                    lb, ub, @(z)nonlcon(z));

h = z(1)/N;
x(1) = 1;
v(1) = 0;
    
for i = 1:1:N
    x(i+1) = x(i) + h*v(i);
    v(i+1) = v(i) + h*(-k*x(i) - a*v(i) + z(i+1));
end
tiempo = 0:fval/(N-1):fval;
% Graficos
figure(1)
plot(tiempo, z(2:N+1))
title("Control u(t) con N = "+int2str(N))
ylabel('u(t)')
xlabel('t')
saveas(gcf,'C:\Users\Felipe\Desktop\Coop backup\P2_u_100.png')

figure(2)
plot(tiempo, x(1:N))
title("Posicion x(t) con N = "+int2str(N))
ylabel('x(t)')
xlabel('t')
saveas(gcf,'C:\Users\Felipe\Desktop\Coop backup\P2_x_100.png')

figure(3)
plot(tiempo, v(1:N))
title("Velocidad v(t) con N = "+int2str(N))
ylabel('v(t)')
xlabel('t')
saveas(gcf,'C:\Users\Felipe\Desktop\Coop backup\P2_v_100.png')


