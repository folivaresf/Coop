% P1
clear all

% Condiciones Iniciales
x0=1;
y0=2;
% EDOs
% DSOLVE
syms x(t) y(t)
edo1 = diff(x) == 2*x - (1/2)*x^2 - x.*y;
edo2 = diff(y) == -y + x*y;
edo = [edo1; edo2];
[xsol, ysol] = dsolve(edo);

% ODE45
f = @(t,X) [(2*X(1) - (1/2)*X(1)^2 - X(1)*X(2)); (-X(2) + X(1)*X(2))]; % derivada de x e y
[t,S] = ode45(f,[0,2], [x0,y0]);

% GRAFICO
figure
subplot(4,2,1:4)
plot(t,S(:,1))
title('x(t)')
xlabel('t')
ylabel('x(t)')

subplot(4,2,5:8)
plot(t,S(:,2))
title('y(t)')
xlabel('t')
ylabel('y(t)')