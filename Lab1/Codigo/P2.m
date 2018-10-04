% P2
clear all

% Constantes
A=[0 1; -3 0]; B=[0; 2]; 

% Condiciones Iniciales
x0=[0;0];

% Definir u
u = @(t) 1; % u Caso 1
% u = @(t) t % u Caso 2
% u = @(t) sin(t) % u Caso 3

% EDO
f = @(t,X) A*[X(1), X(2)]' + B*u(t);

[t,S] = ode45(f,[0,10], x0);

% GRAFICO
figure
subplot(2,1,1)
plot(t,S(:,1))
title('x(t)')
xlabel('t')
ylabel('x(t)')

subplot(2,1,2)
plot(t,S(:,2))
title('y(t)')
xlabel('t')
ylabel('y(t)')

figure
plot(S(:,1),S(:,2))
title('plooot')
xlabel('x(t)')
ylabel('y(t)')
