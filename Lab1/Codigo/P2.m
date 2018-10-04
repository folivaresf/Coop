% P2
clear all

% Constantes
A=[0 1; -3 0]; B=[0; 2]; 

% Condiciones Iniciales
x0=[0;0];

% Definir u
u = @(t) 1; % u Caso 1
% u = @(t) t; % u Caso 2
% u = @(t) sin(t); % u Caso 3

% EDO
f = @(t,X) A*[X(1), X(2)]' + B*u(t);

% ODE45
[t,S] = ode45(f,[0,10], x0);

% Variacion de Parametros
eAt = @(t) expm(A*t);
syms s h
g = eAt(h)*inv(eAt(s))*B*u(s);
varpar = eAt(h)*x0 + int(g,s, [0,h]);
var1 = matlabFunction(varpar(1));
var2 = matlabFunction(varpar(2));

% GRAFICOS 
% ODE45
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
title('Trayectoria numerica')
xlabel('x(t)')
ylabel('y(t)')

% Var de Parametros
time = 0:0.1:10;
figure
plot(var1(time),var2(time))
title('Trayectoria analitica')
xlabel('x(t)')
ylabel('y(t)')

% Comparacion
figure
hold on
plot(S(:,1),S(:,2),'r')
plot(var1(time),var2(time),'b')
title('Comparacion')
xlabel('x(t)')
ylabel('y(t)')
legend(["bla","blabla"])
hold off


