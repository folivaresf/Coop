clear all
% Parametros
m = 0.23;
M = 2.4;
g = 9.8;
f = 0.1;
l = 0.36;

A = [0 1 0 0;...
     0 0 -m*g/M f*m/M;...
     0 0 0 1;...
     0 0 (M+m)*g/(M*l) -f*(m+M)/(M*l)];
 
B = [0; 1/M; 0; -1/(M*l)];
W = [1, 0, 0, 0;...
     0, 0, 0, 0;...
     0, 0, 1, 0;...
     0, 0, 0, 0];

H = -[A, B*transpose(B);...
     W, -transpose(A)];

T = 1.5;
%Cambiar valor de N
N = 100;

x0 = [1, 0, 0, 0, 0, 0, 0, 0];
tspan = [0, T];

f = @(t,x) H*x;

%ODE45
[t,x] = ode45(f, tspan, x0);


% Graficos
figure(1)
subplot(2,2,1)
plot(t, x(:,1))
title('x1(t) vs tiempo')
ylabel('x1(t)')
xlabel('t')

subplot(2,2,2)
plot(t, x(:,2))
title('x2(t) vs tiempo')
ylabel('x2(t)')
xlabel('t')

subplot(2,2,3)
plot(t, x(:,3))
title('x3(t) vs tiempo')
ylabel('x3(t)')
xlabel('t')

subplot(2,2,4)
plot(t, x(:,4))
title('x4(t) vs tiempo')
ylabel('x4(t)')
xlabel('t')







