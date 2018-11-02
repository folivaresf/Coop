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

T = 1.5;
%Cambiar valor de N
N = 100;

x0 = zeros(16,1);
tspan = [0, T];

f = @(t,E) -(vec(W) - vec(A'*unvec(E)) - vec(unvec(E)*A) - vec(unvec(E)*(B*B')*unvec(E)));

%ODE45
[t,x] = ode45(f, tspan, x0);

% E = unvec(x(77,:)');

E = @(t2) unvec(x(floor(1+(1.5-t2)*76/1.5),:)');
% E = @(t2) unvec(x(1 + floor(76*t2/T),:)');

g0 = [0.9050, 0.4405, 1.2827, -7.465];
g = @(t2, x2) A*x2 + B*B'*E(t2)*x2;

[t2, x2] = ode45(g, t, g0);

for i=1:77
    u(i) = B'*E(t2(i))*x2(i,:)';
end

% Graficos
figure(1)
subplot(2,2,1)
plot(t2, x2(:,1))
title('x1(t) vs tiempo')
ylabel('x1(t)')
xlabel('t')

subplot(2,2,2)
plot(t2, x2(:,2))
title('x2(t) vs tiempo')
ylabel('x2(t)')
xlabel('t')

subplot(2,2,3)
plot(t2, x2(:,3))
title('x3(t) vs tiempo')
ylabel('x3(t)')
xlabel('t')

subplot(2,2,4)
plot(t2, x2(:,4))
title('x4(t) vs tiempo')
ylabel('x4(t)')
xlabel('t')

figure(2)
plot(t2,u)
title('control optimo vs tiempo')
ylabel('u(t)')
xlabel('t')




