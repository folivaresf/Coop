clear clc
%% 
% constantes
n = 4;
j = 1;

M = 2.4;
m = 0.23;
g = 9.8;
f = 0.1;
l = 0.36;

A = [0 1 0 0; 0 0 -((m*g)/M) ((f*m)/M); 0 0 0 1; 0 0 (((M+m)*g)/(M*l)) -((f*(m+M))/(M*l))];
B = [0 1/M 0 -1/M]';

T = [0 10];
X0 = [0 0 1 0];

clear m M g f l;
%% Ejercicio 2
% controles
u1 = @(t) -5; % constante
u2 = @(t) t; % lineal
u3 = @(t) sin(t); % sinusoidal
u4 = @(t) 5*(t>5); % bang bang
K = [-2  -4 -160 -30]; % para el control feedback
% no encontré controles de tipo etc.

f1 = @(t,X) A*[X(1), X(2), X(3), X(4)]' + B*u1(t);
f2 = @(t,X) A*[X(1), X(2), X(3), X(4)]' + B*u2(t);
f3 = @(t,X) A*[X(1), X(2), X(3), X(4)]' + B*u3(t);
f4 = @(t,X) A*[X(1), X(2), X(3), X(4)]' + B*u4(t);
f5 = @(t,X) (A-B*K)*[X(1), X(2), X(3), X(4)]' ; %Feedback

[t,x] = ode45(f1,T,X0);

plot(t,x(:,1),'LineWidth',3)
title({'Solución con u constante','u(t)=-5'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_01.png')

plot(t,x(:,3),'LineWidth',3)
title({'Solución con u constante','u(t)=-5'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_02.png')

plot3(x(:,1),x(:,3),t,'LineWidth',3)
grid on
title({'Solución con u constante','u(t)=-5'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_03.png')

[t,x] = ode45(f2,T,X0);

plot(t,x(:,1),'LineWidth',3)
title({'Solución con u lineal','u(t)=t'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_04.png')

plot(t,x(:,3),'LineWidth',3)
title({'Solución con u lineal','u(t)=t'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_05.png')

plot3(x(:,1),x(:,3),t,'LineWidth',3)
grid on
title({'Solución con u lineal','u(t)=t'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_06.png')

[t,x] = ode45(f3,T,X0);

plot(t,x(:,1),'LineWidth',3)
title({'Solución con u sinusoidal','u(t)=sin(t)'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_07.png')

plot(t,x(:,3),'LineWidth',3)
title({'Solución con u sinusoidal','u(t)=sin(t)'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_08.png')

plot3(x(:,1),x(:,3),t,'LineWidth',3)
grid on
title({'Solución con u sinusoidal','u(t)=sin(t)'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_09.png')

[t,x] = ode45(f4,T,X0);

plot(t,x(:,1),'LineWidth',3)
title({'Solución con u Bang-Bang','u(t)=5*(t>5)'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_10.png')

plot(t,x(:,3),'LineWidth',3)
title({'Solución con u Bang-Bang','u(t)=5*(t>5)'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_11.png')

plot3(x(:,1),x(:,3),t,'LineWidth',3)
grid on
title({'Solución con u Bang-Bang','u(t)=5*(t>5)'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_12.png')

[t,x] = ode45(f5,T,X0);

plot(t,x(:,1),'LineWidth',3)
title({'Solución con u Feedback lineal','u(t)=-Kx(t)'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_13.png')

plot(t,x(:,3),'LineWidth',3)
title({'Solución con u Feedback lineal','u(t)=-Kx(t)'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_14.png')

plot3(x(:,1),x(:,3),t,'LineWidth',3)
grid on
title({'Solución con u Feedback lineal','u(t)=-Kx(t)'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_15.png')

clear f1 f2 f3 f4 f5 x t K;

%% Ejercicio 3
%X0 = [0 ]

KAB1 = round([B A*B A^2*B A^3*B],4);
KAB2 = round(ctrb(A,B),4);

if sum(sum(KAB1 == KAB2))==16; disp("La matriz de Kalman calculada y la de la función 'ctrb' son iguales"); end

clear KAB1 KAB2;
%% Ejercicio 4

% y = (x1 x3) = C*X
C = [1 0 0 0; 0 0 1 0];

O1 = round([C ; C*A ; C*A^2; C*A^3],4);
O2 = round(obsv(A,C),4);

if sum(sum(O1 == O2))==32; disp("La matriz de observación calculada y la de la función 'obsv' son iguales"); end

clear O1 O2;
%% Ejercicio 5

p = poly(A);
a = flip(p);
D = 0;

B1 =  round([zeros(n-1,1) eye(n-1); -a(1:n)],4);

sys = ss(A,B,C,D);

canon = canon(sys,'companion');

B2 = round(canon.A',4);

if sum(sum(B1 == B2))==16; disp("La matriz de Brunovski calculada y la de la función 'canon' son iguales"); end

clear p a D B1 B2 sys canon;
%% Ejercicio 6
p  = [-1 -2 -3 -4]';
K1 = place(A,B,p);

Q = eye(n);
R = eye(j);
N = 0;
[K2,~,~] = lqr(A,B,Q,R,N);

lambda1 = eig(A-B*K1); % <0
lambda2 = eig(A-B*K2); % <0

if sum(lambda1<0)==4; disp("los v.p. de la matriz obtenida con 'place' son todos negativos"); end
if sum(lambda2<0)==4; disp("los v.p. de la matriz obtenida con 'lqr' son todos negativos"); end
% son estabilizadores
clear lambda1 lambda2 Q R N p

f6 = @(t,X) (A-B*K1)*[X(1), X(2), X(3), X(4)]' ; %Feedback1
f7 = @(t,X) (A-B*K2)*[X(1), X(2), X(3), X(4)]' ; %Feedback1

[t,x] = ode45(f6,T,X0);

plot(t,x(:,1),'LineWidth',3)
title({'Solución con u Feedback lineal','u(t)=-Kx(t)','obtenido usando place'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej6_01.png')

plot(t,x(:,3),'LineWidth',3)
title({'Solución con u Feedback lineal','u(t)=-Kx(t)','obtenido usando place'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej6_02.png')

plot3(x(:,1),x(:,3),t,'LineWidth',3)
grid on
title({'Solución con u Feedback lineal','u(t)=-Kx(t)','obtenido usando place'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej6_03.png')

[t,x] = ode45(f7,T,X0);

plot(t,x(:,1),'LineWidth',3)
title({'Solución con u Feedback lineal','u(t)=-Kx(t)','obtenido usando lqr'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej6_04.png')

plot(t,x(:,3),'LineWidth',3)
title({'Solución con u Feedback lineal','u(t)=-Kx(t)','obtenido usando lqr'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej6_05.png')

plot3(x(:,1),x(:,3),t,'LineWidth',3)
grid on
title({'Solución con u Feedback lineal','u(t)=-Kx(t)','obtenido usando lqr'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej6_06.png')

clear f6 f7 x t
%% Ejercicio 7
p  = [-1 -2 -3 -4]';
K3 = place(A',-C',p);
L = K3';

eig(A+L*C)


%h1 = @(t,X) (AA)*[X(1), X(2), X(3), X(4),X(5), X(6),X(7),X(8)]' + BB*u1(t)+ LL*C*([X(5),X(6),X(7),X(8)]' - [X(1),X(2),X(3),X(4)]');
%h2 = @(t,X) (AA)*[X(1), X(2), X(3), X(4),X(5), X(6),X(7),X(8)]' + BB*u2(t)+ LL*C*([X(5),X(6),X(7),X(8)]' - [X(1),X(2),X(3),X(4)]');
%h3 = @(t,X) (AA)*[X(1), X(2), X(3), X(4),X(5), X(6),X(7),X(8)]' + BB*u3(t)+ LL*C*([X(5),X(6),X(7),X(8)]' - [X(1),X(2),X(3),X(4)]');
%h4 = @(t,X) (AA)*[X(1), X(2), X(3), X(4),X(5), X(6),X(7),X(8)]' + BB*u4(t)+ LL*C*([X(5),X(6),X(7),X(8)]' - [X(1),X(2),X(3),X(4)]');

h1 = @(t,X) [A*[X(1), X(2), X(3), X(4)]' + B*u1(t);...
             A*[X(5), X(6), X(7), X(8)]' + B*u1(t) + L*(C*[X(5),X(6),X(7),X(8)]' - C*[X(1),X(2),X(3),X(4)]')];
h2 = @(t,X) [A*[X(1), X(2), X(3), X(4)]' + B*u2(t);...
             A*[X(5), X(6), X(7), X(8)]' + B*u2(t) + L*(C*[X(5),X(6),X(7),X(8)]' - C*[X(1),X(2),X(3),X(4)]')];
h3 = @(t,X) [A*[X(1), X(2), X(3), X(4)]' + B*u3(t);...
             A*[X(5), X(6), X(7), X(8)]' + B*u3(t) + L*(C*[X(5),X(6),X(7),X(8)]' - C*[X(1),X(2),X(3),X(4)]')];
h4 = @(t,X) [A*[X(1), X(2), X(3), X(4)]' + B*u4(t);...
             A*[X(5), X(6), X(7), X(8)]' + B*u4(t) + L*(C*[X(5),X(6),X(7),X(8)]' - C*[X(1),X(2),X(3),X(4)]')];

X1 = [0 1 1 0 1000 4 4 1];

[t,x] = ode45(h1,T,X1);

% x constante

plot(t,x(:,1),'LineWidth',3)
title({'Solución de x con','u constante, u(t)=-5'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_01.png')

plot(t,x(:,3),'LineWidth',3)
title({'Solución de \theta con','u constante, u(t)=-5'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_02.png')

plot3(x(:,1),x(:,3),t,'LineWidth',3)
grid on
title({'Solución de (x,\theta) con','u constante, u(t)=-5'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_03.png')

% x gorro constante

plot(t,x(:,5),'LineWidth',3)
title({'Solución de x gorro con','u constante, u(t)=-5'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_04.png')

plot(t,x(:,7),'LineWidth',3)
title({'Solución de \theta gorro con','u constante, u(t)=-5'})
ylabel('\theta(t)')
xlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_05.png')

plot3(x(:,5),x(:,7),t,'LineWidth',3)
grid on
title({'Solución de (x,\theta) gorro con','u constante, u(t)=-5'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_6.png')


[t,x] = ode45(h2,T,X1);

%  x lineal

plot(t,x(:,1),'LineWidth',3)
title({'Solución de x con','u lineal, u(t)=t'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_07.png')

plot(t,x(:,3),'LineWidth',3)
title({'Solución de \theta con','u lineal, u(t)=t'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_08.png')

plot3(x(:,1),x(:,3),t,'LineWidth',3)
grid on
title({'Solución de (x,\theta) con','u lineal, u(t)=t'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_09.png')

% x gorro lineal

plot(t,x(:,5),'LineWidth',3)
title({'Solución de x gorro con','u lineal, u(t)=t'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_10.png')

plot(t,x(:,7),'LineWidth',3)
title({'Solución de \theta gorro con','u lineal, u(t)=t'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_11.png')

plot3(x(:,5),x(:,7),t,'LineWidth',3)
grid on
title({'Solución de (x,\theta) gorro con','u lineal, u(t)=t'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_12.png')

[t,x] = ode45(h3,T,X1);

% x sinusoidal

plot(t,x(:,1),'LineWidth',3)
title({'Solución de x con','u sinusoidal, u(t)=sin(t)'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_13.png')

plot(t,x(:,3),'LineWidth',3)
title({'Solución de \theta con','u sinusoidal, u(t)=sin(t)'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_14.png')

plot3(x(:,1),x(:,3),t,'LineWidth',3)
grid on
title({'Solución de (x,\theta) con','u sinusoidal, u(t)=sin(t)'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_15.png')

%X gorro sinusoidal

plot(t,x(:,5),'LineWidth',3)
title({'Solución de x gorro con','u sinusoidal, u(t)=sin(t)'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_16.png')

plot(t,x(:,7),'LineWidth',3)
title({'Solución de \theta gorro con','u sinusoidal, u(t)=sin(t)'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_17.png')

plot3(x(:,5),x(:,7),t,'LineWidth',3)
grid on
title({'Solución de (x,\theta) gorro con','u sinusoidal, u(t)=sin(t)'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_18.png')

[t,x] = ode45(h4,T,X1);

% x bang bang

plot(t,x(:,1),'LineWidth',3)
title({'Solución de x con','u Bang-Bang, u(t)=5*(t>5)'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_19.png')

plot(t,x(:,3),'LineWidth',3)
title({'Solución de \theta con','u Bang-Bang, u(t)=5*(t>5)'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_20.png')

plot3(x(:,1),x(:,3),t,'LineWidth',3)
grid on
title({'Solución de (x,\theta) con','u Bang-Bang, u(t)=5*(t>5)'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_21.png')

%x gorro bang bang

plot(t,x(:,5),'LineWidth',3)
title({'Solución de x gorro con','u Bang-Bang, u(t)=5*(t>5)'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_22.png')

plot(t,x(:,7),'LineWidth',3)
title({'Solución de \theta gorro con','u Bang-Bang, u(t)=5*(t>5)'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_23.png')

plot3(x(:,5),x(:,7),t,'LineWidth',3)
grid on
title({'Solución de (x,\theta) gorro con','u Bang-Bang, u(t)=5*(t>5)'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_24.png')

clear h1 h2 h3 h4 u1 u2 u3 u4

%% Ejercicio 8

%D = B*K1;

T=[0 20];
X1 = [0 1 1 0 1000 4 4 1];

% con Luenberger
g1 = @(t,X) [A*[X(1), X(2), X(3), X(4)]' - B*K1*[X(5),X(6),X(7),X(8)]';...
             A*[X(5), X(6), X(7), X(8)]' - B*K1*[X(5),X(6),X(7),X(8)]' + L*(C*[X(5),X(6),X(7),X(8)]' - C*[X(1),X(2),X(3),X(4)]')];

g2 = @(t,X) [A*[X(1), X(2), X(3), X(4)]' - B*K2*[X(5),X(6),X(7),X(8)]';...
             A*[X(5), X(6), X(7), X(8)]' - B*K2*[X(5),X(6),X(7),X(8)]' + L*(C*[X(5),X(6),X(7),X(8)]' - C*[X(1),X(2),X(3),X(4)]')];

[t,x] = ode45(g1,T,X1);

hold on
plot(t(1:100),x(1:100,1),'LineWidth',3)
plot(t(1:100),x(1:100,5),'LineWidth',3)
title({'Solución de x y x gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando place'})
xlabel('t')
ylabel('x(t)')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej8_01.png')

hold on
plot(t,x(:,1),'LineWidth',3)
plot(t,x(:,5),'LineWidth',3)
title({'Solución de x y x gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando place'})
xlabel('t')
ylabel('x(t)')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej8_02.png')

hold on
plot(t(1:120),x(1:120,3),'LineWidth',3)
plot(t(1:120),x(1:120,7),'LineWidth',3)
title({'Solución de \theta y \theta gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando place'})
xlabel('t')
ylabel('\theta(t)')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej8_03.png')

hold on
plot(t,x(:,3),'LineWidth',3)
plot(t,x(:,7),'LineWidth',3)
ylabel('\theta(t)')
title({'Solución de \theta y \theta gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando place'})
xlabel('t')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej8_04.png')

plot3(x(:,1),x(:,3),t,'LineWidth',3)
grid on
title({'Solución de (x,\theta) con','u Feedback, u(t)=-Kx(t)','obtenido usando place'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej8_05.png')

plot3(x(:,5),x(:,7),t,'LineWidth',3)
grid on
title({'Solución de (x,\theta) gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando place'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej8_06.png')

[t,x] = ode45(g2,T,X1);

hold on
plot(t(1:100),x(1:100,1),'LineWidth',3)
plot(t(1:100),x(1:100,5),'LineWidth',3)
title({'Solución de x y x gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando lqr'})
xlabel('t')
ylabel('x(t)')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej8_07.png')

hold on
plot(t,x(:,1),'LineWidth',3)
plot(t,x(:,5),'LineWidth',3)
title({'Solución de x y x gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando lqr'})
xlabel('t')
ylabel('x(t)')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej8_08.png')

hold on
plot(t(1:120),x(1:120,3),'LineWidth',3)
plot(t(1:120),x(1:120,7),'LineWidth',3)
title({'Solución de \theta y \theta gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando lqr'})
xlabel('t')
ylabel('\theta(t)')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej8_09.png')

hold on
plot(t,x(:,3),'LineWidth',3)
plot(t,x(:,7),'LineWidth',3)
ylabel('\theta(t)')
title({'Solución de \theta y \theta gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando lqr'})
ylabel('t')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej8_10.png')

plot3(x(:,1),x(:,3),t,'LineWidth',3)
grid on
title({'Solución de (x,\theta) con','u Feedback, u(t)=-Kx(t)','obtenido usando lqr'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej8_11.png')

plot3(x(:,5),x(:,7),t,'LineWidth',3)
grid on
title({'Solución de (x,\theta) gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando lqr'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej8_12.png')

clear g1 g2
%% Ejercicio 9
M = 2.4;
m = 0.23;
g = 9.8;
f = 0.1;
l = 0.36;
T = [0 40];

% X(1)  = x
% X(2)  = x'
% X(3)  = \theta
% X(4)  = \theta'

f1 = @(x,y) -m*g*sin(x)*cos(x) + m*l*y^2*sin(x) + m*f*y*cos(x);
f2 = @(x)    M + m*sin(x)^2;
f3 = @(x,y) (M+m)*(g*sin(x)-f*y) - m*l*y^2*sin(x)*cos(x);
u1  = @(x,y,z,w) K1*[x, y, z, w]';
u2  = @(x,y,z,w) K2*[x, y, z, w]';

% sin luemberger

g1 = @(t,X) [ X(2)                                                               ;...
             (f1(X(3),X(4))+ u1(X(1), X(2), X(3), X(4)))/(f2(X(3))  )            ;...
              X(4)                                                               ;...
             (f3(X(3),X(4))+ u1(X(1), X(2), X(3), X(4))*cos(X(3)))/(f2(X(3))*l)] ;

g2 = @(t,X) [ X(2)                                                               ;...
             (f1(X(3),X(4))+ u2(X(1), X(2), X(3), X(4)))/(f2(X(3))  )            ;...
              X(4)                                                               ;...
             (f3(X(3),X(4))+ u2(X(1), X(2), X(3), X(4))*cos(X(3)))/(f2(X(3))*l)] ;

% con luenberger
LL = L*C;

g3 = @(t,X) [ X(2)                                                               ;...
             (f1(X(3),X(4))+ u1(X(5), X(6), X(7), X(8)))/(f2(X(3))  )            ;...
              X(4)                                                               ;...
             (f3(X(3),X(4))+ u1(X(5), X(6), X(7), X(8))*cos(X(3)))/(f2(X(3))*l)  ;...
             ((X(6))                                                               +(LL(1,:)*([X(5),X(6),X(7),X(8)]' - [X(1),X(2),X(3),X(4)]')));...
             (((f1(X(7),X(8))+ u1(X(5), X(6), X(7), X(8)))/(f2(X(7))  ) )           +(LL(2,:)*([X(5),X(6),X(7),X(8)]' - [X(1),X(2),X(3),X(4)]')));...
             ((X(8))                                                               +(LL(3,:)*([X(5),X(6),X(7),X(8)]' - [X(1),X(2),X(3),X(4)]')));...
             (((f3(X(7),X(8))+ u1(X(5), X(6), X(7), X(8))*cos(X(7)))/(f2(X(7))*l)) + (LL(4,:)*([X(5),X(6),X(7),X(8)]' - [X(1),X(2),X(3),X(4)]'))) ];
         
g4 = @(t,X) [ X(2)                                                               ;...
             (f1(X(3),X(4))+ u2(X(5), X(6), X(7), X(8)))/(f2(X(3))  )            ;...
              X(4)                                                               ;...
             (f3(X(3),X(4))+ u2(X(5), X(6), X(7), X(8))*cos(X(3)))/(f2(X(3))*l)  ;...
             ((X(6)                                                             ) +(LL(1,:)*([X(5),X(6),X(7),X(8)]' - [X(1),X(2),X(3),X(4)]')));...
             (((f1(X(7),X(8))+ u2(X(5), X(6), X(7), X(8))          )/(f2(X(7))  ))+(LL(2,:)*([X(5),X(6),X(7),X(8)]' - [X(1),X(2),X(3),X(4)]')));...
             ((X(8)                                                             ) +(LL(3,:)*([X(5),X(6),X(7),X(8)]' - [X(1),X(2),X(3),X(4)]')));...
             (((f3(X(7),X(8))+ u2(X(5), X(6), X(7), X(8))*cos(X(7)))/(f2(X(7))*l))+(LL(4,:)*([X(5),X(6),X(7),X(8)]' - [X(1),X(2),X(3),X(4)]')))];
         
%[t,x] = ode45(g1,T,X0);
%[t,x] = ode45(g2,T,X0);
         
%[t,x] = ode45(g3,T,X1);

%% GRAFICOS FINALES DE G1 Y G2
X0 = [0 1 1 0];
%X0= [0 1 0.1 0];
[t,x] = ode45(g1,T,X0);

plot(t,x(:,1),'LineWidth',2)
title({'Solución con u Feedback lineal','u(t)=-Kx(t)','obtenido usando place'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_01.png')

plot(t,x(:,3),'LineWidth',2)
title({'Solución con u Feedback lineal','u(t)=-Kx(t)','obtenido usando place'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_02.png')

plot3(x(:,1),x(:,3),t,'LineWidth',1)
grid on
title({'Solución con u Feedback lineal','u(t)=-Kx(t)','obtenido usando place'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_03.png')

[t,x] = ode45(g2,T,X0);

plot(t,x(:,1),'LineWidth',2)
title({'Solución con u Feedback lineal','u(t)=-Kx(t)','obtenido usando lqr'})
ylabel('x(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_04.png')

plot(t,x(:,3),'LineWidth',2)
title({'Solución con u Feedback lineal','u(t)=-Kx(t)','obtenido usando lqr'})
ylabel('\theta(t)')
xlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_05.png')

plot3(x(:,1),x(:,3),t,'LineWidth',1)
grid on
title({'Solución con u Feedback lineal','u(t)=-Kx(t)','obtenido usando lqr'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_06.png')

%% GRAFICOS FINALES DE G3 Y G4
%X1 = [0 1 0.1 0 1 0.5 0.1 0.1];
X1 = [0 1 1 0 1000 4 4 1];
[t,x] = ode45(g3,T,X1);

hold on
plot(t(1:400),x(1:400,1),'LineWidth',1)
plot(t(1:400),x(1:400,5),'LineWidth',1)
title({'Solución de x y x gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando place'})
xlabel('t')
ylabel('x(t)')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_07.png')

hold on
plot(t,x(:,1),'LineWidth',1)
plot(t,x(:,5),'LineWidth',1)
title({'Solución de x y x gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando place'})
xlabel('t')
ylabel('x(t)')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_08.png')

hold on
plot(t(1:1000),x(1:1000,3),'LineWidth',1)
plot(t(1:1000),x(1:1000,7),'LineWidth',1)
title({'Solución de \theta y \theta gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando place'})
xlabel('t')
ylabel('\theta(t)')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_09.png')

hold on
plot(t,x(:,3),'LineWidth',1)
plot(t,x(:,7),'LineWidth',1)
ylabel('\theta(t)')
title({'Solución de \theta y \theta gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando place'})
xlabel('t')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_10.png')

plot3(x(:,1),x(:,3),t,'LineWidth',1)
grid on
title({'Solución de (x,\theta) con','u Feedback, u(t)=-Kx(t)','obtenido usando place'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_11.png')

plot3(x(:,5),x(:,7),t,'LineWidth',1)
grid on
title({'Solución de (x,\theta) gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando place'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_12.png')

[t,x] = ode45(g4,T,X1);

hold on
plot(t(1:500),x(1:500,1),'LineWidth',1)
plot(t(1:500),x(1:500,5),'LineWidth',1)
title({'Solución de x y x gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando lqr'})
xlabel('t')
ylabel('x(t)')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_13.png')

hold on
plot(t,x(:,1),'LineWidth',1)
plot(t,x(:,5),'LineWidth',1)
title({'Solución de x y x gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando lqr'})
xlabel('t')
ylabel('x(t)')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_14.png')

hold on
plot(t(1:1000),x(1:1000,3),'LineWidth',1)
plot(t(1:1000),x(1:1000,7),'LineWidth',1)
title({'Solución de \theta y \theta gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando lqr'})
xlabel('t')
ylabel('\theta(t)')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_15.png')

hold on
plot(t,x(:,3),'LineWidth',1)
plot(t,x(:,7),'LineWidth',1)
ylabel('\theta(t)')
title({'Solución de \theta y \theta gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando lqr'})
ylabel('t')
set(gca,'fontsize',13)
hold off
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_16.png')

plot3(x(:,1),x(:,3),t,'LineWidth',1)
grid on
title({'Solución de (x,\theta) con','u Feedback, u(t)=-Kx(t)','obtenido usando lqr'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_17.png')

plot3(x(:,5),x(:,7),t,'LineWidth',1)
grid on
title({'Solución de (x,\theta) gorro con','u Feedback, u(t)=-Kx(t)','obtenido usando lqr'})
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
set(gca,'fontsize',13)
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej9_18.png')


%% weás que no apañaron

