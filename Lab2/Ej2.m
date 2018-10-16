
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

clear m M g f l
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

plot3(x(:,1),x(:,3),t)
grid on
title('Solución con u constante (u(t)=-5)')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_1.png')

[t,x] = ode45(f2,T,X0);

plot3(x(:,1),x(:,3),t)
grid on
title('Solución con u lineal (u(t)=t)')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_2.png')

[t,x] = ode45(f3,T,X0);

plot3(x(:,1),x(:,3),t)
grid on
title('Solución con u tipo Sinusoidal (u(t)=sin(t))')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_3.png')

[t,x] = ode45(f4,T,X0);

plot3(x(:,1),x(:,3),t)
grid on
title('Solución con u tipo Bang-Bang (u(t)=5*(t>5))')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_4.png')

[t,x] = ode45(f5,T,X0);

plot3(x(:,1),x(:,3),t)
grid on
title('Solución con u tipo Feedback')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej2_5.png')

clear f1 f2 f3 f4 f5 x t K

%% Ejercicio 3
%X0 = [0 ]

KAB1 = round([B A*B A^2*B A^3*B],4);
KAB2 = round(ctrb(A,B),4);

if sum(sum(KAB1 == KAB2))==16; disp("La matriz de Kalman calculada y la de la función 'ctrb' son iguales"); end

clear KAB1 KAB2
%% Ejercicio 4

% y = (x1 x3) = C*X
C = [1 0 0 0; 0 1 0 0];

O1 = round([C ; C*A ; C*A^2; C*A^3],4);
O2 = round(obsv(A,C),4);

if sum(sum(O1 == O2))==16; disp("La matriz de observación calculada y la de la función 'obsv' son iguales"); end

clear O1 O2
%% Ejercicio 5

p = poly(A);
a = flip(p);
D = 0;

B1 =  round([zeros(n-1,1) eye(n-1); -a(1:n)],4);

sys = ss(A,B,C,D);

canon = canon(sys,'companion');

B2 = round(canon.A',4);

if sum(sum(B1 == B2))==16; disp("La matriz de Brunovski calculada y la de la función 'canon' son iguales"); end

clear p a D B1 B2 sys canon
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

plot3(x(:,1),x(:,3),t)
grid on
title('Solución con u tipo Feedback usando place')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej6_1.png')

[t,x] = ode45(f7,T,X0);

plot3(x(:,1),x(:,3),t)
grid on
title('Solución con u tipo Feedback usando lqr')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej6_2.png')

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

plot3(x(:,1),x(:,3),t)
grid on
title('Solución x con control constante')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_1.png')

plot3(x(:,5),x(:,7),t)
grid on
title('Solución a x gorro con control constante')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_2.png')


[t,x] = ode45(h2,T,X1);

plot3(x(:,1),x(:,3),t)
grid on
title('Solución x con control lineal')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_3.png')

plot3(x(:,5),x(:,7),t)
grid on
title('Solución a x gorro con control lineal')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_4.png')


[t,x] = ode45(h3,T,X1);

plot3(x(:,1),x(:,3),t)
grid on
title('Solución x con control sinusoidal')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_5.png')

plot3(x(:,5),x(:,7),t)
grid on
title('Solución a x gorro con control sinusoidal')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_6.png')

[t,x] = ode45(h4,T,X1);

plot3(x(:,1),x(:,3),t)
grid on
title('Solución x con control Bang Bang')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_7.png')

plot3(x(:,5),x(:,7),t)
grid on
title('Solución a x gorro con control Bang Bang')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')
saveas(gcf,'C:\Users\usuario\Desktop\Control óptimo\Lab 2\ej7_8.png')

clear h1 h2 h3 h4 u1 u2 u3 u4

%% Ejercicio 8

%D = B*K1;

T=[0 40];
X1 =[0 1 1 0 1000 400 4 1];

g1 = @(t,X) [A*[X(1), X(2), X(3), X(4)]' - B*K1*[X(5),X(6),X(7),X(8)]';...
             A*[X(5), X(6), X(7), X(8)]' - B*K1*[X(5),X(6),X(7),X(8)]' + L*(C*[X(5),X(6),X(7),X(8)]' - C*[X(1),X(2),X(3),X(4)]')];

g2 = @(t,X) [A*[X(1), X(2), X(3), X(4)]' - B*K2*[X(5),X(6),X(7),X(8)]';...
             A*[X(5), X(6), X(7), X(8)]' - B*K2*[X(5),X(6),X(7),X(8)]' + L*(C*[X(5),X(6),X(7),X(8)]' - C*[X(1),X(2),X(3),X(4)]')];
         
[t,x] = ode45(g1,T,X1);

hold on
plot(t,x(:,1))
plot(t,x(:,5))
hold off

plot3(x(:,1),x(:,3),t)
grid on
title('Solución a x con control Feedback')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')

plot3(x(:,5),x(:,7),t)
grid on
title('Solución a x gorro con control Feedback')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')


%% Ejercicio 9
M = 2.4;
m = 0.23;
g = 9.8;
f = 0.1;
l = 0.36;
T = [0 10];

% X(1)  = x
% X(2)  = x'
% X(3)  = \theta
% X(4)  = \theta'

f1 = @(x,y) -m*g*sin(x)*cos(x) + m*l*y^2*sin(x) + m*f*y*cos(x);
f2 = @(x)    M + m*sin(x)^2;
f3 = @(x,y) (M+m)*(g*sin(x)-f*y) - m*l*y^2*sin(x)*cos(x);
u1  = @(x,y,z,w) K1*[x, y, z, w]';
u2  = @(x,y,z,w) K2*[x, y, z, w]';

g3 = @(t,X) [ X(2)                                                               ;...
             (f1(X(3),X(4))+ u1(X(1), X(2), X(3), X(4)))/(f2(X(3))  )            ;...
              X(4)                                                               ;...
             (f3(X(3),X(4))+ u1(X(1), X(2), X(3), X(4))*cos(X(3)))/(f2(X(3))*l)] ;

g4 = @(t,X) [ X(2)                                                               ;...
             (f1(X(3),X(4))+ u2(X(1), X(2), X(3), X(4)))/(f2(X(3))  )            ;...
              X(4)                                                               ;...
             (f3(X(3),X(4))+ u2(X(1), X(2), X(3), X(4))*cos(X(3)))/(f2(X(3))*l)] ;

[t,x] = ode45(g3,T,X0);

plot3(x(:,1),x(:,3),t)
grid on
title('Solución a x con control Feedback')
xlabel('x(t)')
ylabel('\theta(t)')
zlabel('t')


%% weás que no apañaron

% [V] = odeToVectorField(diff(x1,t) == (x2                                     -  D(1,:)*[z1 z2 z3 z4]'),...
%                        diff(x2,t) == (-m*g/M*x3 + f*m/M*x4                   -  D(2,:)*[z1 z2 z3 z4]'),...
%                        diff(x3,t) == (x4                                     -  D(3,:)*[z1 z2 z3 z4]'),...
%                        diff(x4,t) == ( -(m+M)*g/(M*l)*x3 + f*(m+M)/(M*l)*x4  -  D(4,:)*[z1 z2 z3 z4]'),...
%                        diff(z1,t) == (z2                                     -  D(1,:)*[z1 z2 z3 z4]' + L(1,:)*[z1-x1; z3-x3]),...
%                        diff(z2,t) == (-m*g/M*z3 + f*m/M*z4                   -  D(2,:)*[z1 z2 z3 z4]' + L(2,:)*[z1-x1; z3-x3]),...
%                        diff(z3,t) == (z4                                     -  D(3,:)*[z1 z2 z3 z4]' + L(3,:)*[z1-x1; z3-x3]),...
%                        diff(z4,t) == ( -(m+M)*g/(M*l)*z3 + f*(m+M)/(M*l)*z4  -  D(4,:)*[z1 z2 z3 z4]' + L(4,:)*[z1-x1; z3-x3]) );
