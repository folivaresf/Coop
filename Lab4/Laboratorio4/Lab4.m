%% Problema 3
% IMPORTANTE!
% es necesario cambiar manualmente T y x0 en R
x0 = 2;
T = 3.5;
% calcular p0
p0 = fzero(@(po) R(po),1);

%Resolver problema
z0 = [x0;p0];
tspan = [0 T];
[t,z] = ode45(F,tspan,z0);

p = z(:,2);
u = @(p) (p<(1/sqrt(3)))*(1-(2*p)/(sqrt(p^2+1)));

control = zeros(length(t),1);
for i = 1:length(t)
    control(i) = u(p(i));
end

% Graficar solución
subplot(2,1,1)
plot(t, z(:,1))
xlabel('tiempo t')
ylabel('solucion x(t)')
title('Contaminación')

subplot(2,1,2)
plot(t,control)
xlabel('tiempo t')
ylabel('solucion u(t) óptimo de Pontryaguin')
title('Fertilizante')

%% Problema 6
% Resolver el problema 
% IMPORTANTE!
% es necesario cambiar manualmente x0 e z0 en R2 tambien!
x0 = 2;
z0 = 5;
T = 10; %fijo
% calcular p0

options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000,...
    'Display', 'iter', 'TolX', 1e-15, 'TolFun', 1e-15);
%p0 = fsolve(@(po) R2(po),[0.3; 0.25; 10],options);
p0=[0.3;0.25;12.005];
w0 = [x0;0;z0;p0];
tspan = [0 T];
odeopts = odeset('NonNegative', [1, 2, 3, 4]);
[t,w] = ode45(@(t,w) F2(t,w),tspan,w0,odeopts);


%Graficar solución
subplot(2,3,1)
plot(t, w(:,1))
xlabel('tiempo t')
ylabel('solucion x1(t)')
title('Contaminación')

subplot(2,3,2)
plot(t, w(:,2))
xlabel('tiempo t')
ylabel('solucion x2(t)')
title('Cereal')

subplot(2,3,3)
plot(t, w(:,3))
xlabel('tiempo t')
ylabel('solucion x3(t)')
title('Insectos')
 
% subplot(5,1,2)
% plot(t,control)
% xlabel('tiempo t')
% ylabel('solucion u(t) óptimo')
% title('Control dado por Pontryaguin')


p = w(:,4:6);

controlU = zeros(length(t),1);
for i = 1:length(t)
    controlU(i) = uopt(p(i,:));
end

controlV = zeros(length(t),1);
for i = 1:length(t)
    controlV(i) = vopt(p(i,:));
end

subplot(2,3,4)
plot(t,controlU)
xlabel('tiempo t')
ylabel('solucion u(t) óptimo')
title('Control u por Pontryaguin')


subplot(2,3,5)
plot(t,controlV)
xlabel('tiempo t')
ylabel('solucion v(t) óptimo')
title('Control v por Pontryaguin')

