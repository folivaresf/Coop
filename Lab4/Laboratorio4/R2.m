function devol= R2(p0)
% p0 de dim 3
% correr la solucion con esas condiciones iniciales
x0 = 2; %variable
y0 = 0; %fijo
z0 = 5; %variable
beta = 1; %variable
T = 10; % FIJO
tspan = [0 T];
w0 = [x0;y0;z0;p0(1);p0(2);p0(3)]; % condiciones iniciales
odeopts = odeset('NonNegative', [1, 2, 3, 4]);
[t,w] = ode45(@F2,tspan,w0,odeopts); % resolver EDO
% ARREGLAR
%adjunto solucion
px = w(:,4);
py = w(:,5);
pz = w(:,6); 
% instante final
pxT = px(length(t)); 
pyT = py(length(t)); 
pzT = pz(length(t));
devol = [pxT-beta;pyT+1;pzT];
end

