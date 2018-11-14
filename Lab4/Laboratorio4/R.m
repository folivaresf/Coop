function pT = R(p0)
T = 3.5;
x0 = 2;
tspan = [0 T];
z0 = [x0;p0]; % condiciones iniciales
[t,y] = ode45(@F,tspan,z0); % resolver EDO
p = y(:,2); % p es la segunda coordenada
pT = p(length(t)); % instante final
end

