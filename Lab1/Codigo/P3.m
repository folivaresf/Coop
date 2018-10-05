% P3
clear all

% Datos
A = [-1 1 1; 2 1 5; 7 3 0];
b = [13; 37; 20];
lb = [0, 0, 0]; ub = [];

f = [-2, -7, -3];
optionsIP = optimoptions('linprog','Algorithm','interior-point');
optionsDS = optimoptions('linprog','Algorithm','dual-simplex');

% Resolver Punto-Interior y Dual-Simplex
[xIP, optIP] = linprog(f, A, b, [], [], lb, ub, optionsIP);
[xDS, optDS] = linprog(f, A, b, [], [], lb, ub, optionsDS);

% Comparar con errores porcentuales
x_error = abs(xIP-xDS)/xIP;
opt_error = abs(optIP-optDS)/optIP;