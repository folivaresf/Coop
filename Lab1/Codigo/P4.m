% P4
clear all

f = @(x) -x^3; % max Vol = - min -Vol
A = 5;
b = 97;

x = fmincon(f, 1, A, b)