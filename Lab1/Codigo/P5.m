% P5
clear all

f = @(x) [x(1)^3 + x(2)^3; x(1) - x(2); x(1) + x(2) + 7];
lb = [-1,-1];
ub = [0,0];

[x,fval] = fminimax(@myfun,[0.1,0.1], [], [], [], [], lb, ub);

% Por si acaso
function f = myfun(x)
    f(1) =  x(1)^3 + x(2)^3;
    f(2) = x(1) - x(2);
    f(3) = x(1) + x(2) + 7;
end
