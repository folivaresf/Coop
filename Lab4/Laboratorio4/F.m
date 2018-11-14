function zdot = F(t,z)
alpha = 1;
AA = [-alpha 0; -2 alpha];
p = z(2);
ut =(p<(1/sqrt(3)))*(1-(2*p)/(sqrt(p^2+1))); %control optimo
zdot = AA*z + [ut;0]; % dinamica de (x,p)
end
