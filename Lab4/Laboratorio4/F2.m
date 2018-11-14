function wdot = F2(t,w)
% w tiene dim 6
alpha = 0.12;
b = 2.85;
M=300.83;
m=0.00083;
c = @(t) 0.0288 + (0.0613 - 0.0288)*((t-5)^2)/25;
d = @(t) 0.287 + (0.767 - 0.287)*((t-5)^2)/25;
x = w(1);
y = w(2);
z = w(3);
px = w(4);
py = w(5);
pz = w(6);
xdot = uopt(w(4:6)) + vopt(w(4:6)) - alpha*x;
ydot = -b*y*z + sqrt((M-uopt(w(4:6)))*(m+uopt(w(4:6))));
zdot = z*(c(t)*y -d(t)) - vopt(w(4:6));
pxdot = alpha*px;
pydot = py*b*z - pz*z*c(t);
pzdot = py*b*y - pz*y*c(t);
wdot = [xdot;ydot;zdot;pxdot;pydot;pzdot];
end

