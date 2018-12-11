xmin = -1;
vmin = -1;

xmax = 1;
vmax = 1;

% hx = 0.01;
% hv = 0.01; 

hx = 1;
hv = 1;

nx = (xmax - xmin)/hx;
nv = (vmax - vmin)/hv;

x = xmin:hx:xmax;
v = vmin:hv:vmax;
[X,V] = meshgrid(x,v);

U = ceil(rand(length(v)*length(x), 1)-0.5);

[f1,f2,l]= funcion([X V], U);
figure(1)
ind = 1:10:length(f1);
quiver(X(ind, ind),V(ind, ind),f1(ind, ind),f2(ind, ind))
figure(2)
mesh(X,V,f1)
figure(3)
mesh(X,V,f2)

Q = Qfuncion(f1,f2);

V0 = ones((nx+1)*(nv+1), 1);
u0 = ones((nx+1)*(nv+1), 1);

% [V u] = howard(V0, u0, X, V);
