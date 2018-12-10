xmin = -1;
vmin = -1;

xmax = 1;
vmax = 1;

hx = 0.01;
hv = 0.01; 

nx = (xmax - xmin)/hx;
nv = (vmax - vmin)/hv;

x = xmin:hx:xmax;
v = vmin:hv:vmax;
[X,V] = meshgrid(x,v);

U = ceil(rand(length(v),length(x))-0.5);

[f1,f2,l]= funcion([X V],U);
figure(1)
quiver(X,V,f1,f2)
figure(2)
mesh(X,V,f1)
figure(3)
mesh(X,V,f2)

Q = Qfuncion(f1,f2);