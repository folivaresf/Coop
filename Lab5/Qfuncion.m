function Q = Qfuncion(f1,f2)
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

dim = (nx+1)*(nv+1);

k = 1:dim;

Q=sparse(dim,dim);
for i=1:nx+1
    for j=1:nv+1
        k = (nx+1)*(j-1) + i;
        Q(k,k) = -abs(f1(i,j)/hx)-abs(f2(i,j)/hv);
        if k-1 >0
            Q(k,k-1) = -f1(i,j)*(f1(i,j)<0)/hx;
        end
        if k+1<=dim
            Q(k,k+1) = f1(i,j)*(f1(i,j)>=0)/hx;
        end
        if k-(nx+1)>0
            Q(k,k-(nx+1)) = -f2(i,j)*(f2(i,j)<0)/hv;
        end
        if k+(nx+1)<=dim
            Q(k,k+(nx+1)) = f2(i,j)*(f2(i,j)>=0)/hv;
        end
    end
end
end

