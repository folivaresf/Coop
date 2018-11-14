function valor = vopt(p)
% p de dim 3
V=5;
px = p(1);
pz = p(3);
valor = V*(px < pz);
end

