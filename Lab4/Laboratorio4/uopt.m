function valor = uopt(p)
% p de dim 3
M=300.83;
m=0.00083;
px = p(1);
py = p(2);
aux = ((M-m)/2 - (M+m)*px/(2*sqrt(px^2 + py^2)));
valor = aux*(py <0)*(aux > 0) + M*(py>0)*((px*sqrt(M)) < (py*sqrt(m)));
end

