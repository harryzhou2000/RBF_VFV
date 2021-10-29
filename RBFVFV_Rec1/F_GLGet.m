function [GLNodes,GLWeights] = F_GLGet(ngl)
nlev = ngl+1; %nlev-1 is rank of formula

p = sym('p',[1,nlev]);
syms x;
p0 = 1;
p(1) = x;
p(2) = ((2*1+1)*x*p(1) - 1 * p0) / (1+1);

for n = 2:(nlev-1)
    p(n+1) = (p(n)*x*(2*n+1) - n*p(n-1))/(n+1);
    
end
for n=2:nlev
    p(n) = expand(p(n)); 
end

fp7 = symfun(p(nlev),x);
mfp7 = matlabFunction(fp7);

r = 1;
f = fp7;
f = diff(f);
v = zeros(1,nlev);
for k = 1:nlev
    r = r*k;
    v(k) = f(0)/r;
    f = diff(f);
end
v = double(v);
v = [v(end:-1:1),mfp7(0)];
nodes = sort(roots(v));
f = diff(fp7);
mf = matlabFunction(f);
weights = 2./(1-nodes.^2)./(mf(nodes).^2);

GLWeights = weights;
GLNodes = nodes;