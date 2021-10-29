function L = WBAP_L2(a,n,p)
if nargin == 1
    n = 1;
    p = 4;
elseif nargin == 2
    p = 4;
end
eps = -1e-10;
t = a(:,2:end)./(a(:,1)+eps) + eps;

L = a(:,1).*(n + sum(1./t.^(p-1),2))./(n + sum(1./t.^p,2));