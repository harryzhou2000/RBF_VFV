function x = F_CGSolve(A,b)
down = 1e-3;
downA = 1e-10;
itmax = 20;
x = zeros(size(b));
r = b - A*x;
p = r;
rsqr0 = r'*r;
rsqrt = rsqr0*down^2;
rsqr = rsqr0;
if rsqr<=downA^2
   fprintf('CG Iter = %d, res2 = %.4e\n',0,sqrt(rsqr));
   return;
end
for iter = 1:itmax
   alpha =  r'*r/(p'*(A*p));
   x = x + alpha * p;
   r = r - alpha*A*p;
   rsqrP = rsqr;
   rsqr = r'*r;
   beta = rsqr/rsqrP;
   p = p*beta + r;
   if(rsqr<=rsqrt || rsqr<=downA^2)
       fprintf('CG Iter = %d, res2 = %.4e\n',iter,sqrt(rsqr));
       break;
   end
    
end
if(iter == itmax)
    fprintf('CG Iter = %d, res2 = %.4e\n',iter,sqrt(rsqr));
end