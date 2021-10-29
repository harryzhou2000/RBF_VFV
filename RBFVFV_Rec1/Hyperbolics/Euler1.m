syms u1 u2 u3 p gm lam a rho u
p = (u3-0.5*u2^2/u1)*(gm-1);
F = [u2;u2^2/u1 + p; u2/u1*(u3+p)];
U = [u1;u2;u3];
A = simplify(expand(jacobian(F,U)));
aext = simplify(expand(sqrt((u3/u1-0.5*(u2/u1)^2)*(gm^2-gm))));
pretty(A)
spec = simplify(expand(det(A-eye(3)*lam)));
pretty(spec)
specFs = factor(spec,lam);
pretty(specFs)

[V,D] = eig(A);
simplify(expand(A*V - V*D))

V(:,1) = V(:,1)/V(1,1);
V(:,2) = V(:,2)/V(1,2);
V(:,3) = V(:,3)/V(1,3);
V = subs(simplify(expand(V)),[aext,u1,u2],[a,rho,rho*u]);