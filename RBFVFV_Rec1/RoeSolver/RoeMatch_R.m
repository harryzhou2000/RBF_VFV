clear;
close all;
UL = [1,0.75,1]';
UR = [0.125,0,0.1]';

% u2U(U2u([1;2;3],1.4),1.4)
N = 1000;
xs = linspace(0.00,1,N+1); 
hx = xs(2) - xs(1);
As = abs(xs) * pi * 2; Vs = abs(diff(xs.^2 * pi));
% As = ones(size(xs)); Vs = hx;
xc = (xs(2:end) + xs(1:end-1))/2;

u0 = ones(3,N);
g = 1.4;
% u0(:,xc<=0) = u0(:,xc<=0).* U2u(UL,g);
% u0(:,~(xc<=0)) =u0(:,~(xc<=0)).* U2u(UR,g);
u0(:,:) = [1.4;0;1/(g-1)] .* (0.0*rand(size(xc)) + 1);
u0(:,abs(xc)<=0.01) = u0(:,abs(xc)<=0.01) .* [1;1;10];

%%
dtmax = 0.1/N;
% dtmax = 1e-4;
u = u0;
for iter = 1:round(N)*40
    [dudt1,vm] = TVD1(u,g,hx,As,Vs);
    dt = min(hx/max(vm) * 0.5,dtmax); 
    um1 = u + dt * dudt1;
    um2 = 0.75 * u + 0.25 * (um1 + dt * TVD1(um1,g,hx,As,Vs));
    u = 1/3*u + 2/3*(um2 + dt * TVD1(um2,g,hx,As,Vs));
    
    plot(u2U(u,g)')
    fprintf('iter = %d dt = %e cfl = %e\n',iter ,dt, dt/(hx/max(vm)));
    drawnow;
end
%%
plot(u2U(u,g)')

function u = U2u(U,g)
u = U;
u(2,:) = u(1,:) * u(2,:);
u(3,:) = u(3,:)/(g-1) + 0.5 * u(2,:).^2./u(1,:);
end


function U = u2U(u,g)
U = u;
U(2,:) =  U(2,:) ./ U(1,:);
U(3,:) = (U(3,:) - 0.5 * U(2,:).^2 .* U(1,:))*(g-1);
end


function [dFdt,VM] = TVD1(u,g,h,As,Vs)
% u(:,1) = u(:,2);
delta = [zeros(3,1),u(:,2:end) - u(:,1:end-1),zeros(3,1)];
sigR = delta(:,2:end)/h.*TVDLimiter(delta(:,1:end-1),delta(:,2:end)) * 0.5 * 1;
sigR(:,1) = 0;
[cc,VM]=RoeSolver(...
    [(u(:,2) - sigR(:,2)*h/2).*[1;-1;1],u(:,1:end-1)+sigR(:,1:end-1)*h/2], [u(:,1:end) - sigR(:,1:end)*h/2] ...
    ,g);
dFdt = -diff([cc,zeros(3,1)].*As, ...
    1,2) ./ Vs;
U = u2U(u,g);
dFdt(2,:) = dFdt (2,:) + h * 2 * pi * U(3,:)./Vs;
% dFdt(:,1) = 0;
dFdt(:,end) = 0;
% dFdt(:,1) = dFdt(:,2);
end