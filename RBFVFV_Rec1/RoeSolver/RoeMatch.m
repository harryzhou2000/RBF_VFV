clear;
% close all;
% UL = [1,0.75,1]';
% UR = [0.125,0,0.1]';
UL = [1,0,10]';
UR = [1,0,1]';

% u2U(U2u([1;2;3],1.4),1.4)
N = 1000;
xs = linspace(-0.5,0.5,N+1); 
hx = xs(2) - xs(1);
As = ones(size(xs)); Vs = hx;
xc = (xs(2:end) + xs(1:end-1))/2;

u0 = ones(3,N);
g = 1.4;
u0(:,xc<=0) = u0(:,xc<=0).* U2u(UL,g);
u0(:,~(xc<=0)) =u0(:,~(xc<=0)).* U2u(UR,g);
% u0(:,:) = [1.4;0;1/(g-1)] .* (xc*0 + 1);
% u0(:,abs(xc)<=0.1) = u0(:,abs(xc)<=0.1) .* [1;1;1000];

%%
dtmax = 0.1/N;
% dtmax = 1e-4;
u = u0;
t = 0;
for iter = 1:round(N)*4
    [dudt1,vm] = TVD1(u,g,hx,As,Vs);
    dt = min(hx/max(vm) * 0.5,dtmax); 
    um1 = u + dt * dudt1;
    um2 = 0.75 * u + 0.25 * (um1 + dt * TVD1(um1,g,hx,As,Vs));
    u = 1/3*u + 2/3*(um2 + dt * TVD1(um2,g,hx,As,Vs));
    t = t+dt;
    if(t>=4.83e-2)
        break;
    end
    plot(u2U(u,g)')
    fprintf('iter = %d t = %e dt = %e cfl = %e\n',iter,t ,dt, dt/(hx/max(vm)));
    drawnow;
end
%% 
U = u2U(u,g);
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
delta = [zeros(3,1),u(:,2:end) - u(:,1:end-1),zeros(3,1)];
sigR = delta(:,2:end)/h;
sigR = sigR.*TVDLimiter(delta(:,1:end-1),delta(:,2:end)) * 0.5 * 1;
Umax = max(u(:,[1,1:end-1]),u(:,[2:end,end]));
Umin = min(u(:,[1,1:end-1]),u(:,[2:end,end]));
Umax = max(Umax,u);Umin = min(Umin,u);
sigMax = min((Umax - u)/h*2,(u-Umin)/2*2);
sigCompress = min(sigMax./(abs(sigR) + 1e-10) , 1);
% sigR = sigR .* sigCompress;
UL = u(:,1:end-1)+sigR(:,1:end-1)*h/2;
UR = u(:,2:end) - sigR(:,2:end)*h/2;

[cc,VM]=RoeSolver(...
    UL, UR ...
    ,g);
dFdt = -diff([zeros(3,1),cc,zeros(3,1)].*As, ...
    1,2) ./ Vs;
dFdt(:,1) = 0;
dFdt(:,end) = 0;
end