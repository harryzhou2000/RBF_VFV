
clear;
%%

nx = 200;
ny = 100;
xdim = 2;
ydim = 1;
vB = 0;
vB2 = 1.2;
rhoB = 1;
rhoB2 = 1;
eB = 0.6;
eB2 = 2;
gm = 1.4;
dt = 0.01;

pB = (gm-1)*(eB - 0.5*vB^2);
aB = sqrt(gm*pB/rhoB);
MaB = vB/aB
pB2 = (gm-1)*(eB2 - 0.5*vB2^2);
aB2 = sqrt(gm*pB2/rhoB2);
MaB2 = vB2/aB2
%%
xs = linspace(0,xdim,nx+1);
ys = linspace(0,ydim,ny+1);
xc = 0.5*(xs(1:end-1) + xs(2:end));
yc = 0.5*(ys(1:end-1) + ys(2:end));
[xm,ym] = meshgrid(xc,yc);
sz = size(xm);
dx = xs(2) - xs(1);
dy = ys(2) - ys(1);

u0 = vB*ones(size(xm));
v0 = zeros(size(xm));
rho0 = ones(size(xm))*rhoB;
e0   = ones(size(xm))*eB;
cent = ((xm-0.5).^2 + (ym-0.5).^2) < 0.2^2;
rho0(cent) = 5*rhoB;
e0(cent) = 15*eB;
solid = ((((xm-0.).^2/0.4^2 + (ym-0.7).^2/0.2^2) < 1) |...
    (((xm-0.).^2/0.4^2 + (ym-0.3).^2/0.2^2) < 1) )&...
    (abs(ym-0.5)-0.02)>sqrt(xm+0.1)/5;
solid = false(size(xm));
solid = reshape(solid,1,[]);
US = ones(1,sum(solid)).*[1;0;0;1];


[js,is] = meshgrid(1:nx,1:ny);
ilower = [is(end,:);  is(1:end-1,:)];
iupper = [is(2:end,:);is(1,:)      ];
jleft  = [js(:,end),  js(:,1:end-1)];
jright = [js(:,2:end),js(:,1)      ];

ilower = reshape(ilower,1,[]);
iupper = reshape(iupper,1,[]);
jleft  = reshape(jleft ,1,[]);
jright = reshape(jright,1,[]);
is = reshape(is,1,[]);
js = reshape(js,1,[]);
Ilower = sub2ind(sz,ilower,js);
Iupper = sub2ind(sz,iupper,js);
Ileft  = sub2ind(sz,is,jleft);
Iright = sub2ind(sz,is,jright);


U = [reshape(rho0,1,[]);reshape(u0,1,[]);reshape(v0,1,[]);reshape(e0,1,[])];
U(2,:) = U(2,:) .* U(1,:);
U(3,:) = U(3,:) .* U(1,:);
Uprev = U;
%%
B = false(size(xm));
B(:,1) = true; B(:,end) = true; B(1,:) = true; B(end,:) = true;
B2 = B & xm<0.5 & ym>0.45 & ym<0.55;
B = reshape(B,1,[]);
B2 = reshape(B2,1,[]);
tvd = 1;
see = 10;
cfl = 0.5;
UB = [rhoB;rhoB*vB;0;eB] .* ones(1,sum(B));
UB2 = [rhoB2;rhoB2*vB2;0;eB2] .* ones(1,sum(B2));
for iter = 1:2500
    
    
    %
    %     deltaL = U - U(:,Ileft);
    %     sigL = deltaL/dx.*TVDLimiter(deltaL(:,Iright),deltaL) *tvd;
    %     sigR = deltaL(:,Iright)/dx.*TVDLimiter(deltaL,deltaL(:,Iright)) *tvd;
    %     deltaLo = U - U(:,Ilower);
    %     sigLo = deltaLo/dy.*TVDLimiter(deltaLo(:,Iupper),deltaLo) *tvd;
    %     sigUp = deltaLo(:,Iupper)/dy.*TVDLimiter(deltaLo,deltaLo(:,Iupper)) *tvd;
    %
    %
    %     [leftF,vm1] =  RoeSolver2(U(:,Ileft) + sigR(:,Ileft)*dx/2, U - sigL*dx/2,gm);
    %     [downF,vm2] = RoeSolver2c(U(:,Ilower) + sigLo(:,Ilower)*dy/2,U - sigUp*dy/2,gm);
    %     dta = min(dx/vm1,dy/vm2) * cfl;
    %     inc = (leftF(:,Iright) - leftF) / dx + (downF(:,Iupper) - downF) /dy;
    %     U = U - inc*min(dt,dta);
    %     U(1,U(1,:)<0) = 0;
    tvdc = tvd;
    aim = true;
    while aim
        [dU1,dta1,aim1] = increase(U,B,B2,solid,UB,UB2,US,Ileft,Iright,Iupper,Ilower,tvdc,dx,dy,gm,cfl);
        dtc = min(dt,dta1);
        U1 = U + dtc * dU1;
        [dU2,dta2,aim2] = increase(U1,B,B2,solid,UB,UB2,US,Ileft,Iright,Iupper,Ilower,tvdc,dx,dy,gm,cfl);
        aim = aim1||aim2;
        tvdc = tvdc*0.7;
    end
    Uprev = U;
    U = 0.5*U + 0.5*U1 + 0.5*dtc*dU2;
    
    
    
    
    %     U =U + dU1*dtc;
    
    
    if(mod(iter,see)==0)
        U(:,B) = UB;
        U(:,B2) = UB2;
        U(1:4,solid) = US(1:4,:);
        rhoC = U(1,:);
%         rhoC = U(2,:)./U(1,:);
        rhoC = reshape(rhoC,sz);
        surf(xm,ym,rhoC,'LineStyle','none');
        view([0,0,1]);
        axis equal;
        fprintf('iter = %d\n',iter);
        drawnow;
    end
end

%%
rhoC = U(1,:);
u = U(2,:)./U(1,:);
v = U(3,:)./U(1,:);
rhoC = reshape(rhoC,sz);
u = reshape(u,sz);
v = reshape(v,sz);
% surf(xm,ym,rhoC,'LineStyle','none');
% view([0,0,1]);
% quiver(xm,ym,u,v);
contour(xm,ym,rhoC,32);
% streamslice(xm,ym,u,v,11);
% contour(xm,ym,curl(u,v),32);
axis equal;
fprintf('iter = %d\n',iter);
drawnow;


%%
function [inc,dta,aim] = increase(U,B,B2,solid,UB,UB2,US,Ileft,Iright,Iupper,Ilower,tvd,dx,dy,gm,cfl)
U(:,B) = UB;
U(:,B2) = UB2;
U(1:4,solid) = US(1:4,:);

deltaL = U - U(:,Ileft);
sigL = deltaL/dx.*TVDLimiter(deltaL(:,Iright),deltaL) *tvd;
sigR = deltaL(:,Iright)/dx.*TVDLimiter(deltaL,deltaL(:,Iright)) *tvd;
deltaLo = U - U(:,Ilower);
sigLo = deltaLo/dy.*TVDLimiter(deltaLo(:,Iupper),deltaLo) *tvd;
sigUp = deltaLo(:,Iupper)/dy.*TVDLimiter(deltaLo,deltaLo(:,Iupper)) *tvd;


[leftF,vm1, aim1] =  RoeSolver2(U(:,Ileft) + sigR(:,Ileft)*dx/2, U - sigL*dx/2,gm);
[downF,vm2, aim2] = RoeSolver2c(U(:,Ilower) + sigLo(:,Ilower)*dy/2,U - sigUp*dy/2,gm);
aim = aim1 || aim2;

dta = min(dx/vm1,dy/vm2) * cfl;
inc = (leftF(:,Iright) - leftF) / dx + (downF(:,Iupper) - downF) /dy;
inc = -inc;

end