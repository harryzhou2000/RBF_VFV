
clear;
%%

nx = 64*3; lx = 3;
ny = 64; ly = 1;
tmax = 4;
vB = 0.0001;
vB2 = 0.001;
rhoB = 1.4;
rhoB2 = 2;
pB = 1;
% eB = 7;
eB2 = 10;
gm = 1.4;
dt = 0.01;

eB = pB/(gm-1) + 0.5*vB^2*rhoB;
aB = sqrt(gm*pB/rhoB);
MaB = vB/aB
pB2 = (gm-1)*(eB2 - 0.5*vB2^2);
aB2 = sqrt(gm*pB2/rhoB2);
MaB2 = vB2/aB2

%%
xs = linspace(0,lx,nx+1);
ys = linspace(0,ly,ny+1);
xc = 0.5*(xs(1:end-1) + xs(2:end));
yc = 0.5*(ys(1:end-1) + ys(2:end));
[xm,ym] = meshgrid(xc,yc);
sz = size(xm);
dx = xs(2) - xs(1);
dy = ys(2) - ys(1);

u0 = vB*(xm-0.5).^0*1;

v0 = zeros(size(xm));
rho0 = ones(size(xm))*rhoB;
e0   = ones(size(xm))*eB;
cent = ((xm-0.5).^2 + (ym-0.5).^2) < 0.2^2;
% rho0(cent) = 4*rhoB;
% e0(cent) = 10*eB;
% u0(cent) = 0;
solid = (((xm-0.4).^2/0.1^2 + (ym-0.6).^2/0.01^2) < -1) |...
    (((xm-0.4).^2/0.1^2 + (ym-0.4).^2/0.01^2) < -1) |...
    ym<0.05  | ym>0.95 | (ym<(0.05+0.2*0.9)&xm>1&xm<3);


xb = diff(solid,1,2);
yb = diff(solid,1,1);
solidRight = [xb>0,false(sz(1),1)];
solidLeft  = [false(sz(1),1),xb<0];
solidUpper = [yb>0;false(1,sz(2))];
solidLower = [false(1,sz(2));yb<0];
imagesc(xc,yc,solidLower | solidUpper | solidLeft);
set(gca,'YDir','normal');

solid = reshape(solid,1,[]);
solidRight = reshape(solidRight,1,[]);
solidLeft = reshape(solidLeft,1,[]);
solidUpper = reshape(solidUpper,1,[]);
solidLower = reshape(solidLower,1,[]);
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
t = 0;

rec = true;

recorder = VideoWriter('RoeTest2rec.mp4','MPEG-4');
recorder.FrameRate = 60;
recorder.open();


%%

B = false(size(xm));
B(:,1) = true; B(:,end) = true; B(1,:) = true; B(end,:) = true;

xb = diff(B,1,2);
yb = diff(B,1,1);
BsolidRight = [xb>0,false(sz(1),1)];
BsolidLeft  = [false(sz(1),1),xb<0];
BsolidUpper = [yb>0;false(1,sz(2))];
BsolidLower = [false(1,sz(2));yb<0];
BsolidRight = reshape(BsolidRight,1,[]);
BsolidLeft = reshape(BsolidLeft,1,[]);
BsolidUpper = reshape(BsolidUpper,1,[]);
BsolidLower = reshape(BsolidLower,1,[]);

Bx = xm(B);
By = ym(B);
B = reshape(B,1,[]);
Bx = reshape(Bx,1,[]);
By = reshape(By,1,[]);
tvd = 1;
see = 20;
cfl = 0.5;
sav = 20;
pvt = 0.001;
UB = [rhoB;rhoB*vB;0;eB] .* ones(1,sum(B));
Breg2 = Bx<0.5 & By > 0.45 & By<0.55;
Breg3 = Bx>0.5 & By > 0.45 & By<0.55;
% UB(:,Breg2) = [rhoB2;rhoB2*vB2;0;eB2]*ones(1,sum(Breg2));
% UB(:,Breg3) = [rhoB2;-rhoB2*vB2;0;eB2]*ones(1,sum(Breg3));



% (U,B,B2,solid,UB,UB2,US,Ileft,Iright,Iupper,Ilower,tvdc,dx,dy,gm,cflc)
U = gpuArray(U);
UB = gpuArray(UB);
B = gpuArray(B);
solid = gpuArray(solid);
US = gpuArray(US);
Ileft = gpuArray(Ileft);
Iright = gpuArray(Iright);
Iupper = gpuArray(Iupper);
Ilower = gpuArray(Ilower);
BsolidLeft = gpuArray(BsolidLeft);
BsolidLower = gpuArray(BsolidLower);
BsolidRight = gpuArray(BsolidRight);
BsolidUpper = gpuArray(BsolidUpper);
solidLeft = gpuArray(solidLeft);
solidLower = gpuArray(solidLower);
solidRight = gpuArray(solidRight);
solidUpper = gpuArray(solidUpper);

% dx = gpuArray(dx);
% dy = gpuArray(dy);
% tvd = gpuArray(tvd);
% dt = gpuArray(tvd);
sv1 = U;
sv2 = U;
tic;

% load U1
phi = gpuArray(reshape((xm*vB),1,[]));
%% Initializer
for iter = 1:0  
    phiN = (dy^2*(phi(Ileft) + phi(Iright)) + dx^2*(phi(Iupper) + phi(Ilower)))/(dx^2+dy^2)*0.5;
    phiN(solid) = 0;
    phiN(BsolidLeft(Iright)) = phiN(BsolidLeft) - vB*dx;
    phiN(BsolidRight(Ileft)) = phiN(BsolidRight) + vB*dx;
    phiN(BsolidUpper(Ilower)) = phiN(BsolidUpper);
    phiN(BsolidLower(Iupper)) = phiN(BsolidLower);
    
    phiN(solidLeft(Iright)) = phiN(solidLeft);
    phiN(solidRight(Ileft)) = phiN(solidRight);
    phiN(solidUpper(Ilower)) = phiN(solidUpper);
    phiN(solidLower(Iupper)) = phiN(solidLower);
    
    phiN = phiN - mean(phiN,'all');
    phiO = phi;
    phi = phiN;
    
    if mod(iter,1) ==0
        res = max(phiO-phiN,[],'all');
        fprintf('iter = %d, res = %e\n',iter,res);
        U(2,:) = (phi(Iright)-phi(Ileft))/2.*U(1,:)/dx;
        U(3,:) = (phi(Iupper)-phi(Ilower))/2.*U(1,:)/dy;
        U(2,solid) = 0;
        U(3,solid) = 0;
        U(2,B) = 0;
        U(3,B) = 0;
        subplot(1,2,1)
        imagesc(xc,yc,reshape(U(2,:)./U(1,:),sz));
        colorbar;
        set(gca,'YDir','normal');
        drawnow;
        subplot(1,2,2)
        imagesc(xc,yc,reshape(U(3,:)./U(1,:),sz));
        colorbar;
        set(gca,'YDir','normal');
        drawnow;
        
    end
    
end


%% RUN

B = false(size(xm));
B(:,1) = true; B(:,end) = true; B(1,:) = true; B(end,:) = true;

xb = diff(B,1,2);
yb = diff(B,1,1);
BsolidRight = [xb>0,false(sz(1),1)];
BsolidLeft  = [false(sz(1),1),xb<0];
BsolidUpper = [yb>0;false(1,sz(2))];
BsolidLower = [false(1,sz(2));yb<0];
BsolidRight = reshape(BsolidRight,1,[]);
BsolidLeft = reshape(BsolidLeft,1,[]);
BsolidUpper = reshape(BsolidUpper,1,[]);
BsolidLower = reshape(BsolidLower,1,[]);

Bx = xm(B);
By = ym(B);
B = reshape(B,1,[]);
Bx = reshape(Bx,1,[]);
By = reshape(By,1,[]);
tvd = 1;
tvddim = 0.5^4;
see = 20;
cfl = 0.5;
sav = 20;
pvt = 0.001;
UB = [rhoB;rhoB*vB;0;eB] .* ones(1,sum(B));
Breg2 = Bx<0.5 & By > 0.45 & By<0.55;
Breg3 = Bx>0.5 & By > 0.45 & By<0.55;
UB(:,Breg2) = [rhoB2; rhoB2*vB2;0;eB2]*ones(1,sum(Breg2));
UB(:,Breg3) = [rhoB2;-rhoB2*vB2;0;eB2]*ones(1,sum(Breg3));



% (U,B,B2,solid,UB,UB2,US,Ileft,Iright,Iupper,Ilower,tvdc,dx,dy,gm,cflc)
U = gpuArray(U);
UB = gpuArray(UB);
B = gpuArray(B);
solid = gpuArray(solid);
US = gpuArray(US);
Ileft = gpuArray(Ileft);
Iright = gpuArray(Iright);
Iupper = gpuArray(Iupper);
Ilower = gpuArray(Ilower);
BsolidLeft = gpuArray(BsolidLeft);
BsolidLower = gpuArray(BsolidLower);
BsolidRight = gpuArray(BsolidRight);
BsolidUpper = gpuArray(BsolidUpper);
solidLeft = gpuArray(solidLeft);
solidLower = gpuArray(solidLower);
solidRight = gpuArray(solidRight);
solidUpper = gpuArray(solidUpper);
for iter = 1:100000
    
    
    vB = 0.0001;
    %     vB = min((1+t),3 );
    eB = pB/(gm-1) + 0.5*vB^2*rhoB;
    UB = [rhoB;rhoB*vB;0;eB] .* ones(1,sum(B));
%     UB(:,Breg2) = [rhoB; rhoB*vB2;0;eB]*ones(1,sum(Breg2));
%     UB(:,Breg3) = [rhoB;-rhoB*vB2;0;eB]*ones(1,sum(Breg3));
    %
    %     deltaL = U - U(:,Ileft);
    %     sigL = deltaL/dx.*TVDLimiter(deltaL(:,Iright),deltaL) *tvd;
    %     sigR = deltaL(:,Iright)/dx.*TVDLimiter(deltaL,deltaL(:,I  right)) *tvd;
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
    %     tvdc = tvd*ones(1,size(U,2));
    tvdc = tvd;
    cflc = cfl;
    aim = true;
    ending = false;
    while aim
        [dU1,dta1,aim1] = increase(U,B,solid,UB,US,Ileft,Iright,Iupper,Ilower,...
            solidLeft,solidRight,solidUpper,solidLower,...
            BsolidLeft,BsolidRight,BsolidUpper,BsolidLower,...
            tvdc,dx,dy,gm,cflc,pvt);
        dtc = min(dt,min(dta1));
        if(dtc+t>tmax)
            dtc = tmax-t;
            ending = true;
        end
        U1 = U + dtc .* dU1;
        [dU2,dta2,aim2] = increase(U1,B,solid,UB,US,Ileft,Iright,Iupper,Ilower,...
            solidLeft,solidRight,solidUpper,solidLower,...
            BsolidLeft,BsolidRight,BsolidUpper,BsolidLower,...
            tvdc,dx,dy,gm,cflc,pvt);
        aim = aim1||aim2;
        fprintf('tvdc = %f\n',tvdc);
        tvdc = tvdc*tvddim;
        %         cflc = cflc*0.5;
        if(tvdc < 1e-8*tvd)
            U(:,B) = UB;
            
            U1 = (U(:,Iupper) + U(:,Ilower) +U(:,Ileft) +U(:,Iright))/4;
            U1(:,solidRight) = U(:,solidRight);
            U1(:,solidLeft) = U(:,solidLeft);
            U1(:,solidUpper) = U(:,solidUpper);
            U1(:,solidLower) = U(:,solidLower);
            U = U1*0.5 + U*0.5;
            fprintf('Damping!!\n');
        end
    end
    Uprev = U;
    U = 0.5*U + U1*(1-0.5) + 0.5*dtc.*dU2;
    t = t+min(dtc,[],'all');
    
    
    
    %     U =U + dU1*dtc;
    
    fprintf('iter = %d time = %f\n',iter,toc);
    if(mod(iter,see)==0)
        maxShow = 8 ;
        minShow = 0;
        maxShow2 = 10;
        minShow2 = 0;
        
        U(:,B) = UB;
        
        U(1:4,solid) = US(1:4,:);
        rhoC = U(1,:);
        u = U(2,:)./U(1,:);
        v = U(3,:)./U(1,:);
        p = (U(4,:) - 0.5*(u.^2+v.^2).*rhoC) * (gm-1);
        as = abs(gm*p./rhoC);
        %         rhoC = U(2,:)./U(1,:);
        rhoC = reshape(rhoC,sz);
        p = reshape(p,sz);
        %         as = reshape(as,sz);
        %%% Max min
%         rhoC(1) = maxShow;
%         rhoC(2) = minShow;
%         p(1) = maxShow2;
%         p(2) = minShow2;
        
        subplot(2,1,2);
        %         surf(xm,ym,max(min(p,maxShow2),minShow2),'LineStyle','none');
        imagesc(xc,yc,max(min(reshape(sqrt(u.^2./as),sz),maxShow2),minShow2));
        set(gca,'YDir','normal');
        view([0,0,1]);
        axis equal;
        colorbar;
        
        subplot(2,1,1);
        %         surf(xm,ym,max(min(rhoC,maxShow),minShow),'LineStyle','none');
        imagesc(xc,yc,max(min(rhoC,maxShow),minShow));
        set(gca,'YDir','normal');
        view([0,0,1]);
        axis equal;
        colorbar;
        
        fprintf('## iter = %d time = %f\n',iter,toc);
        tic;
        title(sprintf('t = %e',t));
        frame = getframe(gcf);
        %         frame = imresize(frame2im(frame),[430,544],'Method','bicubic');
        recorder.writeVideo(frame);
        drawnow;
    end
    
    if mod(iter,sav) == 0
        sv2 = sv1;
        sv1 = U;
    end
    
    if ending
       break; 
    end
end
toc;
%%
recorder.close();
rhoC = U(1,:);
u = U(2,:)./U(1,:);
v = U(3,:)./U(1,:);
u = max(min(u,vB*3),-vB*3);
v = max(min(v,vB*3),-vB*3);
rhoC = reshape(rhoC,sz);
u = reshape(u,sz);
v = reshape(v,sz);
% surf(xm,ym,rhoC,'LineStyle','none');
% view([0,0,1]);
quiver(xm,ym,u,v);
% contour(xm,ym,min(sqrt(u.^2+v.^2),3),32);
% contour(xm,ym,rhoC,32);
% streamslice(xm,ym,u,v,11);
% contour(xm,ym,curl(u,v),32);
axis equal;
fprintf('iter = %d\n',iter);
drawnow;


%%
function [inc,dta,aim] = increase(U,B,solid,UB,US,Ileft,Iright,Iupper,Ilower,...
    solidLeft,solidRight,solidUpper,solidLower,...
    BsolidLeft,BsolidRight,BsolidUpper,BsolidLower,...
    tvd,dx,dy,gm,cfl,pvt)


U(:,B) = UB;



U(1:4,solid) = US(1:4,:);


p = (U(4,:) - 0.5*(U(2,:).^2+U(3,:).^2)./U(1,:)) * (gm-1);
maxp = max(p,[],'all');
pviolate = (p<maxp*pvt);
U(4,pviolate) = 0.5*(U(2,pviolate).^2+U(3,pviolate).^2)./U(1,pviolate) + maxp*pvt*2/(gm-1);
maxRho = max(U(1,:),[],'all');
U(1,U(1,:)<maxRho*pvt) = maxRho*pvt;

deltaL = U - U(:,Ileft);
sigL = deltaL/dx.*TVDLimiter(deltaL(:,Iright),deltaL) *tvd;
sigR = deltaL(:,Iright)/dx.*TVDLimiter(deltaL,deltaL(:,Iright)) *tvd;
deltaLo = U - U(:,Ilower);
sigLo = deltaLo/dy.*TVDLimiter(deltaLo(:,Iupper),deltaLo) *tvd;
sigUp = deltaLo(:,Iupper)/dy.*TVDLimiter(deltaLo,deltaLo(:,Iupper)) *tvd;
sigL(B) = 0;
sigR(B)= 0;
sigLo(B) = 0;
sigUp(B) = 0;
% sigL(solidLeft) = 0;
% % sigR(solidLeft) = 0;
% % sigL(solidRight) = 0;
% sigR(solidRight) = 0;
% sigLo(solidLower) = 0;
% % sigUp(solidLower) = 0;
% % sigLo(solidUpper) = 0;
% sigUp(solidUpper) = 0;

sigL(solidLeft) = sigR(solidLeft);
sigR(solidRight) = sigL(solidRight);
sigLo(solidLower) = sigUp(solidLower);
sigUp(solidUpper) = sigLo(solidUpper);

sigL(BsolidLeft) = 0;
sigR(BsolidRight) = 0;
sigLo(BsolidLower) = 0;
sigUp(BsolidUpper) = 0;

uleft = U(:,Ileft) + sigR(:,Ileft)*dx/2;
uright = U - sigL*dx/2;
uleft(:,solidLeft) = uright(:,solidLeft);
uleft(2,solidLeft) = -uright(2,solidLeft);
uright(:,solidRight(Ileft)) = uleft(:,solidRight(Ileft));
uright(2,solidRight(Ileft)) = -uleft(2,solidRight(Ileft));

ulower = U(:,Ilower) + sigLo(:,Ilower)*dy/2;
uupper = U - sigUp*dy/2;
ulower(:,solidLower) = uupper(:,solidLower);
ulower(3,solidLower) = -uupper(3,solidLower);
uupper(:,solidUpper(Ilower)) = ulower(:,solidUpper(Ilower));
uupper(3,solidUpper(Ilower)) = -ulower(3,solidUpper(Ilower));

[leftF,vm1, aim1] =  RoeSolver2(uleft,uright ,gm);
[downF,vm2, aim2] = RoeSolver2c(ulower,uupper,gm);
aim = aim1 || aim2;

dta = min(dx./vm1,dy./vm2) * cfl;
inc = (leftF(:,Iright) - leftF) / dx + (downF(:,Iupper) - downF) /dy;
inc = -inc;

end
