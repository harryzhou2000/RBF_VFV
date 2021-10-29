nx = 100;
ny = 100;
eps_lim = 100; %Òª×ã¹»´ó 1+
N_lim = 2;

dx = 1/nx;
dy = 1/ny;

is = 1:ny;
js = 1:nx;
[jm,im] = meshgrid(js,is);
xc = linspace(0.5/nx,1-0.5/nx,nx);
yc = linspace(0.5/ny,1-0.5/ny,ny);
[xm,ym] = meshgrid(xc,yc);

jmLeft = circshift(jm,1,2);
jmRight = circshift(jm,-1,2);

imLo = circshift(im,1,1);
imUp = circshift(im,-1,1);

figure(99123);
u0 = double(((xm-0.5).^2 + (ym-0.5+0.14).^2 < 0.15^2) |...
    ((xm-0.5).^2 + (ym-0.5-0.14).^2 < 0.15^2));
surf(xm,ym,u0,'LineStyle','none');
view([0,90]);
        

a1 = 1;
a2 = 1;
% dt = 0.1 * sqrt(dx^2+dy^2)/sqrt(a1^2+a2^2);
dt = 0.1 * dx / a1;

iter = nx/a1*1*10;
see = 50;

%% Rotational 
figure(1);
u = u0;
for i = 1:iter
    dudx = (u(is,circshift(js,-1)) - u(is,circshift(js,1)))/(2*dx);
    dudy = (u(circshift(is,-1),js) - u(circshift(is,1),js))/(2*dy);
%     surf(xm,ym,dudx,'LineStyle','none');
    
    dudxLo =    reshape(dudx(circshift(is,1),js),1,[]);
    dudxUp =    reshape(dudx(circshift(is,-1),js),1,[]);
    dudxLeft =  reshape(dudx(is,circshift(js,1)),1,[]);
    dudxRight = reshape(dudx(is,circshift(js,-1)),1,[]);
    dudyLo =    reshape(dudy(circshift(is,1),js),1,[]);
    dudyUp =    reshape(dudy(circshift(is,-1),js),1,[]);
    dudyLeft =  reshape(dudy(is,circshift(js,1)),1,[]);
    dudyRight = reshape(dudy(is,circshift(js,-1)),1,[]);
    
    dudx = reshape(dudx,1,[]);
    dudy = reshape(dudy,1,[]);
    
    dudxi = limit_t([dudx;dudy], ...
        reshape([dudxLo;dudyLo;dudxLeft;dudyLeft;dudxUp;dudyUp;dudxRight;dudyRight],...
        2,4,[]),eps_lim,N_lim);
    dudx = reshape(dudxi(1,:),ny,nx);
    dudy = reshape(dudxi(2,:),ny,nx);

%     dudx = limit_t(dudx,reshape([dudxLo;dudxLeft;dudxUp;dudxRight],1,4,[]),...
%         eps_lim, N_lim);
%     dudy = limit_t(dudy,reshape([dudyLo;dudyLeft;dudyUp;dudyRight],1,4,[]),...
%         eps_lim, N_lim);
%     dudx = reshape(dudx,ny,nx);
%     dudy = reshape(dudy,ny,nx);
    
    
    
    uright = 0.5*dx *dudx + u;
    uup = 0.5*dy *dudy + u;
    u = u + (uup(circshift(is,1),js)*dx * a1 + uright(is,circshift(js,1))*dy * a2...
        -( uup(is,js)*dx * a1 + uright(is,js)*dy * a2))/(dx*dy)*dt;
    
    if(mod(i,see)==0)
        surf(xm,ym,u,'LineStyle','none');
        view([0,90]);
        drawnow;
    end
    
end

u_cc = u;

%% original
figure(2);
u = u0;
for i = 1:iter
    dudx = (u(is,circshift(js,-1)) - u(is,circshift(js,1)))/(2*dx);
    dudy = (u(circshift(is,-1),js) - u(circshift(is,1),js))/(2*dy);
%     surf(xm,ym,dudx,'LineStyle','none');
    
    dudxLo =    reshape(dudx(circshift(is,1),js),1,[]);
    dudxUp =    reshape(dudx(circshift(is,-1),js),1,[]);
    dudxLeft =  reshape(dudx(is,circshift(js,1)),1,[]);
    dudxRight = reshape(dudx(is,circshift(js,-1)),1,[]);
    dudyLo =    reshape(dudy(circshift(is,1),js),1,[]);
    dudyUp =    reshape(dudy(circshift(is,-1),js),1,[]);
    dudyLeft =  reshape(dudy(is,circshift(js,1)),1,[]);
    dudyRight = reshape(dudy(is,circshift(js,-1)),1,[]);
    
    dudx = reshape(dudx,1,[]);
    dudy = reshape(dudy,1,[]);
    

    dudx = limit_t(dudx,reshape([dudxLo;dudxLeft;dudxUp;dudxRight],1,4,[]),...
        eps_lim, N_lim);
    dudy = limit_t(dudy,reshape([dudyLo;dudyLeft;dudyUp;dudyRight],1,4,[]),...
        eps_lim, N_lim);
    dudx = reshape(dudx,ny,nx);
    dudy = reshape(dudy,ny,nx);
    
    
    
    uright = 0.5 * dx * dudx + u;
    uup = 0.5 * dy * dudy + u;
    u = u + (uup(circshift(is,1),js)*dx * a1 + uright(is,circshift(js,1))*dy * a2...
        -( uup(is,js)*dx * a1 + uright(is,js)*dy * a2))/(dx*dy)*dt;
    
    if(mod(i,see)==0)
        surf(xm,ym,u,'LineStyle','none');
        view([0,90]);
        drawnow;
    end
    
end


%%
figure(3);
u_cr = u;
cut = round(nx*0.5);
plot(ym(:,cut),u_cc(:,cut),ym(:,cut),u_cr(:,cut));