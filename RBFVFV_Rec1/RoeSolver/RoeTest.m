clear;
k = 5;

xs = linspace(0,1,1001);
h = xs(2)-xs(1);
xc = (xs(1:end-1) + xs(2:end))/2;
us = ones(3,numel(xc));
us(1,xc<0.5) = 0.445;
us(1,xc>=0.5) = 0.5;
us(2,xc<0.5) = 0.698;
us(2,xc>=0.5) = 0;
us(2,:) = us(2,:).*us(1,:);
us(3,xc<0.5) = 3.528;
us(3,xc>=0.5) = 0.571;
us(3,:) = us(3,:)/(1.4-1) + us(2,:).^2 ./ us(1,:) *.5;
dt = 0.0001;
right = (1:numel(xc))+1;
right(end) = 1;
left = (1:numel(xc))-1;
left(1) = numel(xc);
figure(2);
ax = gca;
for iter = 1:0.1/dt
    delta = us - us(:,left);
    sigL = delta/h.*TVDLimiter(delta(right),delta) *0;
    sigR = delta(right)/h.*TVDLimiter(delta,delta(right)) *0;
    
    fs = RoeSolver(us(:,left)+sigR(:,left)*h/2,us-sigL*h/2,1.4);
    incs = -[diff(fs,1,2),fs(:,1)-fs(:,end)]/h;
    incs(:,1) = 0;
    incs(:,end) = 0;
    us = us + dt*incs;
    plot(ax,xc,us');
    drawnow;
     fprintf('Iter = %d\n',iter);
end








