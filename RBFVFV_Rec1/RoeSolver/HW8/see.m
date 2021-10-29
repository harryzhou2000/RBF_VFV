
dats = {...
    'R1_100_0.5_HY.mat',...
    'R1_500_0.5_HY.mat',...
    'R1_1000_0.5_HY.mat'};
Ns = [100 500 1000];

clf(1);
for it = 1:3
    load(dats{it});
    subplot(1,2,1);
    U = u2U(u,g);
    hold on;
    plot(xc,U(1,:),'.-','DisplayName',sprintf('N=%d',Ns(it)));
    subplot(1,2,2);
    hold on;
    plot(xc,U(2,:),'.-','DisplayName',sprintf('N=%d',Ns(it)));
end
subplot(1,2,1);
legend;xlim([-0.1,0.4])
grid on;grid minor;title('\rho, 1st Order, t = 0.1');
subplot(1,2,2);
legend;xlim([-0.1,0.4])
grid on;grid minor;title('U, 1st Order, t = 0.1');


%%
dats = {...
    'R2_100_0.5_HY.mat',...
    'R2_500_0.5_HY.mat',...
    'R2_1000_0.5_HY.mat'};
Ns = [100 500 1000];

clf(1);
for it = 1:3
    load(dats{it});
    subplot(1,2,1);
    U = u2U(u,g);
    hold on;
    plot(xc,U(1,:),'.-','DisplayName',sprintf('N=%d',Ns(it)));
    subplot(1,2,2);
    hold on;
    plot(xc,U(2,:),'.-','DisplayName',sprintf('N=%d',Ns(it)));
end
subplot(1,2,1);
legend;xlim([-0.1,0.4])
grid on;grid minor;title('\rho, 2nd Order, t = 0.1');
subplot(1,2,2);
legend;xlim([-0.1,0.4])
grid on;grid minor;title('U, 2nd Order, t = 0.1');
%%
clf;
load R1_500_0.5_HY.mat
hold on;
subplot(1,2,1);
U = u2U(u,g);
hold on;
plot(xc,U(1,:),'.-','DisplayName',sprintf('1st order'));
subplot(1,2,2);
hold on;
plot(xc,U(2,:),'.-','DisplayName',sprintf('1st order'));
load R2_500_0.5_HY.mat
hold on;
subplot(1,2,1);
U = u2U(u,g);
hold on;
plot(xc,U(1,:),'.-','DisplayName',sprintf('2nd order'));
subplot(1,2,2);
hold on;
plot(xc,U(2,:),'.-','DisplayName',sprintf('2nd order'));
XL = [0.1,0.16];
subplot(1,2,1);
legend;xlim(XL)
grid on;grid minor;title('\rho, N=500, t = 0.1');
subplot(1,2,2);
legend;xlim(XL)
grid on;grid minor;title('U, N=500, t = 0.1');
ylim([0,inf]);
%%%%%%%%%%%%%%%%%%
%%
clf;
load R2_1000_0.5_NHY.mat
hold on;
subplot(1,2,1);
U = u2U(u,g);
hold on;
plot(xc,U(1,:),'.-','DisplayName',sprintf('No entropy fix'));
subplot(1,2,2);
hold on;
plot(xc,U(2,:),'.-','DisplayName',sprintf('No entropy fixr'));
load R2_1000_0.5_HY.mat
hold on;
subplot(1,2,1);
U = u2U(u,g);
hold on;
plot(xc,U(1,:),'.-','DisplayName',sprintf('With entropy fix'));
subplot(1,2,2);
hold on;
plot(xc,U(2,:),'.-','DisplayName',sprintf('With entropy fix'));
XL = [-0.01,0.01];
subplot(1,2,1);
legend;xlim(XL)
grid on;grid minor;title('\rho, N=1000, 2nd order, t = 0.1');
subplot(1,2,2);
legend;xlim(XL)
grid on;grid minor;title('U, N=1000, 2nd order, t = 0.1');
% ylim([0,inf]);u
%%%%%%%%%%%%%%%%%%

function U = u2U(u,g)
U = u;
U(2,:) =  U(2,:) ./ U(1,:);
U(3,:) = (U(3,:) - 0.5 * U(2,:).^2 .* U(1,:))*(g-1);
end

