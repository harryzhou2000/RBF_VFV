kappas = linspace(0,pi,10000000);
km1 = (1-exp(-1j*kappas))/1j;
km2 = (3-4*exp(-1j*kappas)+exp(-2*1j*kappas))/2j;
km3 = (2*exp(1j*kappas) + 3 - 6*exp(-1j*kappas) + exp(-2*1j*kappas))/6j;

plot(kappas,real(km1),kappas,real(km2),kappas,real(km3),kappas,kappas)
legend('\kappa_m_1','\kappa_m_2','\kappa_m_3','\kappa','Location','northwest');
xlabel('\kappa');
ylabel('Real(\kappa_m)');

axis equal;
grid on;
grid minor;
xlim([0,pi]);
ylim([0,inf]);


plot(kappas,imag(km1),kappas,imag(km2),kappas,imag(km3),kappas,kappas*0)
legend('\kappa_m_1','\kappa_m_2','\kappa_m_3','0','Location','southwest');
xlabel('\kappa');
ylabel('Imag(\kappa_m)');

axis equal;
grid on;
grid minor;
xlim([0,pi]);
ylim([-inf,inf]);

km1re = abs((real(km1)-kappas)./kappas);
km2re = abs((real(km2)-kappas)./kappas);
km3re = abs((real(km3)-kappas)./kappas);
kappa1g = kappas(km1re<=0.01);
kappa2g = kappas(km2re<=0.01 );
kappa3g = kappas(km3re<=0.01);
plot(kappa2g)
k1top = max(kappa1g);
k2top = max(kappa2g);
k3top = max(kappa3g);
%%
alpha = 0.0537;
beta = 1;
realk5 = alpha*(sin(3*kappas)-4*sin(2*kappas)+5*sin(kappas)) - 1/6*sin(2*kappas) + 4/3 * sin(kappas);
imagk5 = beta*(cos(3*kappas)-6*cos(2*kappas)+15*cos(kappas)-10);
plot(kappas,(realk5-kappas)./kappas);
max((realk5-kappas)./kappas)
max(kappas(abs((realk5-kappas)./kappas)<0.01))
%%
plot(kappas,imagk5)