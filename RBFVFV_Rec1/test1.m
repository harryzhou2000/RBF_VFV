figure(1);clf;figure(2);clf;
%%

bp = 1/3;
hdc = 2.2;
w = 6;
options = Poly3_w(5); Name = sprintf('Poly3 w=%.1f',5);
% options = Mq3_bp_hdc_w(bp, hdc ,w); Name = sprintf('Mq3-%.2f-%.2f-%.1f',bp,hdc,w);
options = Mq5R_bp_hdc_w(1/3, 1.8 ,5); Name = sprintf('Mq5R-%.2f-%.2f-%.1f',bp,hdc,w);
options = Mq5_bp_hdc_w(1/3, 1.8 ,5); Name = sprintf('Mq5-%.2f-%.2f-%.1f',bp,hdc,w);
% options = Mq3R_bp_hdc_w(bp, hdc ,w); Name = sprintf('Mq3R-%.2f-%.2f-%.1f',bp,hdc,w);
%%%%%%%%%%%%%% 1/3 2 5
%%%%%%%%%%%%%% 1/3 2.5 10
% bp 在1/3时，相对色散都是负，更大时有正，色散较优
%%%%%%%%%%%%%%
[kappas, kappaAs] = F_test_Fourier(options);
%
figure(1)
hold on;
% plot(kappas,abs(real(kappaAs)-1))
% plot(kappas,(real(kappaAs)-1))
plot(kappas,(real(kappaAs)).*kappas,'DisplayName',Name)
xlim([0,pi]);
axis equal;
grid on;
grid minor;
legend;
%
figure(2)
hold on;
plot(kappas,imag(kappaAs).*kappas,'DisplayName',Name)
xlim([0,pi]);
grid on;
grid minor;
legend;
%
figure(1)
%%
figure(1)
% title('色散曲线');
xlabel('\kappa');
ylabel('Real(\kappa'')');
figure(2)
% title('耗散曲线');
xlabel('\kappa');
ylabel('Imag(\kappa'')');