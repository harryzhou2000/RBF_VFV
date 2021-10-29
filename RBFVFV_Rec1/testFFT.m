
xs = linspace(0,2*pi*(1-1/1024),1024);
ys = exp(100i*xs);
plot(abs(fft(ys)));

%%
ysfft = zeros(1,128);
xs = linspace(0,2*pi*(1-1/128),128);
k = 1;
ysfft(k+1) = 128;
ys = ifft(ysfft);
plot(imag(ys));
hold on;
plot(sin(k*xs));
hold off;