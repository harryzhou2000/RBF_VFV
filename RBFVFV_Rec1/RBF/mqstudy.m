h = 1;
c = h*1.1;
bw = 1/3;

f1 = @(x) sqrt(((x+bw)/c).^2+1);
f2 = @(x) sqrt(((x-bw)/c).^2+1);
f0 = @(x) x;
xs = linspace(-0.5,0.5,100);
ca = [0 10 10];
ys = ca(2)*f1(xs) + ca(3)*f2(xs) + ca(1)*f0(xs);
plot(xs,ys);
%%
syms x
fmqa = sqrt((x/c)^2 + 1);
fplot(diff(fmqa,x))
xlim([-1,1])
