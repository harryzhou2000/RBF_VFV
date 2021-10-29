function options = Poly3_w5()
options.kp = 3; %多项式阶数
options.kr = 0; %rbf阶数
options.nb = 2; %rbf基点数
options.hdc = 1.14; %h/c c是rbf尺度 condA optimized:2 color opt: 1.14
options.rbftype = 1; %1=mq，2=高斯
options.w0 = 5;%零阶和一阶导数的权重
options.basepoints = [-0.2 0.2];
options.extradirs  = 1; %泛函所用界面导数阶数比自由度多的数目
options.extend = 1;
options.constbound = false;
end