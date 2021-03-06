function options = Mq5_bp_hdc_w(bp,hdc,w)
options.kp = 3; %多项式阶数
options.kr = 1; %rbf阶数
options.nb = 2; %rbf基点数
options.hdc = hdc; %h/c c是rbf尺度 condA optimized:2 color opt: 1.14
options.rbftype = 1; %1=mq，2=高斯
options.w0 = w;%零阶和一阶导数的权重
options.basepoints = [-bp bp];
options.extradirs  = 0; %泛函所用界面导数阶数比自由度多的数目
options.extend = 1;
options.constbound = false;
end