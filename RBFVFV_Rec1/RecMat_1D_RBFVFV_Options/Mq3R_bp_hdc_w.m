function options = Mq3R_bp_hdc_w(bp,hdc,w)
options.kp = 1; %����ʽ����
options.kr = 2; %rbf����
options.nb = 1; %rbf������
options.hdc = hdc; %h/c c��rbf�߶� condA optimized:2 color opt: 1.14
options.rbftype = 1; %1=mq��2=��˹
options.w0 = w;%��׺�һ�׵�����Ȩ��
options.basepoints = [0];
options.extradirs  = 0; %�������ý��浼�����������ɶȶ����Ŀ
options.extend = 1;
options.constbound = false;
end