function options = Poly1()
options.kp = 1; %����ʽ����
options.kr = 0; %rbf����
options.nb = 2; %rbf������
options.hdc = 1.14; %h/c c��rbf�߶� condA optimized:2 color opt: 1.14
options.rbftype = 1; %1=mq��2=��˹
options.w0 = 5;%��׺�һ�׵�����Ȩ��
options.basepoints = [-0.2 0.2];
options.extradirs  = 1; %�������ý��浼�����������ɶȶ����Ŀ
options.extend = 1;
options.constbound = false;
end