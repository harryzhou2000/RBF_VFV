% [GLNodes,GLWeights]= F_GLGet(20);
kp = 1; %多项式阶数
kr = 2; %rbf阶数
nb = 2; %rbf基点数
hdc = 2; %h/c c是rbf尺度
rbftype = 1; %1=mq，2=高斯
ndof = kp + kr*nb; %每个单元的重构自由度
ndir = ndof; %重构泛函所用导数阶数
basepoints = [-0.2 0.2];
weights = ones(1,ndir+1);
weights(1) = 1;
weights(2) = 1;



dfrbfs = cell(ndof,1);
frbfs  = cell(ndof,1);
frbf0s = cell(ndof,1);


for it = 1:kp
    [dfrbfs{it},frbfs{it},frbf0s{it}] = F_GetRbfBasis(0,it,0,1,ndir);
end
for it = 1:kr
    for jt = 1:nb
        [dfrbfs{kp + jt+(it-1)*kr},frbfs{kp + jt+(it-1)*kr},frbf0s{kp + jt+(it-1)*kr}] =...
            F_GetRbfBasis(rbftype,it-1,basepoints(jt),hdc,ndir);
    end
end

%%
h = 0.01;
rs = -h/2:0.001:h/2;
plot(rs,frbfs{2}(rs,h));
frbf0s{2}(h)