function [elems,basis,levels] = RecMat_1D_RBFVFV(elems,options)
nGLNode = 21;
[GLNodes,GLWeights]= F_GLGet(nGLNode-1);
if(nargin == 2)
    kp = options.kp; %多项式阶数
    kr = options.kr; %rbf阶数
    nb = options.nb; %rbf基点数
    hdc = options.hdc; %h/c c是rbf尺度 condA optimized:2 color opt: 1.14
    rbftype = options.rbftype; %1=mq，2=高斯
    w0 = options.w0;%零阶和一阶导数的权重
    basepoints = options.basepoints;
    extradirs  = options.extradirs; %泛函所用界面导数阶数比自由度多的数目
    extend = options.extend;
else %默认
    kp = 3; %多项式阶数
    kr = 0; %rbf阶数
    nb = 2; %rbf基点数
    hdc = 1.14; %h/c c是rbf尺度 condA optimized:2 color opt: 1.14
    rbftype = 1; %1=mq，2=高斯
    w0 = 5;%零阶和一阶导数的权重
    basepoints = [-0.2 0.2];
    extradirs  = 1; %泛函所用界面导数阶数比自由度多的数目
    extend = 1;
end
levels = [1:kp,kp+(1:kr)*nb];
% levels = [1:kp,kp+(1:kr*nb)];
ndof = kp + kr*nb; %每个单元的重构自由度
ndir = ndof+extradirs ; %重构泛函所用导数阶数

weights = ones(1,ndir+1);
weights(1) = w0;
if(ndir>0)
    weights(2) = w0;
end

%% Symbolic Calculation of VFV
dfrbfs = cell(ndof,1);
frbfs  = cell(ndof,1);
frbf0s = cell(ndof,1);

dfrbfsL = cell(ndof,1);
frbfsL  = cell(ndof,1);
frbf0sL = cell(ndof,1);
dfrbfsR = cell(ndof,1);
frbfsR  = cell(ndof,1);
frbf0sR = cell(ndof,1);

for it = 1:kp
    [dfrbfs{it},frbfs{it},frbf0s{it}] = F_GetRbfBasis(0,it,0,1,ndir);
    [dfrbfsL{it},frbfsL{it},frbf0sL{it}] = F_GetRbfBasis_Bias(0,it,0,1,ndir,-1);
    [dfrbfsR{it},frbfsR{it},frbf0sR{it}] = F_GetRbfBasis_Bias(0,it,0,1,ndir,1);
end
for it = 1:kr
    for jt = 1:nb
        [dfrbfs{kp + jt+(it-1)*nb},frbfs{kp + jt+(it-1)*nb},frbf0s{kp + jt+(it-1)*nb}] =...
            F_GetRbfBasis(rbftype,it-1,basepoints(jt),hdc,ndir);
        [dfrbfsL{kp + jt+(it-1)*nb},frbfsL{kp + jt+(it-1)*nb},frbf0sL{kp + jt+(it-1)*nb}] =...
            F_GetRbfBasis_Bias(rbftype,it-1,basepoints(jt),hdc,ndir,-1);
        [dfrbfsR{kp + jt+(it-1)*nb},frbfsR{kp + jt+(it-1)*nb},frbf0sR{kp + jt+(it-1)*nb}] =...
            F_GetRbfBasis_Bias(rbftype,it-1,basepoints(jt),hdc,ndir,1);
    end
end

basis.ndof = ndof;
basis.dfrbfs = dfrbfs;
basis.frbfs = frbfs;
basis.frbf0s = frbf0s;
basis.dfrbfsL = dfrbfsL;
basis.frbfsL = frbfsL;
basis.frbf0sL = frbf0sL;
basis.dfrbfsR = dfrbfsR;
basis.frbfsR = frbfsR;
basis.frbf0sR = frbf0sR;

%% Integrals of VFV
for it = 1:numel(elems)
    elems(it).dofRec = ndof;
    h = elems(it).hx;
    left = elems(it).neighbours(1);
    right = elems(it).neighbours(2);
    hl = elems(left).hx;
    hr = elems(right).hx;
    dleft = (elems(left).hx + elems(it).hx)/2;
    dright = (elems(right).hx + elems(it).hx)/2;
    
    rbfszbRight = nan(ndof,1);
    rbfszbLeft  = nan(ndof,1);
    drbfszbRight = nan(ndof,ndir);
    drbfszbLeft  = nan(ndof,ndir);
    for mt = 1:ndof
        rbfszbRight(mt,1) = frbfs{mt}( h/2,h) - frbf0s{mt}(h);
        rbfszbLeft(mt,1)  = frbfs{mt}(-h/2,h) - frbf0s{mt}(h);
        drbfszbRight(mt,:) = dfrbfs{mt}( h/2,h);
        drbfszbLeft(mt,:)  = dfrbfs{mt}(-h/2,h);
    end
    dsRight = [rbfszbRight,drbfszbRight];
    dsLeft  = [rbfszbLeft ,drbfszbLeft ];
    % 此时：dsleft(m,n): 线元左侧界面，基函数（0均值）m 的 n阶导，x正为正
    % 由于一维，界面上无需积分，按照matlab风格计算矩阵后就是积分值
    % 高维（如用Gauss点积分则在每个积分点上操作后Gauss加和）
    factLeft = dleft.^((0:ndir) -0.5)./factorial(0:ndir);
    factRight = dright.^((0:ndir) -0.5)./factorial(0:ndir);
    elems(it).MLeft = dsLeft.*factLeft.*weights;
    elems(it).MRight = dsRight.*factRight.*weights;
    
%     % left bias
%     rbfszbRightL = nan(ndof,1);
%     rbfszbLeftL  = nan(ndof,1);
%     drbfszbRightL = nan(ndof,ndir);
%     drbfszbLeftL  = nan(ndof,ndir);
%     for mt = 1:ndof
%         rbfszbRightL(mt,1) = frbfsL{mt}( h/2,h,hl) - frbf0sL{mt}(h,hl);
%         rbfszbLeftL(mt,1)  = frbfsL{mt}(-h/2,h,hl) - frbf0sL{mt}(h,hl);
%         drbfszbRightL(mt,:) = dfrbfsL{mt}( h/2,h,hl);
%         drbfszbLeftL(mt,:)  = dfrbfsL{mt}(-h/2,h,hl);
%     end
%     dsRightL = [rbfszbRightL,drbfszbRightL];
%     dsLeftL  = [rbfszbLeftL ,drbfszbLeftL ];
%     elems(it).MLeftL = dsLeftL.*factLeft.*weights;
%     elems(it).MRightL = dsRightL.*factRight.*weights;
%     %%%%%
%     
%     % right bias
%     rbfszbRightR = nan(ndof,1);
%     rbfszbLeftR  = nan(ndof,1);
%     drbfszbRightR = nan(ndof,ndir);
%     drbfszbLeftR  = nan(ndof,ndir);
%     for mt = 1:ndof
%         rbfszbRightR(mt,1) = frbfsR{mt}( h/2,h,hr) - frbf0sR{mt}(h,hr);
%         rbfszbLeftR(mt,1)  = frbfsR{mt}(-h/2,h,hr) - frbf0sR{mt}(h,hr);
%         drbfszbRightR(mt,:) = dfrbfsR{mt}( h/2,h,hr);
%         drbfszbLeftR(mt,:)  = dfrbfsR{mt}(-h/2,h,hr);
%     end
%     dsRightR = [rbfszbRightR,drbfszbRightR];
%     dsLeftR  = [rbfszbLeftR ,drbfszbLeftR ];
%     elems(it).MLeftR = dsLeftR.*factLeft.*weights;
%     elems(it).MRightR = dsRightR.*factRight.*weights;
%     %%%%%
    
    elems(it).A = elems(it).MLeft*transpose(elems(it).MLeft) + elems(it).MRight*transpose(elems(it).MRight);
    elems(it).Ainv = inv(elems(it).A);
    elems(it).bw = [rbfszbLeft/dleft,rbfszbRight/dright]*weights(1)^2;%通过左右均值跳量计算b的权重（注意有权重）
    
%     %left
%     elems(it).AL = elems(it).MLeftL*transpose(elems(it).MLeftL) + elems(it).MRightL*transpose(elems(it).MRightL);
%     elems(it).AinvL = inv(elems(it).AL);
%     elems(it).bwL = [rbfszbLeftL/dleft,rbfszbRightL/dright]*weights(1)^2;%通过左右均值跳量计算b的权重（注意有权重）
%     %right
%     elems(it).AR = elems(it).MLeftR*transpose(elems(it).MLeftR) + elems(it).MRightR*transpose(elems(it).MRightR);
%     elems(it).AinvR = inv(elems(it).AR);
%     elems(it).bwR = [rbfszbLeftR/dleft,rbfszbRightR/dright]*weights(1)^2;%通过左右均值跳量计算b的权重（注意有权重）
end



for it = 1:numel(elems)
    elems(it).dofRec = ndof;
%     h = elems(it).hx;
    left = elems(it).neighbours(1);
    right = elems(it).neighbours(2);
    elems(it).B{1} = elems(it).MLeft  * transpose(elems(left).MRight);
    elems(it).B{2} = elems(it).MRight * transpose(elems(right).MLeft);
%     elems(it).BL{1} = elems(it).MLeftL  * transpose(elems(left).MRight);
%     elems(it).BL{2} = elems(it).MRightL * transpose(elems(right).MLeft);
%     elems(it).BR{1} = elems(it).MLeftR  * transpose(elems(left).MRight);
%     elems(it).BR{2} = elems(it).MRightR * transpose(elems(right).MLeft);
%     elems(it).B{1} = elems(it).B{1};
%     elems(it).B{2} = elems(it).B{2};
end

switch extend
    case 1
        %% Coping with the LSQ matrices %李万爱
        for it = 1:numel(elems)
            elems(it).dofRec = ndof;
            h = elems(it).hx;
            left = elems(it).neighbours(1);
            right = elems(it).neighbours(2);
            hl = elems(left).hx;
            hr = elems(right).hx;
            dleft = (elems(left).hx + elems(it).hx)/2;
            dright = (elems(right).hx + elems(it).hx)/2;
            
            
            BlsqL =  nan(nGLNode,ndof);
            BlsqR =  nan(nGLNode,ndof);
            Alsq = nan(nGLNode,ndof);
            for dof = 1:ndof
                BlsqL(:,dof) = (frbfsL{dof}(GLNodes/2*h,h,hl)-frbf0sL{dof}(h,hl)) .*GLWeights;
                BlsqR(:,dof) = (frbfsR{dof}(GLNodes/2*h,h,hr)-frbf0sR{dof}(h,hr)) .*GLWeights;
                Alsq(:,dof) = (frbfs{dof}(GLNodes/2*h,h)-frbf0s{dof}(h)) .*GLWeights;
            end
            elems(it).TlsqL = BlsqL\Alsq;
            elems(it).TlsqR = BlsqR\Alsq;
            %             elems(it).TlsqL = inv(elems(it).TlsqL);
            %             elems(it).TlsqR = inv(elems(it).TlsqR);
        end
    case 2
        %% Coping with the LSQ matrices %李万爱++
        for it = 1:numel(elems)
            elems(it).dofRec = ndof;
            h = elems(it).hx;
            left = elems(it).neighbours(1);
            right = elems(it).neighbours(2);
            hl = elems(left).hx;
            hr = elems(right).hx;
            dleft = (elems(left).hx + elems(it).hx)/2;
            dright = (elems(right).hx + elems(it).hx)/2;
            
            
            AlsqL =  nan(nGLNode,ndof);
            AlsqR =  nan(nGLNode,ndof);
            Blsq = nan(nGLNode,ndof);
            AMlsqL = nan(nGLNode,ndof);
            AMlsqR = nan(nGLNode,ndof);
            for dof = 1:ndof
                AlsqL(:,dof) = (frbfs{dof}(GLNodes/2*hl + dleft, hl)-frbf0s{dof}(hl)) .*GLWeights;
                AlsqR(:,dof) = (frbfs{dof}(GLNodes/2*hr - dright,hr)-frbf0s{dof}(hr)) .*GLWeights;
                Blsq(:,dof) =  (frbfs{dof}(GLNodes/2*h,h)-frbf0s{dof}(h)) .*GLWeights;
                AMlsqL(:,dof) = (frbfs{dof}(GLNodes/2*hl,hl)-frbf0s{dof}(hl)) .*GLWeights;
                AMlsqR(:,dof) = (frbfs{dof}(GLNodes/2*hr,hr)-frbf0s{dof}(hr)) .*GLWeights;
            end
            
            elems(it).TlsqL = [AlsqL\Blsq,AMlsqL\AMlsqL];
            elems(it).TlsqR = [AlsqR\Blsq,AMlsqR\AMlsqR];
            elems(it).UlsqL = AlsqL\GLWeights;
            elems(it).UlsqR = AlsqR\GLWeights;
        end
    case 3
        %% Coping with the LSQ matrices %李万爱 inverse
        for it = 1:numel(elems)
            elems(it).dofRec = ndof;
            h = elems(it).hx;
            left = elems(it).neighbours(1);
            right = elems(it).neighbours(2);
            hl = elems(left).hx;
            hr = elems(right).hx;
            dleft = (elems(left).hx + elems(it).hx)/2;
            dright = (elems(right).hx + elems(it).hx)/2;
            
            
            BlsqL =  nan(nGLNode,ndof);
            BlsqR =  nan(nGLNode,ndof);
            Alsq = nan(nGLNode,ndof);
            for dof = 1:ndof
                BlsqL(:,dof) = (frbfsL{dof}(GLNodes/2*h,h,hl)-frbf0sL{dof}(h,hl)) .*GLWeights;
                BlsqR(:,dof) = (frbfsR{dof}(GLNodes/2*h,h,hr)-frbf0sR{dof}(h,hr)) .*GLWeights;
                Alsq(:,dof) = (frbfs{dof}(GLNodes/2*h,h)-frbf0s{dof}(h)) .*GLWeights;
            end
            elems(it).TlsqL = BlsqL\Alsq;
            elems(it).TlsqR = BlsqR\Alsq;
            elems(it).TlsqL = inv(elems(it).TlsqL);
            elems(it).TlsqR = inv(elems(it).TlsqR);
        end
        
    case 4
        %% Coping with the LSQ matrices %李万爱 inverse2
        for it = 1:numel(elems)
            elems(it).dofRec = ndof;
            h = elems(it).hx;
            left = elems(it).neighbours(1);
            right = elems(it).neighbours(2);
            hl = elems(left).hx;
            hr = elems(right).hx;
            dleft = (elems(left).hx + elems(it).hx)/2;
            dright = (elems(right).hx + elems(it).hx)/2;
            
            
            BlsqL =  nan(nGLNode,ndof);
            BlsqR =  nan(nGLNode,ndof);
            Alsq = nan(nGLNode,ndof);
            for dof = 1:ndof
                BlsqL(:,dof) = (frbfs{dof}(GLNodes/2*h + dleft,hl)-frbf0s{dof}(hl)) .*GLWeights;
                BlsqR(:,dof) = (frbfs{dof}(GLNodes/2*h - dright,hr)-frbf0s{dof}(hr)) .*GLWeights;
                Alsq(:,dof) = (frbfs{dof}(GLNodes/2*h,h)-frbf0s{dof}(h)) .*GLWeights;
            end
            elems(it).TlsqL = BlsqL\Alsq;
            elems(it).TlsqR = BlsqR\Alsq;
            elems(it).TlsqL = inv(elems(it).TlsqL);
            elems(it).TlsqR = inv(elems(it).TlsqR);
        end
end
if options.constbound
    elems(1).bw =  elems(1).bw*0;
    elems(end).bw =  elems(end).bw*0;
    elems(1).Ainv =  elems(1).Ainv*0;
    elems(end).Ainv =  elems(end).Ainv*0;
%     for i = 1:2
%         elems(1).B{i} =  elems(1).B{i} *0;
%         elems(end).B{i}  =  elems(end).B{i} *0;
%     end
end
elems = rmfield(elems,'MLeft');
elems = rmfield(elems,'MRight');
elems = rmfield(elems,'stencils');
end