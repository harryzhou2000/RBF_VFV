function [A,Bw] = F_BuildMat(elems)
% A = sparse(elems(1).dofrec*numel(elems), elems(1).dofrec*numel(elems));
nelem = elems(1).dofRec^2*numel(elems)*3;
indi = nan(nelem,1);
indj = nan(nelem,1);
v = nan(nelem,1);
is = (1:elems(1).dofRec)'.*ones(1,elems(1).dofRec);
js = is';
kstart = 1;
sz = elems(1).dofRec^2;
szs = elems(1).dofRec;


for it = 1:numel(elems)
    v(kstart:kstart+sz-1) = elems(it).A(:);
    indi(kstart:kstart+sz-1) = is(:) + szs*(it-1);
    indj(kstart:kstart+sz-1) = js(:) + szs*(it-1);
    kstart = kstart + sz;
    for n = 1:numel(elems(it).neighbours)
        nb=  elems(it).neighbours(n);
        v(kstart:kstart+sz-1) = -elems(it).B{n}(:);
        indi(kstart:kstart+sz-1) = is(:) + szs*(it-1);
        indj(kstart:kstart+sz-1) = js(:) + szs*(nb-1);
        kstart = kstart + sz;
    end
end



A = sparse(indi,indj,v);

nelem = elems(1).dofRec*numel(elems)*2;
indi = nan(nelem,1);
indj = nan(nelem,1);
v = nan(nelem,1);
is = 1:elems(1).dofRec;
js = ones(size(is));
kstart = 1;
sz = elems(1).dofRec;
szs = elems(1).dofRec;


for it = 1:numel(elems)
    v(kstart:kstart+sz-1) = -sum(elems(it).bw,2);
    indi(kstart:kstart+sz-1) = is(:) + szs*(it-1);
    indj(kstart:kstart+sz-1) = js(:) + (it-1);
    kstart = kstart + sz;
    for n = 1:numel(elems(it).neighbours)
        nb=  elems(it).neighbours(n);
        v(kstart:kstart+sz-1) = elems(it).bw(:,n);
        indi(kstart:kstart+sz-1) = is(:) + szs*(it-1);
        indj(kstart:kstart+sz-1) = js(:) + (nb-1);
        kstart = kstart + sz;
    end
end

Bw = sparse(indi,indj,v);