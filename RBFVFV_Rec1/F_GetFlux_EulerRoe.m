function [f1,f2,f3,dtm] = F_GetFlux_EulerRoe(u1_elems,u2_elems,u3_elems,basis,gamma)
N = numel(u1_elems);

uleft  = nan(3,N); %left of ith start
uright = nan(3,N);
for it = 1:N
    uleft(:,it) = [u1_elems(it).value;u2_elems(it).value;u3_elems(it).value];
    for kt = 1:u1_elems(it).dofRec
        uleft(:,it) = uleft(:,it) +...
            [u1_elems(it).recvalue(kt)*(basis.frbfs{kt}(u1_elems(it).hx/2 ,u1_elems(it).hx) - basis.frbf0s{kt}(u1_elems(it).hx));...
             u2_elems(it).recvalue(kt)*(basis.frbfs{kt}(u2_elems(it).hx/2 ,u2_elems(it).hx) - basis.frbf0s{kt}(u2_elems(it).hx));...
             u3_elems(it).recvalue(kt)*(basis.frbfs{kt}(u3_elems(it).hx/2 ,u3_elems(it).hx) - basis.frbf0s{kt}(u3_elems(it).hx))];
    end
    uright(:,it) = [u1_elems(it).value;u2_elems(it).value;u3_elems(it).value];
    for kt = 1:u1_elems(it).dofRec
        uright(:,it) = uright(:,it) +...
            [u1_elems(it).recvalue(kt)*(basis.frbfs{kt}(-u1_elems(it).hx/2 ,u1_elems(it).hx) - basis.frbf0s{kt}(u1_elems(it).hx));...
             u2_elems(it).recvalue(kt)*(basis.frbfs{kt}(-u2_elems(it).hx/2 ,u2_elems(it).hx) - basis.frbf0s{kt}(u2_elems(it).hx));...
             u3_elems(it).recvalue(kt)*(basis.frbfs{kt}(-u3_elems(it).hx/2 ,u3_elems(it).hx) - basis.frbf0s{kt}(u3_elems(it).hx))];
    end
end
uleft = [uleft(:,end),uleft(:,1:end-1)]; %[rho rhou e]' left
[f,vmax] = RoeSolver(uleft,uright,gamma);
%%%%
% incU = uright - uleft;
% a2 = (gamma-1)./asmid .* ((Hmid - qsmid ) .* incU(1,:) + umid.*incU(2,:) - incU(3,:));
% a1 = 0.5 * ( incU(1,:)-a2 - (incU(2,:) - umid.*incU(1,:))./sqrt(asmid));
% a3 = 0.5 * ( incU(1,:)-a2 - (incU(2,:) + umid.*incU(1,:))./sqrt(asmid));
% 
% pleft = (uleft(3,:) - 0.5 * uleft(2,:).^2./uleft(1,:))*(gamma-1);
% fleft = [uleft(2,:);uleft(2,:).^2./uleft(1,:)+pleft;uleft(2,:)./uleft(1,:).*(pleft+uleft(3,:))];
% pright = (uright(3,:) - 0.5 * uright(2,:).^2./uright(1,:))*(gamma-1);
% fright = [uright(2,:);uright(2,:).^2./uright(1,:)+pright;uright(2,:)./uright(1,:).*(pright+uright(3,:))];
% f = 0.5*(fleft + fright) - 0.5*([ones(size(umid));Lambda1;Hmid - umid.*sqrt(asmid)] .* a1.*abs(Lambda1) + ...
%     [ones(size(umid));umid;0.5*qsmid] .* a2.*abs(umid) +...
%     [ones(size(umid));Lambda3;Hmid + umid.*sqrt(asmid)] .* a3.*abs(umid));

%%%%
f = f./[u1_elems.hx];
f = [-diff(f,1,2),f(:,end)-f(:,1)];
dt = [u1_elems.hx]./vmax;
dtm = min(dt,[],'all');


f(:,1) = 0;
f(:,end) = 0;
f1 = f(1,:);
f2 = f(2,:);
f3 = f(3,:);



% f1 = -f1;
% f2 = -f2;
% f3 = -f3;
% f1(1) = 0;
% f2(1) = 0;
% f3(1) = 0;

