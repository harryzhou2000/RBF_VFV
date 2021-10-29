function [f,vmax] = RoeSolver(uleft,uright,gamma)

wleft = U2W(uleft); %sqrt(rho) [1 u H]' left
wright = U2W(uright);
wa = (wleft + wright)*0.5;
umid = wa(2,:)./wa(1,:); 
Hmid = wa(3,:)./wa(1,:);
qsmid = umid.^2;
asmid = ((gamma-1)*(Hmid - 0.5 *qsmid));
if(sum(asmid<0.0000)>0)
   error('imag a'); 
%     asmid(asmid<0) = abs(asmid(asmid<0.0001));
end
% asmid(asmid<0) = 0;

Lambda1 = umid - sqrt(asmid);
Lambda3 = umid + sqrt(asmid); % Lambda2 = umid
vmax = max(abs([Lambda1;Lambda3;umid]),[],1);
% Harten-Yee
deltaEntropy = 0.2*(sqrt(qsmid)+sqrt(asmid));
L1Fix = abs(Lambda1)<deltaEntropy;
L3Fix = abs(Lambda3)<deltaEntropy;
Lambda1(L1Fix) = sign(Lambda1(L1Fix)).*(Lambda1(L1Fix).^2+deltaEntropy(L1Fix).^2)./deltaEntropy(L1Fix)/2;
Lambda3(L3Fix) = sign(Lambda3(L3Fix)).*(Lambda3(L3Fix).^2+deltaEntropy(L3Fix).^2)./deltaEntropy(L3Fix)/2;

% a2Left = (gamma-1)./asmid .* ((Hmid - qsmid ) .* uleft(1,:) + umid.*uleft(2,:) - uleft(3,:));
% a1Left = 0.5 * ( uleft(1,:)-a2Left - (uleft(2,:) - umid.*uleft(1,:))./sqrt(asmid));
% a3Left = 0.5 * ( uleft(1,:)-a2Left - (uleft(2,:)=*hyh + umid.*uleft(1,:))./sqrt(asmid));
% 
% a2Right = (gamma-1)./asmid .* ((Hmid - qsmid ) .* uright(1,:) + umid.*uright(2,:) - uright(3,:));
% a1Right = 0.5 * ( uright(1,:)-a2Right - (uright(2,:) - umid.*uright(1,:))./sqrt(asmid));
% a3Right = 0.5 * ( uright(1,:)-a2Right - (uright(2,:) + umid.*uright(1,:))./sqrt(asmid));
% 
% a1 = a1Left.*(Lambda1>0) + a1Right.*(Lambda1<=0);
% a2 = a2Left.*(umid>0) + a2Right.*(umid<=0);
% a3 = a3Left.*(Lambda3>0) + a3Right.*(Lambda3<=0);
% 
% f1 = a1.*Lambda1 + a2.*umid + a3.*Lambda3;
% f2 = a1.*Lambda1.^2 + a2.*umid.^2 + a3.*Lambda3.^2;
% f3 = a1.*(Hmid - sqrt(asmid).*umid).*Lambda1 + a2.*0.5.*qsmid.*umid + a3.*(Hmid + sqrt(asmid).*umid).*Lambda3;
% 
% f = [f1;f2;f3];
%%%%
incU = uright - uleft;
a2 = (gamma-1)./asmid .* ((Hmid - qsmid ) .* incU(1,:) + umid.*incU(2,:) - incU(3,:));
a1 = 0.5 * ( incU(1,:)-a2 - (incU(2,:) - umid.*incU(1,:))./sqrt(asmid));
a3 = incU(1,:) - a1 - a2;

pleft = (uleft(3,:) - 0.5 * uleft(2,:).^2./uleft(1,:))*(gamma-1);
fleft = [uleft(2,:);uleft(2,:).^2./uleft(1,:)+pleft;uleft(2,:)./uleft(1,:).*(pleft+uleft(3,:))];
pright = (uright(3,:) - 0.5 * uright(2,:).^2./uright(1,:))*(gamma-1);
fright = [uright(2,:);uright(2,:).^2./uright(1,:)+pright;uright(2,:)./uright(1,:).*(pright+uright(3,:))];
f = 0.5*(fleft + fright) - 0.5*([ones(size(umid));umid - sqrt(asmid);Hmid - umid.*sqrt(asmid)] .* a1.*abs(Lambda1) + ...
    [ones(size(umid));umid;0.5*qsmid] .* a2.*abs(umid) +...
    [ones(size(umid));umid + sqrt(asmid);Hmid + umid.*sqrt(asmid)] .* a3.*abs(Lambda3));

%%%

function w = U2W(u)
    w = [sqrt(u(1,:));...
        u(2,:)./sqrt(u(1,:)) ;...
        (u(3,:) + (u(3,:) - 0.5*u(2,:).^2./u(1,:)) *(gamma-1))./sqrt(u(1,:))];
end
end