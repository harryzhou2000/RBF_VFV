function [dfrbf,frbf,frbf0] = F_GetRbfBasis(type,p,basepoint,hdc,ndir)

syms r c h

switch type
    case 0
        frbf = ((r-basepoint*h)/c)^p;
    case 1
        frbf = ((r-basepoint*h)/c)^p*sqrt(1 + (r-basepoint*h)^2/c^2);
    case 2
        frbf = ((r-basepoint*h)/c)^p*exp(-((r-basepoint*h)/c)^2);
end
dfrbf(1) = diff(frbf,r);
for it = 2:ndir
    dfrbf(it) = diff(dfrbf(it-1),r);
end
frbf0 = int(frbf,r);
frbf0 = (subs(frbf0,r,h/2) - subs(frbf0,r,-h/2))/h;
frbf0 = subs(frbf0,c,h/hdc);

dfrbf = subs(dfrbf,c,h/hdc);
frbf = subs(frbf,c,h/hdc);

dfrbf = matlabFunction(symfun(dfrbf,[r,h]));
frbf = matlabFunction(symfun(frbf,[r,h]));
frbf0 = matlabFunction(symfun(frbf0,h));





