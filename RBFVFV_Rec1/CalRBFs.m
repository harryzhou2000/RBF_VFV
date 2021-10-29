syms r c h

frbf = sqrt(1 + r^2/c^2);

%%
hdc = 2;
dfrbf(1) = diff(frbf,r);
for it = 2:8
    dfrbf(it) = diff(dfrbf(it-1),r);
end
dfrbfRight = subs(dfrbf,r,hdc*c/2);
dfrbfRight = subs(dfrbfRight,c,h/hdc);
dfrbfLeft = subs(dfrbf,r,-hdc*c/2);
dfrbfLeft = subs(dfrbfLeft,c,h/hdc);

frbfRight = subs(frbf,r,hdc*c/2);
frbfRight = subs(frbfRight,c,h/hdc);
frbfLeft = subs(frbf,r,-hdc*c/2);
frbfLeft = subs(frbfLeft,c,h/hdc);

frbf0 = int(frbf);
frbf0 = subs(frbf0,r,hdc*c/2) - subs(frbf0,r,-hdc*c/2);
frbf0 = subs(frbf0,c,h/hdc);

