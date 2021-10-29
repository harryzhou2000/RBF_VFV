function F_plotRecs(elems,basis)
N = numel(elems);
clf;
hold on;
for it = 1:N
    xp = linspace(-0.5,0.5,20)'*elems(it).hx + elems(it).centx;
    vrec = ones(size(xp,1),1) .* elems(it).value;
    for jt = 1:elems(it).dofRec
        vrec = vrec + elems(it).recvalue(jt)...
            *(basis.frbfs{jt}(xp - elems(it).centx ,elems(it).hx) - basis.frbf0s{jt}(elems(it).hx));
    end
    %     xps = [xps;xp];
    %     vrecs = [vrecs;vrec];
    plot(xp,real(vrec));
    plot([xp(1),xp(end)],real([vrec(1),vrec(end)]),'ok');
end
plot([elems.centx],real([elems.value]),'d');
grid on;
grid minor;
hold off;