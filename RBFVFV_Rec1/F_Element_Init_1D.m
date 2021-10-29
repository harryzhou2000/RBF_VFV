function elems = F_Element_Init_1D(elems,xends)

for it = 1:numel(elems)
    elems(it).nodexs = [xends(it),xends(it+1)];
    elems(it).centx  = (max(elems(it).nodexs)+min(elems(it).nodexs))/2;
    elems(it).hx     = max(elems(it).nodexs)-min(elems(it).nodexs);
    elems(it).neighbours = [mod(it-2,numel(elems))+1,mod(it,numel(elems))+1];
    elems(it).stencils = [mod(it-2,numel(elems))+1,mod(it,numel(elems))+1];
end