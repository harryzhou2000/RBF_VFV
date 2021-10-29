function fnum = F_GetFlux_LinConv(elems,basis)
N = numel(elems);
fnum  = nan(1,N);
for it = 1:N
    left = elems(it).neighbours(1);
    fL = elems(left).value;
    for kt = 1:elems(left).dofRec
        fL = fL + elems(left).recvalue(kt)...
            *(basis.frbfs{kt}(elems(left).hx/2 ,elems(left).hx) - basis.frbf0s{kt}(elems(left).hx));
    end
    fR = elems(it).value;
    for kt = 1:elems(it).dofRec
        fR = fR + elems(it).recvalue(kt)...
            *(basis.frbfs{kt}(elems(it).hx/2 ,elems(it).hx) - basis.frbf0s{kt}(elems(it).hx));
    end
    fnum(it) = (fR-fL);
end
fnum = fnum./[elems.hx];