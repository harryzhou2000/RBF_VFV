function elems = ReconstructInit(elems)
for it = 1:numel(elems)
    elems(it).recvalue = zeros(elems(it).dofRec,1);
end
end