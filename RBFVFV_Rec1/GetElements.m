function elems = GetElements(m)

neighbours = cell(m,1);
stencils = cell(m,1);
A = cell(m,1);
B = cell(m,1);
bw = cell(m,1); %weights for b
hx = cell(m,1);
nodexs = cell(m,1);
centx = cell(m,1);
dofRec = cell(m,1);
value = cell(m,1);
recvalue = cell(m,1);
elems = struct('neighbours',neighbours,'stencils',stencils,'A',A,'B',B,'bw',bw,'hx',hx,'nodexs',nodexs,'centx',centx,...
    'dofRec',dofRec,'value',value,'recvalue',recvalue);

