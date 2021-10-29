[GLNodes,GLWeights]= F_GLGet(20);

%%
func = @(x) sin(x).^2;
integ = F_GLInt(func,0,pi,GLNodes,GLWeights);
integ-pi/2
