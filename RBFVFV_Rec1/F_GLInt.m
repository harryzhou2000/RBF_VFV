function integ = F_GLInt(f,x0,x1,Nodes,Weights)
fs = f(Nodes/2*(x1-x0) + (x1+x0)/2);
integ = Weights'*fs*(x1-x0)/2;