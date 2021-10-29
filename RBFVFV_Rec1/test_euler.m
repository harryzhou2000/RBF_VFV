clear;
Ne = 1000;
CFL = 0.5;
xs = linspace(0,1,Ne+1);
h = xs(2) - xs(1);
dt = 0.001;

gamma = 1.4;
% options = Poly3_w(5);
options = Mq3_bp_hdc_w(1/3,2.2,6);
options.extend = 3;
options.constbound = true;

u1_elems = GetElements(Ne);
u1_elems = F_Element_Init_1D(u1_elems,xs);
[u1_elems,basis,levels] = RecMat_1D_RBFVFV(u1_elems, options);
u2_elems = GetElements(Ne);
u2_elems = F_Element_Init_1D(u2_elems,xs);
[u2_elems,basis,levels] = RecMat_1D_RBFVFV(u2_elems, options);
u3_elems = GetElements(Ne);
u3_elems = F_Element_Init_1D(u3_elems,xs);
[u3_elems,basis,levels] = RecMat_1D_RBFVFV(u3_elems, options);
[A,Bw] = F_BuildMat(u1_elems);
K = 1;
u1_init = (([u1_elems.centx]>0.5)*0.5   +   ([u1_elems.centx]<=0.5)*0.445);
u2_init = (([u1_elems.centx]>0.5)*0   +     ([u1_elems.centx]<=0.5)*0.698) .* u1_init * 0;
u3_init = (([u1_elems.centx]>0.5)*0.571   + ([u1_elems.centx]<=0.5)*3.528)./(gamma-1) + 0.5 * u2_init.^2 ./ u1_init;

% u_init = sin([u_elems.centx] * 10 * pi);
for it = 1:numel(u1_elems)
   u1_elems(it).value = u1_init(it); 
   u2_elems(it).value = u2_init(it); 
   u3_elems(it).value = u3_init(it); 
end
u1_elems = ReconstructInit(u1_elems);
u1_elems = Reconstruct_Sparse(u1_elems,A,Bw,options,levels);
u2_elems = ReconstructInit(u2_elems);
u2_elems = Reconstruct_Sparse(u2_elems,A,Bw,options,levels);
u3_elems = ReconstructInit(u3_elems);
u3_elems = Reconstruct_Sparse(u3_elems,A,Bw,options,levels);
% figure(1);
% 
% F_plotRecs(u1_elems,basis);
% figure(2);
% F_plotRecs(u2_elems,basis);
% figure(3);
% F_plotRecs(u3_elems,basis);
% plot([u1_elems.centx],[u1_elems.value],[u2_elems.centx],[u2_elems.value]./[u1_elems.value]);
ue{1} = u1_elems;
ue{2} = u2_elems;
ue{3} = u3_elems;
%%
% 
t = 0;
tm = 0.1;
ifend = false;
% for iter = 1:10000
%     clc;
%     [k1_1,k1_2,k1_3,dtm] = F_GetFlux_EulerRoe(u1_elems,u2_elems,u3_elems,basis,gamma);
%     dt = dtm*CFL;
%     if t + dt >=tm
%         dt = tm-dt;
%         ifend = true;
%     end
%     u2_1 = Reconstruct( F_Element_Increment(u1_elems,dt*k1_1/2),levels);
%     u2_2 = Reconstruct( F_Element_Increment(u2_elems,dt*k1_2/2),levels);
%     u2_3 = Reconstruct( F_Element_Increment(u3_elems,dt*k1_3/2),levels);
%     
%     [k2_1,k2_2,k2_3] = F_GetFlux_EulerRoe(u2_1,u2_2,u2_3,basis,gamma);
%     u2_1 = Reconstruct( F_Element_Increment(u1_elems,dt*k2_1/2),levels);
%     u2_2 = Reconstruct( F_Element_Increment(u2_elems,dt*k2_2/2),levels);
%     u2_3 = Reconstruct( F_Element_Increment(u3_elems,dt*k2_3/2),levels);
%     
%     [k3_1,k3_2,k3_3] = F_GetFlux_EulerRoe(u2_1,u2_2,u2_3,basis,gamma);
%     u2_1 = Reconstruct( F_Element_Increment(u1_elems,dt*k3_1  ),levels);
%     u2_2 = Reconstruct( F_Element_Increment(u2_elems,dt*k3_2  ),levels);
%     u2_3 = Reconstruct( F_Element_Increment(u3_elems,dt*k3_3  ),levels);
%     
%     [k4_1,k4_2,k4_3] = F_GetFlux_EulerRoe(u2_1,u2_2,u2_3,basis,gamma);
%     u1_elems = Reconstruct( F_Element_Increment(u1_elems, dt/6*(k1_1 + 2*k2_1 + 2*k3_1 + k4_1)),levels);
%     u2_elems = Reconstruct( F_Element_Increment(u2_elems, dt/6*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2)),levels);
%     u3_elems = Reconstruct( F_Element_Increment(u3_elems, dt/6*(k1_3 + 2*k2_3 + 2*k3_3 + k4_3)),levels);
%     
%     t = t+dt;
%     plot([u1_elems.centx],[u1_elems.value],[u2_elems.centx],[u2_elems.value]./[u1_elems.value]);
%     title(sprintf('t=%e',t));
%     fprintf('Iter = %d, at = %f\n',iter,t);
%     drawnow;
%     if(ifend)
%         break;
%     end
% end

for iter = 1:10000
    [k4{1},k4{2},k4{3},dtm] = F_GetFlux_EulerRoe(ue{1},ue{2},ue{3},basis,gamma);
    dt = dtm*CFL;
    if t + dt >=tm
        dt = tm-t;
        ifend = true;
    end
    parfor it = 1:3
        ue{it} = Reconstruct_Sparse( F_Element_Increment(ue{it}, dt*( k4{it})),A,Bw,options,levels);
    end
    
    t = t+dt;
    plot([ue{1}.centx],[ue{1}.value],[ue{2}.centx],[ue{2}.value]./[ue{1}.value]);
    title(sprintf('t=%e',t));
    fprintf('Iter = %d, at = %f\n',iter,t);
    drawnow;
    if(ifend)
        break;
    end
end
