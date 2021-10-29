clear;
Ne = 100;
CFL = 0.5;
a = -1;
xs = linspace(0,1,Ne+1);
h = xs(2) - xs(1);
dt = abs(a)*CFL * h;
options = Poly3_w(5);
options = Mq3_bp_hdc_w(0.3333,2,5);
% options = Mq5_bp_hdc_w(0.01,2,5);
% options = Mq3R_bp_hdc_w(1/3,10,5);

% options = Mq3_bp_hdc_w(1/3,2.5,10);
% options = Mq3_bp_hdc_w(1/3,2.2,6);
% options = Mq3R_bp_hdc_w(1/3,1.8,5);
% options = Poly1();
options.extend = 3;

u_elems = GetElements(Ne);
u_elems = F_Element_Init_1D(u_elems,xs);
[u_elems,basis,levels] = RecMat_1D_RBFVFV(u_elems, options);
% levels = [];
u_init = double(([u_elems.centx]>0.25) & ([u_elems.centx]<0.75));
% u_init = sin([u_elems.centx] * 2 * pi * 3);
for it = 1:numel(u_elems)
    u_elems(it).value = u_init(it);
end

u_elems = ReconstructInit(u_elems);
u_elems = Reconstruct(u_elems,options,levels);
u_elems = Reconstruct(u_elems,options,levels);
u_elems = Reconstruct(u_elems,options,levels);
F_plotRecs(u_elems,basis);
%%
[A,Bw] = F_BuildMat(u_elems);
for iter = 1:round(Ne/CFL*1)
    clc;
    k1 = a * F_GetFlux_LinConv(u_elems,basis);
    %     u2 = Reconstruct( F_Element_Increment(u_elems,dt*k1/2),options, levels);
    u2 = Reconstruct_Sparse( F_Element_Increment(u_elems,dt*k1/2),A,Bw,options, levels);
    k2 = a * F_GetFlux_LinConv(u2,basis);
    %     u2 = Reconstruct( F_Element_Increment(u_elems,dt*k2/2),options, levels);
    u2 = Reconstruct_Sparse( F_Element_Increment(u_elems,dt*k2/2),A,Bw,options, levels);
    k3 = a * F_GetFlux_LinConv(u2,basis);
    %     u2 = Reconstruct( F_Element_Increment(u_elems,dt*k3  ),options, levels);
    u2 = Reconstruct_Sparse( F_Element_Increment(u_elems,dt*k3  ),A,Bw,options, levels);
    k4 = a * F_GetFlux_LinConv(u2,basis);
    %     u_elems = Reconstruct( F_Element_Increment(u_elems, dt/6*(k1 + 2*k2 + 2*k3 + k4)),options, levels);
    u_elems = Reconstruct_Sparse( F_Element_Increment(u_elems, dt/6*(k1 + 2*k2 + 2*k3 + k4)),A,Bw,options, levels);
    
    %     k1 = a * F_GetFlux_LinConv(u_elems,basis);
    %     u2 = Reconstruct( F_Element_Increment(u_elems,dt*k1/2),options);
    %     k2 = a * F_GetFlux_LinConv(u2,basis);
    %     u2 = Reconstruct( F_Element_Increment(u_elems,dt*k2/2),options);
    %     k3 = a * F_GetFlux_LinConv(u2,basis);
    %     u2 = Reconstruct( F_Element_Increment(u_elems,dt*k3  ),options);
    %     k4 = a * F_GetFlux_LinConv(u2,basis);
    %     u_elems = Reconstruct( F_Element_Increment(u_elems, dt/6*(k1 + 2*k2 + 2*k3 + k4)),options);
    
    plot([u_elems.centx],[u_elems.value],'o-');
    fprintf('Iter = %d, at = %f\n',iter,iter*dt);
    drawnow;
end
%%
hold on;
plot([u_elems.centx],u_init);
hold off;
ylim([0,1])

% for iter = 1:10000
%     k1 = a * F_GetFlux_LinConv(u_elems,basis);
%     u_elems = Reconstruct(F_Element_Increment(u_elems, dt*k1));
%
%     plot([u_elems.centx],[u_elems.value]);
%     fprintf('Iter = %d, at = %f\n',iter,iter*dt);
%     drawnow;
% end
