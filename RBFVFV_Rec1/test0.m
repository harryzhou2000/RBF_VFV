clc;
clear;
N = 100;
xr = [0,1];
bp = 1/3;
hdc = 1.8;
w = 5;

xends = linspace(xr(1),xr(2),N+1);
elems = GetElements(N);
elems = F_Element_Init_1D(elems,xends);
% options = Mq3R_bp_hdc_w(bp, hdc,w);
% options = Poly1;
options = Poly3_w(5);
% options = Mq3_bp_hdc_w(bp, hdc,w);
options.extend = 1;

[elems,basis,levels] = RecMat_1D_RBFVFV(elems,options);
condA = cond(elems(1).A)
[A,Bw] = F_BuildMat(elems);

%% INITIAL VALUE
km = 2.345*2*pi;
for it = 1:N
    elems(it).value = real([((-exp(1j*elems(it).nodexs(1)*km) + exp(1j*elems(it).nodexs(2)*km))/1j/km/elems(it).hx)]);
    %         elems(it).value = elems(it).centx^2;
    %     elems(it).value = double(elems(it).centx > 0.5);
end
elems = ReconstructInit(elems);
% for it = 1:2
%     elems = Reconstruct(elems);
% end
b = Bw*([elems.value]');
rec = A\b;

tic;
for it =1:1000
    elems = Reconstruct(elems,options);
end
toc
% elems = Reconstruct(elems,options,levels);
% elems = Reconstruct(elems,options,3);
tic;
for it =1:1000
    elems = Reconstruct_Sparse(elems,A,Bw,options);
end
toc
% elems = Reconstruct_Sparse(elems,A,Bw,options,levels);
% elems = Reconstruct(elems,options,3);

%% plot rec
figure(1)
clf;
hold on;
for it = 1:N
    xp = linspace(-0.5,0.5,20)'*elems(it).hx + elems(it).centx;
    vrec = ones(size(xp,1),1) .* elems(it).value;
    for jt = 1:elems(it).dofRec
        vrec = vrec + elems(it).recvalue(jt)...
            *(basis.frbfs{jt}(xp - elems(it).centx ,elems(it).hx) - basis.frbf0s{jt}(elems(it).hx));
        %         plot(xp,real(vrec));
        %          plot([xp(1),xp(end)],real([vrec(1),vrec(end)]),'ok');
    end
    %     xps = [xps;xp];
    %     vrecs = [vrecs;vrec];
    plot(xp,real(vrec));
    plot([xp(1),xp(end)],real([vrec(1),vrec(end)]),'ok');
end
plot([elems.centx],real([elems.value]),'d');
grid on;
grid minor;
xlim([0,0.1]);

switch options.extend
    case {1,3,4}
        %% plot recL
        figure(2)
        clf;
        hold on;
        for it = 1:N
            xp = linspace(-1.5,1.5,20)'*elems(it).hx + elems(it).centx;
            vrec = ones(size(xp,1),1) .* elems(it).value;
            for jt = 1:elems(it).dofRec
                vrec = vrec + elems(it).recvalueL(jt)...
                    *(basis.frbfsL{jt}(xp - elems(it).centx ,elems(it).hx,elems(it).hx) - basis.frbf0sL{jt}(elems(it).hx,elems(it).hx));
                %         plot(xp,real(vrec));
                %         plot([xp(1),xp(end)],real([vrec(1),vrec(end)]),'ok');
            end
            %     xps = [xps;xp];
            %     vrecs = [vrecs;vrec];
            plot(xp,real(vrec));
            plot([xp(1),xp(end)],real([vrec(1),vrec(end)]),'ok');
        end
        plot([elems.centx],real([elems.value]),'d');
        grid on;
        grid minor;
        xlim([0,0.1]);
        %% plot recR
        figure(3)
        clf;
        hold on;
        for it = 1:N
            xp = linspace(-1.5,1.5,20)'*elems(it).hx + elems(it).centx;
            vrec = ones(size(xp,1),1) .* elems(it).value;
            
            for jt = 1:elems(it).dofRec
                vrec = vrec + elems(it).recvalueR(jt)...
                    *(basis.frbfsR{jt}(xp - elems(it).centx ,elems(it).hx,elems(it).hx) - basis.frbf0sR{jt}(elems(it).hx,elems(it).hx));
                %          plot(xp,real(vrec));
                %          plot([xp(1),xp(end)],real([vrec(1),vrec(end)]),'ok');
                %drawnow;
            end
            %     xps = [xps;xp];
            %     vrecs = [vrecs;vrec];
            plot(xp,real(vrec));
            plot([xp(1),xp(end)],real([vrec(1),vrec(end)]),'ok');
            %drawnow;
        end
        plot([elems.centx],real([elems.value]),'d');
        grid on;
        grid minor;
        xlim([0,0.1]);
        
    case 2
        %% plot recL
        figure(2)
        clf;
        hold on;
        for it = 1:N
            left = elems(it).neighbours(1);
            right = elems(it).neighbours(2);
            xp = linspace(-1.5,0.5,20)'*elems(it).hx + elems(it).centx;
            vrec = ones(size(xp,1),1) .* elems(left).value;
            for jt = 1:elems(it).dofRec
                vrec = vrec + elems(it).recvalueL(jt)...
                    *(basis.frbfs{jt}(xp - elems(it).centx +elems(it).hx,elems(it).hx) - basis.frbf0s{jt}(elems(it).hx));
                %         plot(xp,real(vrec));
                %         plot([xp(1),xp(end)],real([vrec(1),vrec(end)]),'ok');
            end
            %     xps = [xps;xp];
            %     vrecs = [vrecs;vrec];
            plot(xp,real(vrec));
            plot([xp(1),xp(end)],real([vrec(1),vrec(end)]),'ok');
        end
        plot([elems.centx],real([elems.value]),'d');
        grid on;
        grid minor;
        xlim([0,0.1]);
        %% plot recR
        figure(3)
        clf;
        hold on;
        for it = 1:N
            left = elems(it).neighbours(1);
            right = elems(it).neighbours(2);
            xp = linspace(-0.5,1.5,20)'*elems(it).hx + elems(it).centx;
            vrec = ones(size(xp,1),1) .* elems(right).value;
            for jt = 1:elems(it).dofRec
                vrec = vrec + elems(it).recvalueR(jt)...
                    *(basis.frbfs{jt}(xp - elems(it).centx -elems(it).hx,elems(it).hx) - basis.frbf0s{jt}(elems(it).hx));
                %         plot(xp,real(vrec));
                %         plot([xp(1),xp(end)],real([vrec(1),vrec(end)]),'ok');
            end
            %     xps = [xps;xp];
            %     vrecs = [vrecs;vrec];
            plot(xp,real(vrec));
            plot([xp(1),xp(end)],real([vrec(1),vrec(end)]),'ok');
        end
        plot([elems.centx],real([elems.value]),'d');
        grid on;
        grid minor;
        xlim([0,0.1]);
end
%%

figure(1)
