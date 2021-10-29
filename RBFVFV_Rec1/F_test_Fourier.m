function [kappas,kappaAs] = F_test_Fourier(options)
h = 0.01;
kappas = (0:100)/100 * pi;
Nmin = 100;

Np = 5*ones(size(kappas));
Ncal = 2*pi./kappas;
Npmin = Nmin./Ncal;
Np = round(max(Npmin,Np));

Nreal = round(2*pi./kappas.*Np);
kappas = 2*pi*Np./Nreal;
kappaAs = nan(size(kappas));
Nreal(1) = 0;
kappaAs(1) = 1;

%%
[N,jmax] = max(Nreal,[],'all','linear');
xr = [0,N*h];
km = kappas(jmax)/h;

xends = linspace(xr(1),xr(2),N+1);
elem0s = GetElements(N);
elem0s = F_Element_Init_1D(elem0s,xends);
[elem0s,basis] = RecMat_1D_RBFVFV(elem0s,options);
condA = cond(elem0s(1).A);

for jt = 2:numel(kappas)
    
    N = round(2*pi/kappas(jt)* Np(jt)) ;
    xr = [0,N*h];
    km = kappas(jt)/h;
    xends = linspace(xr(1),xr(2),N+1);
    elems = elem0s(1:N);
    elems = F_Element_Init_1D(elems,xends);
    
    % INITIAL VALUE
    for it = 1:N
        elems(it).value = ((-exp(1j*elems(it).nodexs(1)*km) + exp(1j*elems(it).nodexs(2)*km))/1j/km/elems(it).hx);
    end
    % Reconstruct
    elems = ReconstructInit(elems);
    elems = Reconstruct(elems,options);

    % Flux Calculation
    freal = (exp(1j*xends(2:end) *km) - exp(1j*xends(1:end-1) *km));
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
    fftfnum = fft(fnum);
    fftfreal = fft(freal);    
%     ks = fnum./freal;
    ks = fftfnum(Np(jt)+1)/fftfreal(Np(jt)+1);
    kappaAs(jt) = mean(ks);
    
    fprintf('---jt = %d, kappaA = %e\n',jt,kappaAs(jt));
end 