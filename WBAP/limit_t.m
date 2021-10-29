function ulim = limit_t(u,uS,eps,N)
    dim = size(u,1);
    Nu = sqrt(sum(u.^2,1));

    theta = cat(2,reshape(u,dim,1,[]),reshape(uS,dim,size(uS,2),[]))./reshape(Nu,1,1,[]);
    Btheta = theta./sqrt(eps^2 + dot(theta,theta,1));
    omega = 1./(1e-10+dot(Btheta,Btheta,1));
    omega(1,1,:) = N* omega(1,1,:);
    omega = omega./sum(omega,2);
    Bc = sum(Btheta.*omega,2);
    ulim = reshape(Nu,1,1,[]) .* eps .* Bc ./ sqrt(1-dot(Bc,Bc,1));
    ulim(:,:,Nu<1e-15) = 0;
    ulim = reshape(ulim,dim,[]);
end