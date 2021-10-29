
% d2,d1      d2,d1 1    0,d1
% d2,d1 3    d2,d1 c    0,d1 4
% d2,d1      d2,d1 2    0,d1

d2 = 100;
d1 = 1;


uc = [d2,d1]';
u1 = [d2,d1]';
u2 = [d2,d1]';
u3 = [d2,d1]';
u4 = [1,d1]';

Rot = rotz(56);
Rot = Rot(1:2,1:2);

epsC = 1;
N = 1;

P1=[limit_t(uc.*[0;1],[u1,u2,u3,u4].*[0;1],epsC,N),...
    limit_t(uc.*[1;0],[u1,u2,u3,u4].*[1;0],epsC,N),...
    limit_t(uc,[u1,u2,u3,u4],epsC,N)];

P1c = [[P1(1,2);P1(2,1)],P1(:,3)]

Ru = Rot* uc;
Rus = Rot* [u1,u2,u3,u4];

P2 = [limit_t(Ru.*[0;1],Rus.*[0;1],epsC,N),...
    limit_t(Ru.*[1;0],Rus.*[1;0],epsC,N),...
    limit_t(Ru,Rus,epsC,N)];

P2c = Rot\[[P2(1,2);P2(2,1)],P2(:,3)]





function ulim = limit_t(u,uS,eps,N)
    Nu = norm(u);
    if(Nu < 1e-10)
       ulim = u * 0;
       return;
    end
    theta = [u,uS]/Nu;
    Btheta = theta./sqrt(eps^2 + dot(theta,theta,1));
    omega = 1./(1e-10+dot(Btheta,Btheta,1));
    omega(1) = N* omega(1);
    omega = omega/sum(omega);
    Bc = sum(Btheta.*omega,2);
    ulim = Nu * eps * Bc / sqrt(1-dot(Bc,Bc,1));
end