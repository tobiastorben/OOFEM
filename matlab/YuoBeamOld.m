syms N1 N2 r Hv1 Hv2 Hw1 Hw2 Ht1 Ht2 Hp1 Hp2 ...
    Gv1 Gv2 Gw1 Gw2 Gt1 Gt2 Gp1 Gp2 x L ay az ...
    by bz E Iy Iz k G A L u1 v1 w1 u2 v2 w2 f1 t1 ...
    p1 f2 t2 p2 u3 v3 w3 f3 t3 p3 u(x) v(x) w(x) f(x) t(x) p(x) real

r = x/L;
ay = 12*E*Iy/(k*G*A*L^2);
by = 1/(1-ay);
az = 12*E*Iz/(k*G*A*L^2);
bz = 1/(1-az);

% 
% N1 = 1-r;
% N2 = r;
% Hv1 = by*(2*r^3 -3*r^2 + ay*r + 1 -ay);
% Hv2 = by*(-2*r^3 + 3*r^2 - ay*r);
% Hw1 = bz*(2*r^3 -3*r^2 + az*r + 1 -az);
% Hw2 = bz*(-2*r^3 + 3*r^2 - az*r);
% Ht1 = L*by*(r^3 + ((1/2)*ay-2)*r^2 + (1- (1/2)*ay)*r);
% Ht2 = L*by*(r^3 -(1+ (1/2)*ay)*r^2 + (ay/2)*r);
% Hp1 = L*bz*(r^3 + ((1/2)*az-2)*r^2 + (1- (1/2)*az)*r);
% Hp2 = L*bz*(r^3 -(1+ (1/2)*az)*r^2 + (az/2)*r);
% Gv1 = (6*by/L)*(r^2-r);
% Gv2 = (6*by/L)*(-r^2+r);
% Gw1 = (6*bz/L)*(r^2-r);
% Gw2 = (6*bz/L)*(-r^2+r);
% Gt1 = by*(3*r^2 + (ay-4)*r + 1 - ay);
% Gt2 = by*(3*r^2 -(ay+2)*r);
% Gp1 = bz*(3*r^2 + (az-4)*r + 1 - az);
% Gp2 = bz*(3*r^2 -(az+2)*r);

uvec = [u1, v1, w1, f1, p1, t1, u2, v2, w2, f2, p2, t2]';
 
  eq1 = E*A*diff(diff(u)) == 0;
  eq2 = k*G*A*(diff(diff(v))-diff(t))==0;
  eq3 = k*G*A*(diff(diff(w))+diff(p))==0;
  eq4 = E*Iy*diff(diff(t))-k*G*A*(diff(v)-t) == 0;
  eq5 = E*Iz*diff(diff(p))+k*G*A*(diff(w)+p) == 0;
  eq6  = k*G*(Iy+Iz)*diff(diff(f)) == 0;
  
  c1 = u(0) == u1;
  c2 = v(0) == v1;
  c3 = w(0) == w1;
  c4 = f(0) == f1;
  c5 = p(0) == p1;
  c6 = t(0) == t1;
  c7 = u(L) == u2;
  c8 = v(L) == v2;
  c9 = w(L) == w2;
  c10 = f(L) == f2;
  c11 = p(L) == p2;
  c12 = t(L) == t2;
S = dsolve(eq1,eq2,eq3,eq4,eq5,eq6,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12);
N = equationsToMatrix([S.u S.v S.w S.f S.p S.t]' == zeros(6,1), uvec);
 

% N = [N1 0 0 0 0 0 N2 0 0 0 0 0;
%      0  Hv1 0 0 Ht1 0 0 Hv2 0 0 Ht2 0;
%      0  0 Hw1 0 0 Hp1 0 0 Hw2 0 0 Hp2;
%      0  0 0 N1 0 0 0 0 0 N2 0 0;
%      0  Gv1 0 0 Gt1 0 0 Gv2 0 0 Gt2 0;
%      0  0 Gw1 0 0 Gp1 0 0 Gw2 0 0 Gp2];
     

B =[diff(N(1,:),x);
    diff(N(2,:),x)-N(6,:);
    diff(N(3,:),x)+N(5,:);
    diff(N(4,:),x);
   -diff(N(5,:),x);
    diff(N(6,:),x)];

D =[ E*A 0 0 0 0 0;
     0 k*G*A 0 0 0 0;
     0 0 k*G*A 0 0 0;
     0 0 0 k*G*(Iy+Iz) 0 0;
     0 0 0 0 E*Iz 0;
     0 0 0 0 0 E*Iy];
 
 integrand = B.'*D*B;
 
 K = int(integrand,x,0,L);
 
 
 
 
%  M =[0 0 0 0 0 0 1 0 0 0 0 0;
%      1 0 0 0 0 0 0 0 0 0 0 0;
%      0 0 0 0 0 0 0 0 1 0 0 0;
%      0 0 0 0 1 0 0 0 0 0 0 0;
%      0 1 0 0 0 0 0 0 0 0 0 0;
%      0 0 0 0 0 0 0 0 0 -1 0 0;
%      0 0 0 0 0 0 1 L 0 0 0 0;
%      1 L L^2 L^3 0 0 0 0 0 0 0 0;
%      0 0 0 0 0 0 0 0 1 L L^2 L^3;
%      0 0 0 0 1 L 0 0 0 0 0 0;
%      0 1 2*L 3*L^2 0 0 0 0 0 0 0 0;
%      0 0 0 0 0 0 0 0 0 -1 -2*L -3*L^2;];
%  
%     dy = (6*E*Iy)/(k*G*A);
%     dz = (6*E*Iz)/(k*G*A);
% 
%  y = [u1, v1, w1, f1, t1+dy, p1-dz, u2, v2, w2, f2, t2+dy, p2-dz]';
%  
%  c = M\y;
%  
%  Nutest =[c(7) + c(8)*x;
%      c(1) + c(2)*x + c(3)*x^2 + c(4)*x^3;
%      c(9) + c(10)*x + c(11)*x^2 + c(12)*x^3;
%      c(5) + c(6)*x;
%      c(2)-dy + 2*c(3)*x + 3*c(4)*x^2;
%      -c(10) + dz - 2*c(11)*x - 3*c(12)*x^2];
%  
%  [matrixform,bla] = equationsToMatrix(Nutest == [u3 v3 w3 f3 t3 p3]',[u1, v1, w1, f1, t1, p1, u2, v2, w2, f2, t2, p2]);
 
  