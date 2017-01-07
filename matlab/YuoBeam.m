syms  x L  E Iy Iz k G A L u1 v1 w1 u2 v2 w2 f1 t1 ...
    p1 f2 t2 p2 u3 v3 w3 f3 t3 p3 u(x) v(x) w(x) f(x) t(x) p(x) real

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

code = ccode(K);
