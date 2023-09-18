function f=testfsolve(w)
im=-1i;
% p6=1;
% p5=0;
p4=4;
p3=0;
p2=1;
p1=w;
p0=2;
gamma = roots([ p4 p3 p2 p1 p0]);
%% analytical solution for third degree polynomial on \gamma^2
% a1=p4/p6;
% a2=p2/p6;
% a3=p0/p6;
% %eq_caracteristica = [p6 p5 p4 p3 p2 p1 p0];
% %
% % Raízes da equação característica
% %rr = roots(eq_caracteristica);
% Q= (3*a2-a1^2)/9;
% R=(9*a1*a2 -27*a3-2*a1^3)/54;
% S1=(R+(Q^3+R^2)^(1/2))^(1/3);
% S2=(R-(Q^3+R^2)^(1/2))^(1/3);
% r1=S1+S2-a1/3;
% r2=-(S1+S2)/2-a1/3+im*sqrt(3)*(S1-S2)/2;
% r3=-(S1+S2)/2-a1/3-im*sqrt(3)*(S1-S2)/2;
% r1=(-p1+sqrt(p1^2-4*p2*p0))/2/p2;
% r2=(-p1+sqrt(p1^2+4*p2*p0))/2/p2;
% 
% gamma1= sqrt(r1);gamma2= -gamma1; 
% gamma3= sqrt(r2);gamma4= -gamma3;
% %gamma5= sqrt(r3);gamma6= -gamma5;
% gamma = [gamma1;gamma2;gamma3;gamma4];
f = gamma(1)^2 + w;
end