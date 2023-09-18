function k_B2nz=function_SEM_homogenius_TBeam(w,rho,A,E,Iz,Kap,L,G)
% builds dynamic stiffness matriz for a thimoshenko beam spectral
% element (homogenius), for a given frequency (corresponding to k_local(w).
%------------------------------------------------------------------------
% 2-noded beam element Dynamic Spectral Stiffness Matrix calculation.
%
%          **********   TIMOSHENKO BEAM   **********
% 
%    k_B2nz  --> element dynamic matrix for 2-noded beam (z-axis)
%    w       --> frequency [rad]
%    rho     --> member density
%    A       --> member x-sectional area
%    E       --> member Young's modulus
%    Iz      --> moment of inertia of member about z-axis
%    Kap     --> x-sectional constant (Kap=5/6 for rectangular x-section)
%    L       --> element length
%    G       --> Shear modulus
%
 

im=sqrt(-1);

Cs=sqrt(G*A*Kap/rho/A);

Q=sqrt(rho*Iz/rho/A);
c0q=sqrt(E*Iz/rho/A);

k1=+sqrt(0.5*((1/Cs)^2+(Q/c0q)^2)*w^2+sqrt((w/c0q)^2+0.25*((1/Cs)^2-(Q/c0q)^2)^2*w^4));
k2=-sqrt(0.5*((1/Cs)^2+(Q/c0q)^2)*w^2-sqrt((w/c0q)^2+0.25*((1/Cs)^2-(Q/c0q)^2)^2*w^4));

xsi1=k1*L;
xsi2=k2*L;

z11=1-exp(-im*xsi1)*exp(-im*xsi2);
z12=exp(-im*xsi1)-exp(-im*xsi2);
z21=exp(-im*xsi1)+exp(-im*xsi2);
z22=1+exp(-im*xsi1)*exp(-im*xsi2);

P1=im/(G*A*Kap*k1)*(E*Iz*k1^2+G*A*Kap-rho*Iz*w^2);
P2=im/(G*A*Kap*k2)*(E*Iz*k2^2+G*A*Kap-rho*Iz*w^2);

r1=(P1-P2)*z11;
r2=(P1+P2)*z12;

dett=r1^2-r2^2;

k11=(xsi2^2-xsi1^2)*(r1*z22+r2*z21)*L/dett;
k12=(-im*xsi2*(r1*z11+r2*z12)+im*xsi1*(r1*z11-r2*z12))*L^2/dett;
k13=(xsi1^2-xsi2^2)*(r1*z21+r2*z22)*L/dett;
k14=(-im*xsi1*(r1*z12-r2*z11)-im*xsi2*(r1*z12+r2*z11))*L^2/dett;

k21=k12;
k22=(-im*xsi1*P2+im*xsi2*P1)*(r1*z22-r2*z21)*L^2/dett;
k23=-k14;
k24=(im*xsi1*P2-im*xsi2*P1)*(r1*z21-r2*z22)*L^2/dett;

k31=k13;
k32=k23;
k33=k11;
k34=-k12;

k41=k14;
k42=k24;
k43=k34;
k44=k22;


k_B2nz=E*Iz/L^3*[ k11  k12  k13  k14
                  k21  k22  k23  k24
                  k31  k32  k33  k34
                  k41  k42  k43  k44 ];
end
                  