%% clear 
clear all 
%close all 
clc 
%% lumped-parameter system set up
k=1; m1=1; m2=m1;
%k2=1; m2=2;
%r1=0; r2=0; %range of non-reciprocal coupling, r=0 means next-neighborhoud
%couple
r=0;
eta=0.01;
e='m';
% boundary == 0 : free-free (default)
% bundary == 1: fixed-fixed
% boundary == 2 : PBC
boundary=2;
% angular frequency vector limit
flim=3.5e3;
wv=2*pi*(2:2:flim);
fv=wv/2/pi;
Nd=100*2*pi*flim;
% gains
kp=-1;
ki=kp;
kd=kp;
gamma_c=[1*kp 0*ki 0*kd];
s=tf('s');
% qusiperiodicity parameters
theta=0;phi=0;
%% FEM
ncell=10;
output = function_buildFEM_xl(m1,m2,k, eta ,ncell,flim,'boundary',boundary);
%M=output.M;C=output.C;K=output.K;
Ms=output.Ms;Cs=output.Cs;Ks=output.Ks;
x=output.x;ndof=output.ndof;%% State Space of interconnection (S,H) 
%% S
output=function_buildSS_xl(Ms,Cs,Ks,ndof,e,r);
Ass=output.Ass;nx=size(Ass,1);
% PS: (Ms,Cs,Ks) were previously multiplied by B*rho,
% so we need to multiply matrixes F (inside B1ss) and T(inside B2ss) also 
B1ss=output.B1ss;nw=size(B1ss,2);
B2ss=output.B2ss;nu=size(B2ss,2);
C1ss=output.C1ss;nz=size(C1ss,1);
C2ss=output.C2ss;ny=size(C2ss,1);
D11ss=output.D11ss;
D12ss=output.D12ss;
D21ss=output.D21ss;
%% Internal stalibility test
%  output.feas=0 = > internally unstable
% output.feas=1 = > internally stable
 %output = function_stability_c(Ass)
B2ss_contr=ones(nx,1);D12ss_contr=zeros(nx,1);C1ss_contr=eye(nx);
% output = function_hinf_norm_ct(Ass,B1ss,C1ss,D11ss)
output = function_Hinf_Optimal_Controller_c(Ass,B1ss,B2ss_contr,C1ss_contr,D12ss_contr,D12ss_contr)
 %%
 D22ss=zeros(ny,nu);
sys=ss(Ass,B2ss,C2ss,D22ss);
ts = 0.0004;
sys_d=c2d(sys,ts);
[Ad,Bd,Cd,Dd]=ssdata(sys_d);
% % note that time discretization does not change matrixes B2,D12,C2 and D21
output=function_H2_Dynamic_Controller_dt(Ad,B1ss,Bd,C1ss,Cd,D11ss,D12ss,D21ss)