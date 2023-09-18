%% clear 
clear all 
%close all 
clc 
%%
% Created on : 17/05/2021 by Danilo Braghini
flim=3.5e3;
wv=2*pi*(2:2:flim);
fv=wv/2/pi;
global feedback rho L1 L2 Lc A1 A2 gamma_c c eta kL Nd
%feedback = input('enter the type of feddback\n(string: d = derivative, i = integrative and p = proportional)\n','s');
feedback='i';%PC properties
s=tf('s');
%%
%PC properties
% acoustic metamaterial
eta=0.035;
Lc = 50e-2;
r1 = 2e-2;
A1 = pi*r1^2;
r2 = r1;
A2 = A1;
%xs = Lc/4;
%xi = xs+Lc/2;
xs = Lc/4;
xi = xs+Lc/2;
rho =1.225;
c=343;
L1= xs;
L2 = xi- xs;
% Feedback control law and sensor locality paramater 
if feedback=='i'
    gamma_c=1.5e-3;
   % gamma_c=0;
elseif feedback=='d'
    gamma_c=2e-10;
elseif feedback=='p'
    gamma_c=2e-14;
elseif feedback=='lead-lag'
    gamma_c=1.5e-3;
end
%gamma_c = -.01*rho*A/Lc;
a=0;
%% Dispersion Relation
%
%-------------------------------------------------------------------------%
%                               SEM   Model                               %  
%-------------------------------------------------------------------------%
% use this section to define the initial values of fsolve optimization

[kL_sem_PB,kL_sem_SB] = function_SEM_Acoustic(L1,L2,rho,A1,A2,c,wv,eta,'gain',gamma_c,'feedback',feedback,'Nd',Nd);
figure
plot(real(fv)/1000,abs(kL_sem_PB)/pi,'k--',real(fv)/1000,-abs(kL_sem_SB)/pi,'r--');grid on
hold on

[kL_sem_PB_passive,kL_sem_SB_passive] = function_SEM_Acoustic(L1,L2,rho,A1,A2,c,wv,eta,'feedback',feedback);
plot(real(fv)/1000,abs(kL_sem_PB_passive)/pi,'k',real(fv)/1000,-abs(kL_sem_SB_passive)/pi,'r');grid on
xlabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
ylabel('$k L_c/ \pi$', 'interpreter', 'latex', 'fontsize', 15)


%% initial conditions set up for w(k) method
%output = function_find_mean_frequncies(kL_sem_PB,fv);
% if output.nPBs < 4
%     pb = output.nPBs;
% else
%     pb = 4;
% end
pb=10;
df=fv(2)-fv(1);
if feedback=='d'
    f1=find(abs(fv-195.7) <= df/2);f2=find(abs(fv-540.7) <= df/2);
    m=zeros(1,4);
    m(1)=f1; m(2)=f2;
    f3=find(abs(fv-880.3) <= df/2);f4=find(abs(fv-1184) <= df/2);
    %
    f5=find(abs(fv-1580) <= df/2);f6=find(abs(fv-1860) <= df/2);
     f7=find(abs(fv-2190) <= df/2);f8=find(abs(fv-2620) <= df/2);
       f9=find(abs(fv-2950) <= df/2);f10=find(abs(fv-3200) <= df/2);
    m(3)=f3;m(4)=f4;
    %
    m(5)=f5;m(6)=f6;m(7)=f7;m(8)=f8;m(9)=f9;m(10)=f10;
elseif feedback=='i'
    
    f1=find(abs(fv-232) <= df/2);f2=find(abs(fv-488) <= df/2);
    m=zeros(1,4);
    m(1)=f1; m(2)=f2;
    f3=find(abs(fv-886) <= df/2);f4=find(abs(fv-1168) <= df/2);
    %
    f5=find(abs(fv-1576) <= df/2);f6=find(abs(fv-1856) <= df/2);
     f7=find(abs(fv-2232) <= df/2);f8=find(abs(fv-2574) <= df/2);
       f9=find(abs(fv-2918) <= df/2);f10=find(abs(fv-3258) <= df/2);
    m(3)=f3;m(4)=f4;
    %
    m(5)=f5;m(6)=f6;m(7)=f7;m(8)=f8;m(9)=f9;m(10)=f10;
    
end
m_passive=m;
%% Inverse method via transfer matriz and SEM formulation
kLv=-pb*pi:8*pi/400:pb*pi;
N=length(kLv); 
fun= @caracteristic_equation_T_Acoustic;
options = optimoptions('fsolve','Display','off');
options.StepTolerance =1e-20;
options.OptimalityTolerance = 1e-20;
options.MaxIterations = 100;
w_TM=zeros(pb,N);w_TM_passive=w_TM;

% This loop runs through the vector kLv.
% PS: remember to adjust the initial conditions depending on the PC. You
% will need to adjust wv range also.
        
gamma_c_aux =gamma_c;
for n =1:N
    kL=kLv(n);
    gamma_c = gamma_c_aux.*exp(1i*a*kL);
    for i=1:pb
        w_TM(i,n) = fsolve(fun,wv(m(i)),options);
    end

end
gamma_c = gamma_c_aux;
%passive system
gamma_c_aux =gamma_c;
for n =1:N
    kL=kLv(n);
    gamma_c = 0;
    for i=1:pb
        w_TM_passive(i,n) = fsolve(fun,wv(m_passive(i)),options);
    end

end
gamma_c = gamma_c_aux;

% Detach frequency and damping factor for each bulk band:
for i=1:pb
    f_real_TM{i} = real(w_TM(i,:))/2/pi/1000;
    f_imag_TM{i} = imag(w_TM(i,:))/2/pi/1000;
    f_real_TM_passive{i} = real(w_TM_passive(i,:))/2/pi/1000;
    f_imag_TM_passive{i} = imag(w_TM_passive(i,:))/2/pi/1000;
end
%%
% Indexing the boundaries between the Brillouin Zones
dk = kLv(2)-kLv(1);
kind=zeros(1,2*pb);
j=1;
for i=1:2:2*pb
   kind(i)=find(abs(kLv + j*pi) <= dk/2);
   kind(i+1)=find(abs(kLv - j*pi) <= dk/2);
   j=j+1;
end
%% FEM
% ncell:number of cells that make the structure (multiple of 3 to be divided between the 3 segments)
% number of elements per cell
ne_cell=45;
%number of cells on the structure
ncell=18;
output = function_buildFEM(ne_cell,ncell,flim,xs,xi,'boundary',1);

M=output.M;C=output.C;K=output.K;
Ms=output.Ms;Cs=output.Cs;Ks=output.Ks;
x=output.x;ndof=output.ndof;
n1=output.n1;n2=output.n2;
%% State Space of (S,H) interconnection 
% S
output=function_buildSS(Ms,Cs,Ks,ndof,ncell,n1,n2);
Ass=output.Ass;
B1ss=output.B1ss;
B2ss=output.B2ss;
C1ss=output.C1ss;
C2ss=output.C2ss;
D11ss=output.D11ss;
D12ss=output.D12ss;
D21ss=output.D21ss;
% H 
if feedback=='i'
     %H=gamma_c*1/s;
     H=7020/(0.1551*s^3+1600.5*s^2+495000*s+3125000);
elseif feedback=='d'
    Nd=2*pi*flim;
    %derivative aproximation
    H=gamma_c*Nd/(1+Nd/s);
elseif feedback=='p'
    %output static compensator L
    L = gamma_c.*eye(ncell);
elseif feedback=='lead-lag'
    H=gamma_c*(s-Nd/10)/(s-Nd);
end
if feedback~='p'
%     figure
%     bode(H)
    [num,dem]=tfdata(H,'v');
    [Ac1,Bc1,Cc1,Dc1]=tf2ss(num,dem);
    nx=size(Ac1,1);
    Acss=zeros(ncell*nx,ncell*nx);
    Bcss=zeros(ncell*nx,ncell);
    Ccss=zeros(ncell,ncell*nx);
    % building state space representation of H (Acss,Bcss,Ccss,Dcss)
    % by blocks, since each cell has a decoupled controller
    for k=1:ncell
        lines=(k-1)*nx+1:k*nx;
        Acss(lines,lines)=Ac1;
        Bcss(lines,k)=Bc1;
        Ccss(k,lines)=Cc1;
    end
    %Acss=Ac1*eye(ncell);   Bcss=Bc1*eye(ncell);
    %Ccss=Cc1*eye(ncell);  
    Dcss=Dc1*eye(ncell);
    %closed loop (S,H) interconnection
    Ass_cl=[Ass+B2ss*Dcss*C2ss B2ss*Ccss
            Bcss*C2ss           Acss];
    Bss_cl=[B1ss+B2ss*Dcss*D21ss
                Bcss*D21ss];
    Css_cl=[C1ss+D12ss*Dcss*C2ss     D12ss*Ccss];
    Dss_cl=D11ss+D12ss*Dcss*D21ss;
else
    %closed loop (S,L) interconnection
    Ass_cl=Ass+B2ss*L*C2ss;
    Bss_cl=B1ss+B2ss*L*D21ss;
    Css_cl=C1ss+D12ss*L*C2ss;
    Dss_cl=D11ss+D12ss*L*D21ss;
end
sys_wz = ss(Ass_cl,Bss_cl,Css_cl,Dss_cl);

%%
[V, Lambda] = eig(Ass_cl);
Lambda = diag(Lambda);
%Vn = Css_cl*V;
 wn = -1i*Lambda;
 [wn,ind] = sort(wn,'ComparisonMethod','real');
% Vn = Vn(:,ind); 
%%
%-------------------------------------------------------------------------%
%          Unfolded Dispersion Diagram Plots by SEM        %
%-------------------------------------------------------------------------%
% 
function_unfolded_dispersion_diagrams(pb,kind,kLv,f_real_TM,f_imag_TM,f_real_TM_passive,f_imag_TM_passive);


%% complex plane
figure
    
plot(f_real_TM{1}, f_imag_TM{1},'k')
hold on
plot(f_real_TM{2}, f_imag_TM{2},'k')
plot(f_real_TM{3}, f_imag_TM{3},'k')
plot(f_real_TM{4}, f_imag_TM{4},'k')
plot(f_real_TM{5}, f_imag_TM{5},'k')
plot(f_real_TM{6}, f_imag_TM{6},'k')
plot(f_real_TM{7}, f_imag_TM{7},'k')
plot(f_real_TM{8}, f_imag_TM{8},'k')
plot(f_real_TM{9}, f_imag_TM{9},'k')
plot(f_real_TM{10}, f_imag_TM{10},'k')
plot(f_real_TM_passive{1}, f_imag_TM_passive{1},'r')
plot(f_real_TM_passive{2}, f_imag_TM_passive{2},'r')
plot(f_real_TM_passive{3}, f_imag_TM_passive{3},'r')
plot(f_real_TM_passive{4}, f_imag_TM_passive{4},'r')
plot(f_real_TM_passive{5}, f_imag_TM_passive{5},'r')
plot(f_real_TM_passive{6}, f_imag_TM_passive{6},'r')
plot(f_real_TM_passive{7}, f_imag_TM_passive{7},'r')
plot(f_real_TM_passive{8}, f_imag_TM_passive{8},'r')
plot(f_real_TM_passive{9}, f_imag_TM_passive{9},'r')
plot(f_real_TM_passive{10}, f_imag_TM_passive{10},'r')
xlabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
ylabel('$\Im(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
zlabel('$k L_c/ \pi$', 'interpreter', 'latex', 'fontsize', 15)

scatter(real(wn)/2/pi/1000,imag(wn)/2/pi/1000,'k')

l_max=max(real(Lambda))
