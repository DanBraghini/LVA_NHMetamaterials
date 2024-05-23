% Created on : 27/09/2021 by Danilo Braghini
% Bulk-Boundary correspondence or  metamaterial(PBC)-metastructure(OBC) correspondence
%% clear 
clear all 
%close all 
clc 
%% Macros
global dmodel rho L1 L2 Lc A1 A2 c eta gamma_c H_pv kp ki kd ideal_filter kL T1 T2 Nd
%feedback = input('enter the type of feddback\n(string: 
feedback='i';%PC properties
kp=1e-2;
%ki=1.5e-2;
ki=-1.5e-3;
kd=2e-10;
%dmodel = input('enter the selected damping model\n(v = viscous, s = structural, n= none)')
dmodel='s';
ideal_filter=0;
% Laplace variable
s=tf('s');
% angular frequency vector
flim=1.5e3;%in kHz
% parameter used for time derivative aproximation
Nd=100*2*pi*flim;
%Nd=100*2*pi*1.5;
wv=2*pi*(2:1:flim);
%wv=2*pi*(-flim:2:flim);
fv=wv/2/pi;
% pb =input('enter the number of pass bands / bulk bands\n(any integer number between 1 and 10)')
pb=4;
% a=coupling range (in cells)
a=0;
% filter = 0 : no filter
% filter = 1 : butterworth filter
% filter = 2 : RC filter
% filter = 3 : filter built for the experiment
filter=0;
plotfilter=0;
% Constant transducer models:
% T1 = sensor, T2 = actuator
%T1=1e-2; T2=T1;
T1=1;T2=1;
function_metasetup(filter,plotfilter,feedback);
% labelplot(1): wavenumber x real frequency
% labelplot(2): wavenumber x imaginary frequency
% labelplot(3): complex frequency plane
% labelplot(4): 3D plot
labelplot=[0 0 1 0];
% periodic=0;
% passive=0;
% colors=0;
%Boundary conditions of the finite structure:
%            boundary == 0 : free-free (default)
%            boundary == 1 : fixed-fixed
%            boundary == 2 : periodic(infinity system)
boundary=2;
% number of cells on the finite structure
ncell=50;
%FEM
ne_cell=7*3;
xs = Lc/4;
xi = xs+Lc/2;
% excitation point for zero-state response
e='m';

%% initial conditions set up for w(k) method
m=function_initialconditions_setup(pb,fv);
%% Inverse method via transfer matrix and SEM formulation
kLv=-pb*pi:8*pi/400:pb*pi;
N=length(kLv); 
fun= @caracteristic_equation_T_Acoustic;
options = optimoptions('fsolve','Display','off');
options.StepTolerance =1e-21;
options.OptimalityTolerance = 1e-21;
options.MaxIterations = 100;
w_TM=zeros(pb,N);w_TM_passive=w_TM;

% This loop runs through the vector kLv.
% PS: remember to adjust the initial conditions depending on the PC. You
% will need to adjust wv range also via flim.
        
gamma_c_aux =gamma_c;
for n =1:N
    kL=kLv(n);
    % for non-local feedback, introduce spatial delay of a cells
    gamma_c = gamma_c_aux.*exp(1i*a*kL);
    for i=1:pb
        w_TM(i,n) = fsolve(fun,wv(m(i)),options);
    end

end
gamma_c = gamma_c_aux;
% Detach frequency and damping factor for each bulk band and converting to kHz:
for i=1:pb
    w_real_TM{i} = real(w_TM(i,:))/2/pi/1000;
    w_imag_TM{i} = imag(w_TM(i,:))/2/pi/1000;
end

%% FEM
output = function_buildFEM(rho, L1, L2, Lc, A1, A2, c, eta,ne_cell,ncell,flim,xs,xi...
    ,'dmodel',dmodel,'boundary',boundary);
M=output.M;C=output.C;K=output.K;
Ms=output.Ms;Cs=output.Cs;Ks=output.Ks;
x=output.x;ndof=output.ndof;
n1=output.n1;n2=output.n2;
%% State Space of interconnection (S,H) 
%% S
output=function_buildSS(Ms,Cs,Ks,ndof,ncell,n1,n2,e,a);
Ass=output.Ass;nx=size(Ass,1);
% PS: (Ms,Cs,Ks) were previously multiplied by B*rho,
% so we need to multiply matrixes F (inside B1ss) and T(inside B2ss) also 
B=rho*c^2;
B1ss=(B*rho).*output.B1ss;nw=size(B1ss,2);
B2ss=(B*rho).*output.B2ss;nu=size(B2ss,2);
C1ss=output.C1ss;nz=size(C1ss,1);
C2ss=output.C2ss;ny=size(C2ss,1);
D11ss=output.D11ss;
D12ss=output.D12ss;
D21ss=output.D21ss;
%% H 
% Feedback law in pressure by volume velocity
H_v=gamma_c*H_pv;
% Correction to volume acceleration
H=H_v*Nd/(1+Nd/s);

%%  closed loop (S,H)
if  norm(gamma_c,1)~=0
%     figure
%     bode(H)
%     title('bode diagram of the feedback law')
    [num,dem]=tfdata(H,'v');
    if size(num,1) > size(dem,1)
        disp('improper feedback law\n');
    end
    %% controlable state space realization of H
    [Ac1,Bc1,Cc1,Dc1]=tf2ss(num,dem);
    nxc=size(Ac1,1);
    Acss=zeros((ncell-a)*nxc,(ncell-a)*nxc);
    Bcss=zeros((ncell-a)*nxc,ncell-a);
    Ccss=zeros(ncell-a,(ncell-a)*nxc);
    %% building state space realization of H (Acss,Bcss,Ccss,Dcss)
    % by blocks, since each cell has a decoupled feedback law
    for k=1:ncell-a
        lines=(k-1)*nxc+1:k*nxc;
        Acss(lines,lines)=Ac1;
        Bcss(lines,k)=Bc1;
        Ccss(k,lines)=Cc1;
    end
    %Acss=Ac1*eye(ncell);   Bcss=Bc1*eye(ncell);
    %Ccss=Cc1*eye(ncell);  
    Dcss=Dc1*eye(ncell-a);
    %% closed loop (S,H) interconnection
    Ass_cl=[Ass+B2ss*Dcss*C2ss B2ss*Ccss
            Bcss*C2ss           Acss];
    Bss_cl=[B1ss+B2ss*Dcss*D21ss
                Bcss*D21ss];
    Css_cl=[C1ss+D12ss*Dcss*C2ss     D12ss*Ccss];
    Dss_cl=D11ss+D12ss*Dcss*D21ss;
    sys_wz = ss(Ass_cl,Bss_cl,Css_cl,Dss_cl);
    Css=Css_cl;
    nx_cl=size(Ass_cl,1);
else 
    % only logical case other then the previous is gamma_c==0, meaning open
    % loop system, u =0
    sys_wz = ss(Ass,B1ss,C1ss,D11ss);
    Css=C1ss;
end
%% verify internal stability of closed loop system
if  norm(gamma_c,1)~=0
    [V, Lambda] = eig(Ass_cl);
else
    [V, Lambda] = eig(Ass);
end
Lambda = diag(Lambda);
Vn = Css*V;
wn = -1i*Lambda;
[wn,ind] = sort(wn,'ComparisonMethod','real');
Vn = Vn(:,ind); 
ind_zero = find(wn == 0);
% getting away negative eigenfrequencies and their modes
if ~isempty(ind_zero)
    wn(1:ind_zero(1)) = [];
    Vn(:,1:ind_zero(1)) = [];
end
%%
%-------------------------------------------------------------------------%
%          Unfolded Dispersion Diagram Plots by SEM        %
%-------------------------------------------------------------------------%
% labelplot=[1 1 0 1];
% function_unfolded_dispersion_diagrams(pb,kLv,w_real_TM,w_imag_TM,w_real_TM,w_imag_TM...
%     ,'setplot',labelplot,'periodic',1);
% zlim([-1 1])
%%
%-------------------------------------------------------------------------%
%                   Bulk-Boundary Correspondence                          %
%-------------------------------------------------------------------------%
function_unfolded_dispersion_diagrams(pb,kLv,w_real_TM,w_imag_TM,w_real_TM,w_imag_TM...
   ,'setplot',labelplot);
zlim([-1 1])
hold on
scatter(imag(wn/2/pi/1000),real(wn/2/pi/1000),'k')
 ylim([0 1.35])
 % xlim([-0.3 0.3]);
  xlim([-2e-2 2e-2]);
%% Eigenmodes Animation
% eigenfrequency search
% fnv=wn/2/pi/1000;
% [fi_imag,fi_real]=ginput(1)
% dfi=1/2;
% j=1;
% % searching the eigenfrequency by its real and imaginary parts
% ind_real=find(abs(real(fnv)-fi_real)<=dfi);
% ind_imag=find(abs(imag(fnv)-fi_imag)<=dfi);
% ind=intersect(ind_real,ind_imag);
% dfi=dfi/10;
% while length(ind) ~= 1 && j <=100
%     ind_real=find(abs(real(fnv)-fi_real)<=dfi);
%     ind_imag=find(abs(imag(fnv)-fi_imag)<=dfi);
%     ind=intersect(ind_real,ind_imag);
%     % updates accuracy until the only eigenfrequency find is the selected
%     if isempty(ind)
%         dfi=2*dfi;
%     else
%         dfi=dfi/10;
%     end
%     
%     j=j+1;
% end
% % norm1 = max(abs(real(Vn(:,mode1))));
% disp('ind [k Hz]')
% fn=fnv(ind);
% display(fn)
% %modo=sin(x)+1j*cos(x)/2;
% mode=Vn(:,ind);
% norm = max(abs(mode));
% % numero de frames
% Nf=100;
% % numero de ciclos
% nc=5; 
% teta=0:2*nc*pi/(Nf-1):2*nc*pi; 
% figure
% % l=length(x);
% % x = [x x(l).*ones(1,l-1)+x(2:l)];
% for i=1:Nf
%     y=abs(mode).*cos(angle(mode)+teta(i));
%     plot(x(1:ndof),y'/norm)
% %    y=[y;y];
% %    plot(x(1:end-1),y'/norm)
%     axis([min(x),max(x), -1,1]);
%     xlabel('$x$ [m]', 'interpreter', 'latex', 'fontsize', 15)
%     ylabel('$ \hat{P}(x,\omega_n)/ ||\hat{P}(x,\omega_n)||_{\infty}$', 'interpreter', 'latex', 'fontsize', 15)     
%     set(gca,'TickLabelInterpreter','Latex','fontsize',15);
%     box on
%     set(gcf, 'Color', 'w');
%     legend(['\omega_n = ' num2str(fn,'%0.3f') 'kHz']);
%     drawnow
% end