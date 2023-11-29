% Created on : 27/09/2021 by Danilo Braghini
%% clear 
clear all 
%close all 
clc 
%% Macros
global dmodel rho L1 L2 Lc A1 A2 c eta gamma_c H_pa kp ki kii kd kdd ideal_filter kL Nd
%feedback = input('enter the type of feddback\n(string): 
feedback=1;%PC properties
if feedback==1
    kp=1.5e-2;ki=0;kd=0;kdd=0;kii=0;
elseif feedback==2
     kp=0;ki=1.5;kd=0;kdd=0;kii=0;
elseif feedback==3
     kp=0;ki=0;kd=1e-5;kdd=0;kii=0;
elseif feedback==4
     kp=0;ki=0;kd=0;kdd=5e-10;kii=0;
elseif feedback==5
     kp=0;ki=0;kd=0;kdd=0;kii=1.5;
end
%dmodel = input('enter the selected damping model
%\n(v = viscous, s = structural, n= none)')
dmodel='v';
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
% a=coupling range (# of cells)
a=0;
% filter = 0 : no filter
% filter = 1 : butterworth filter
% filter = 2 : RC filter
% filter = 3 : filter built for the experiment
filter=0;
filterinfo.w1=1000*2*pi;
filterinfo.w2=1500*2*pi;
%lp= 'low pass', bp = 'band pass', 'hp' ='high pass'
filterinfo.type='hp';
plotfilter=0;
function_metasetup(filter,plotfilter,filterinfo);
% labelplot(1): wavenumber x real frequency
% labelplot(2): wavenumber x imaginary frequency
% labelplot(3): complex frequency plane
% labelplot(4): 3D plot
labelplot=[0 0 0 1];
periodic=0;
plotpassive=0;
colors=0;
%Boundary conditions of the finite structure:
%            boundary == 0 : clamped-clamped (default)/ Neumann
%            boundary == 1 : open-open/ Dirichlet
%            boundary == 2 : periodic(infinity system)
boundary=0;
% number of cells on the finite structure
ncell=18;
%FEM
%ne_cell=25*3;
ne_cell=3*3;
xs = Lc/4;
xi = xs+Lc/2;
% excitation point for zero-state response
e='m';
%% choose function
% unfolded plot of the dispersion diagrams w(k) imposing real k (traveling waves)
unfoldedDD=1;
% Bulk-Boundary correspondence or metamaterial(PBC)-metastructure(OBC) correspondence
bulkboundary=1;
% Forced frequency response for metastructure(OBC)
FRF=1;
% Transient response of the metastructure due to tone-burst kind input on
% the middle
transient=0;
% Eigenmodes simulation of the free wave modes
eigenmodes=0;
%% Numerical models
if unfoldedDD==1 || bulkboundary==1 || eigenmodes==1
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
    % w_TM = function_Inverse_Bloch_Dispersion(fun,options, kLv, wv,pb, a);   
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
    %Hermitian counterpart(passive system)
    if plotpassive==1
        gamma_c_aux =gamma_c;
        for n =1:N
            kL=kLv(n);
            gamma_c=[0 0 0];
            for i=1:pb
                w_TM_passive(i,n) = fsolve(fun,wv(m(i)),options);
            end

        end
        gamma_c = gamma_c_aux;
    end
    % Detach frequency and damping factor for each bulk band and converting to kHz:
    for i=1:pb
        w_real_TM{i} = real(w_TM(i,:))/2/pi/1000;
        w_imag_TM{i} = imag(w_TM(i,:))/2/pi/1000;
        if plotpassive==1
            w_real_TM_passive{i} = real(w_TM_passive(i,:))/2/pi/1000;
            w_imag_TM_passive{i} = imag(w_TM_passive(i,:))/2/pi/1000;
        end
    end
end
if bulkboundary == 1 || transient == 1 || eigenmodes == 1
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
    % PS: (Ms,Cs,Ks) were previously multiplied by A*rho on FEM, following 
    % the formulation from Neurodinne&Atalla, section 4.2. The initial 
    % dimension of external inputs to the PDE was [B/A*(volume acceleration)].
    % So we need to multiply matrixes F (inside B1ss) and T(inside B2ss) by
    % (B*rho) so that unitary generalized forces correspond to volume
    % acceleration
    B=rho*c^2;
    B1ss=(B*rho).*output.B1ss;nw=size(B1ss,2);
    B2ss=(B*rho).*output.B2ss;nu=size(B2ss,2);
    C1ss=output.C1ss;nz=size(C1ss,1);
    C2ss=output.C2ss;ny=size(C2ss,1);
    D11ss=output.D11ss;
    D12ss=output.D12ss;
    D21ss=output.D21ss;
    %% H 
    % Feedback law in pressure by volume acceleration
    H=gamma_c*H_pa;

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
        [V, Lambda] = eigs(Ass_cl,nx_cl);
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
end

%% Plots
%%
%-------------------------------------------------------------------------%
%          Unfolded Dispersion Diagram Plots by SEM        %
%-------------------------------------------------------------------------%
if unfoldedDD==1
    clear options;
    options{1}='setplot';
    options{2}=labelplot;
    options{3}='periodic';
    options{4}=periodic;
    options{5}='plotpassive';
    options{6}=plotpassive;
    if plotpassive==1
        options{7}='w_real_passive';options{8}={w_real_TM_passive};
        options{9}='w_imag_passive';options{10}={w_imag_TM_passive};
    end
    function_unfolded_dispersion_diagrams(pb,kLv,w_real_TM,w_imag_TM,options{:});
    %zlim([-1 1])
end
%%
% figure
% for j=1:pb
%   plot(w_real_TM{j}, w_imag_TM{j},'k')
%   hold on
% end
% xlabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
% ylabel('$\Im(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
% set(gca,'TickLabelInterpreter','Latex','fontsize',15);
% box on
% set(gcf, 'Color', 'w');
% grid on

%%
%-------------------------------------------------------------------------%
%                   Bulk-Boundary Correspondence                          %
%-------------------------------------------------------------------------%
%%
if bulkboundary==1 || eigenmodes==1
     clear options;
    options{1}='setplot';
    options{2}=[0 0 1 0];
    options{3}='periodic';
    options{4}=periodic;
    options{5}='plotpassive';
    options{6}=plotpassive;
    function_unfolded_dispersion_diagrams(pb,kLv,w_real_TM,w_imag_TM,options{:});
    zlim([-1 1])
    hold on
 %    figure
    scatter(imag(wn/2/pi/1000),real(wn/2/pi/1000),'k')
    ylim([0 flim/1000])
  %   xlim([-0.3 0.3]);
  % xlim([-2e-2 2e-2]);
%        ylabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
%      xlabel('$\Im(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
%     set(gca,'TickLabelInterpreter','Latex','fontsize',15);
%      box on
%      grid on
%      set(gcf, 'Color', 'w');
end
  %% FRF on pressure by volume velocity
if FRF==1
    ncell=18;
    F=1;
    output=function_FRF_Ac(wv,ncell,e,boundary,F);
    P_SEM=output.ps_v.';
    figure
    plot(fv/1000,20*log10(abs(P_SEM(:,1))),'r','MarkerSize',8,'LineWidth',1)
    hold on
    plot(fv/1000,20*log10(abs(P_SEM(:,end))),'k','MarkerSize',8,'LineWidth',1)
    if plotpassive==1
        P_passive = output.ps_passive.';
        plot(fv/1000,20*log10(abs(P_passive(:,1))),'r--','MarkerSize',8,'LineWidth',1)
        plot(fv/1000,20*log10(abs(P_passive(:,end))),'k--','MarkerSize',8,'LineWidth',1)
    end
    ylabel('$ p / G $ [dB]', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    %legend('L','R','L filtered','R Filtered','passive')
    xlim([0 flim/1000])
    box on
    set(gcf, 'Color', 'w');
    grid on
end
%%
if transient==1
    %% Transient Response: applied force
    % Nc is the number of circles from central frequency fc.
    % T2 define the envelope frequency f2, which has Nc circles within.
    % Nt is the number of entries on time vector t. The heavside function
    % defines a window on half a period of the envelope.
    % time discretization period dt depends on N, in such manner that
    % the total time of the simulation T can be modified if nedeed.
    Nt=4*4096*2;
    % T=0.01; % total analysis time
    Ap=1;
    fc=250;
    %fc=1300;
    wc = 2*pi*fc;
    Nc=15; 
    Tc=1/fc;
    T2=Nc*Tc*2; 
    f2=1/T2; w2=2*pi*f2;
    dt=T2/2/4096;
    %dt = T2/Nt;
    t=0:dt:(Nt-1)*dt;
    w=Ap*sin(wc*t).*sin(w2*t).^2.*heaviside(T2/2-t); 
    %% zero state response of the closed loop for input w
    if norm(gamma_c,1)~=0
        ni=nx_cl;
    else
        ni=nx;
    end
    z=lsim(sys_wz,w,t,zeros(ni,1));
    %%
    Nplots=100;
    tend=length(t)/5;
    %tend=length(t);
%     DI=round(tend/Nplots);
%     figure
%     for i=1:DI:tend
%         plot3(x(1:ndof),(t(i).*10^3)*ones(1,ndof),z(i,:),'b');
%         hold on
%     end
%     xlabel('Length [m]'), ylabel('Time [ms]'), zlabel('Pressure p(x,t)[Pa]')
%     hold off
%     %ylim([0 100])
%     box on
%     set(gcf, 'Color', 'w');
    %%
    figure
    tend=round(tend);
    for i=1:1:tend
        normm = max(abs(z(i,:)));
        z(i,:)=abs(z(i,:))/normm;
     end
    mesh(t(100:tend).*1e3,x(1:ndof),z(100:tend,1:ndof).');
    set(gcf, 'Color', 'w');
    colormap jet
    h=colorbar ;
    ylabel(h,'$|p(x,t)|$','interpreter', 'latex','FontSize',15)
    hold on
    ylabel('$x [m]$','interpreter', 'latex','FontSize',15); 
    xlabel('$t [ms]$','interpreter', 'latex','FontSize',15');
    view([0,90])
    tlim=t(tend).*1e3;
    xlim([0 tlim])
    box on
end
%% Eigenmodes Animation
if eigenmodes==1
%   eigenfrequency search
    fnv=wn/2/pi/1000;
    [fi_imag,fi_real]=ginput(1)
    dfi=1/2;
    j=1;
    % searching the eigenfrequency by its real and imaginary parts
    ind_real=find(abs(real(fnv)-fi_real)<=dfi);
    ind_imag=find(abs(imag(fnv)-fi_imag)<=dfi);
    ind=intersect(ind_real,ind_imag);
    dfi=dfi/10;
    while length(ind) ~= 1 && j <=100
        ind_real=find(abs(real(fnv)-fi_real)<=dfi);
        ind_imag=find(abs(imag(fnv)-fi_imag)<=dfi);
        ind=intersect(ind_real,ind_imag);
        % updates accuracy until the only eigenfrequency find is the selected
        if isempty(ind)
            dfi=2*dfi;
        else
            dfi=dfi/10;
        end

        j=j+1;
    end
    % norm1 = max(abs(real(Vn(:,mode1))));
    disp('ind [k Hz]')
    fn=fnv(ind);
    display(fn)
    %modo=sin(x)+1j*cos(x)/2;
    mode=Vn(:,ind);
    norm = max(abs(mode));
    % numero de frames
    Nf=100;
    % numero de ciclos
    nc=1; 
    teta=0:2*nc*pi/(Nf-1):2*nc*pi; 
    figure
    % l=length(x);
    % x = [x x(l).*ones(1,l-1)+x(2:l)];
    for i=1:Nf
        y=abs(mode).*cos(angle(mode)+teta(i));
        plot(x(1:ndof),y'/norm)
    %    y=[y;y];
    %    plot(x(1:end-1),y'/norm)
        axis([min(x),max(x), -1,1]);
        xlabel('$x$ [m]', 'interpreter', 'latex', 'fontsize', 15)
        ylabel('$ \hat{P}(x,\omega_n)/ ||\hat{P}(x,\omega_n)||_{\infty}$', 'interpreter', 'latex', 'fontsize', 15)     
        set(gca,'TickLabelInterpreter','Latex','fontsize',15);
        box on
        set(gcf, 'Color', 'w');
        legend(['\omega_n = ' num2str(fn,'%0.3f') 'kHz']);
        drawnow
    end
    %%
        figure
        y=abs(mode).*sin(angle(mode));
        %y=imag(mode);
        plot(x(1:ndof),y'/norm,'LineWidth',2)
    %    y=[y;y];
    %    plot(x(1:end-1),y'/norm)
        axis([min(x),max(x), -1,1]);
        xlabel('$x$ [m]', 'interpreter', 'latex')
        ylabel('$ \hat{P}(x,\omega_n)/ ||\hat{P}(x,\omega_n)||_{\infty}$', 'interpreter', 'latex')     
        set(gca,'TickLabelInterpreter','Latex','fontsize',20);
        box on
        set(gcf, 'Color', 'w');
        legend(['\omega_n = ' num2str(fn,'%0.3f') 'kHz']);
end