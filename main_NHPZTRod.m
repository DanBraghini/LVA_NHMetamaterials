% Created on : 27/09/2021 by Danilo Braghini
%% clear 
clear all 
close all 
clc 
%% macros
% labelplot(1): wavenumber x real frequency
% labelplot(2): wavenumber x imaginary frequency
% labelplot(3): complex frequency plane
% labelplot(4): 3D plot
labelplot=[1 1 1 1];
periodic=0;
plotpassive=0;
colors=0;
%Boundary conditions of the finite structure:
%            boundary == 0 : clamped-clamped (default)/ Neumann
%            boundary == 1 : open-open/ Dirichlet
%            boundary == 2 : periodic(infinity system)
boundary=0;
% number of cells on the finite structure
ncell=40;
%FEM
ne_cell=15*2;
% excitation point for zero-state response
e='m';
%pass bands displayed
pb=4;
%% choose function
% unfolded plot of the dispersion diagrams w(k) imposing real k (traveling waves)
unfoldedDD=1;
% Bulk-Boundary correspondence or metamaterial(PBC)-metastructure(OBC) correspondence
bulkboundary=0;
% Forced frequency response for metastructure(OBC)
FRF=0;
% Transient response of the metastructure due to tone-burst kind input on
% the middle
transient=0;
% Eigenmodes simulation of the free wave modes
eigenmodes=0;
%% Crystal Set up
% Created on : 15/08/2020 by Danilo Braghini
% piezoelectric material
global  E_p rho_p L_s L_a  A_s A_a Gamma_c B1 e_p alpha_p  K_a 
% PZT-5H 
% piezoelectric material properties
E_p_ = 117e9;
% viscous damping model
%  eta =0;
eta = 0.001;
E_p = E_p_*(1+1i*eta);
e_p = 23.3;
rho_p = 7500;
alpha_p = 13.02e-9;
% Sensor and actuator geometrical properties
L_s = 5e-3;
A_s = ((5e-3)^2)*pi;
L_a =L_s;
% A_p = (0.25e-3)*pi;
A_a =A_s;
Lc = L_s+L_a;
%% Control parameters
% K_g = e_p/alpha + 10*E_r*A_r/Lc;
kappa_g = -1;
K_g = kappa_g*e_p/alpha_p;
a =0;
%% Electrical Boundary Conditions define the expressions of B
% 1) Electric-open
B1 = 0; 
% Value of K_g that opens the shunted circuit on the actuator, making the
% system Hermitian again : K_g = e_p/alpha;
% 2) Applied electric capacity
% C_p = 1.02*10^-12;
% B2 = -((e_p^2)*C_p*A_p)/((alpha^2)*A_p + alpha*C_p*L_p);
% 3) Electric-short
% B3 = -((e_p^2)*A_p)/(alpha*L_p);         
% Value of K_g that make a short circuit on the actuator: K_g = 0;
% 4) Applied feedback control
% B4 = -(e_p*A_a/L_a)*(e_p/alpha - K_g); 

Gamma_c = (e_p*A_a*K_g/L_a).*[1 -1;
                               -1 1];

K_a = -(1/L_a)*(A_a*e_p^2/alpha_p)*[1 -1
                                   -1 1 ]; 

% angular frequency vector
flim=1000e3;
wv=2*pi*(10^-4:100:flim);
fv=wv/2/pi;
%% FEM
%if bulkboundary == 1 || transient == 1 || eigenmodes == 1
    % Optimal Rayleigh viscous damping model
    [aM,aK] = function_Calibration_Rayleigh_Damping_non_local(Gamma_c,K_a,B1,L_s,L_a,...
        E_p,A_s,A_a,rho_p,e_p,alpha_p,eta,Lc,wv(10:end),ncell);
    aC=[aM aK];
    output = function_buildFEM_PZTRod(rho_p, L_s, L_a, Lc, A_s, A_a, e_p,alpha_p,E_p_,aC,Gamma_c,ne_cell,ncell,...
        'boundary',boundary,'plotpassive',plotpassive);
    M=output.M;C=output.C;K=output.K;
    Ms=output.Ms;Cs=output.Cs;Ks=output.Ks;
    if plotpassive==1
        Ks_passive=output.Ks_passive;
    end
    x=output.x;ndof=output.ndof;
    n1=output.n1;n2=output.n2;
    %% State Space of interconnection (S,H) 
    %% S
%     output=function_buildSS(Ms,Cs,Ks,ndof,ncell,n1,n2,e,a);
%     Ass=output.Ass;nx=size(Ass,1);
%     % PS: (Ms,Cs,Ks) were previously multiplied by A*rho on FEM, following 
%     % the formulation from Neurodinne&Atalla, section 4.2. The initial 
%     % dimension of external inputs to the PDE was [B/A*(volume acceleration)].
%     % So we need to multiply matrixes F (inside B1ss) and T(inside B2ss) by
%     % (B*rho) so that unitary generalized forces correspond to volume
%     % acceleration
%     B=rho*c^2;
%     B1ss=(B*rho).*output.B1ss;nw=size(B1ss,2);
%     B2ss=(B*rho).*output.B2ss;nu=size(B2ss,2);
%     C1ss=output.C1ss;nz=size(C1ss,1);
%     C2ss=output.C2ss;ny=size(C2ss,1);
%     D11ss=output.D11ss;
%     D12ss=output.D12ss;
%     D21ss=output.D21ss;
%     %% H 
%     % Feedback law in pressure by volume velocity
%     H_v=gamma_c*H_pv;
%     % Correction to volume acceleration
%     H=H_v*Nd/(1+Nd/s);
% 
%     %%  closed loop (S,H)
%     if  norm(gamma_c,1)~=0
%     %     figure
%     %     bode(H)
%     %     title('bode diagram of the feedback law')
%         [num,dem]=tfdata(H,'v');
%         if size(num,1) > size(dem,1)
%             disp('improper feedback law\n');
%         end
%         %% controlable state space realization of H
%         [Ac1,Bc1,Cc1,Dc1]=tf2ss(num,dem);
%         nxc=size(Ac1,1);
%         Acss=zeros((ncell-a)*nxc,(ncell-a)*nxc);
%         Bcss=zeros((ncell-a)*nxc,ncell-a);
%         Ccss=zeros(ncell-a,(ncell-a)*nxc);
%         %% building state space realization of H (Acss,Bcss,Ccss,Dcss)
%         % by blocks, since each cell has a decoupled feedback law
%         for k=1:ncell-a
%             lines=(k-1)*nxc+1:k*nxc;
%             Acss(lines,lines)=Ac1;
%             Bcss(lines,k)=Bc1;
%             Ccss(k,lines)=Cc1;
%         end
%         %Acss=Ac1*eye(ncell);   Bcss=Bc1*eye(ncell);
%         %Ccss=Cc1*eye(ncell);  
%         Dcss=Dc1*eye(ncell-a);
%         %% closed loop (S,H) interconnection
%         Ass_cl=[Ass+B2ss*Dcss*C2ss B2ss*Ccss
%                 Bcss*C2ss           Acss];
%         Bss_cl=[B1ss+B2ss*Dcss*D21ss
%                     Bcss*D21ss];
%         Css_cl=[C1ss+D12ss*Dcss*C2ss     D12ss*Ccss];
%         Dss_cl=D11ss+D12ss*Dcss*D21ss;
%         sys_wz = ss(Ass_cl,Bss_cl,Css_cl,Dss_cl);
%         Css=Css_cl;
%         nx_cl=size(Ass_cl,1);
%     else 
%         % only logical case other then the previous is gamma_c==0, meaning open
%         % loop system, u =0
%         sys_wz = ss(Ass,B1ss,C1ss,D11ss);
%         Css=C1ss;
%     end
    
    Ass_cl=[zeros(ndof,ndof) eye(ndof,ndof)
    -inv(Ms)*Ks -inv(Ms)*Cs];

    B=[zeros(ndof,1)
         inv(Ms)*[zeros(m,1); 1; zeros(ndof-m-1,1)]];

    C1=[eye(ndof,ndof) zeros(ndof,ndof)];
    C2=[zeros(ndof,ndof) eye(ndof,ndof)];
    D=0;
    %% verify internal stability of closed loop system
    [V, Lambda] = eig(Ass_cl);
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
%end
if unfoldedDD ==1 || bulkboundary ==1 || eigenmodes==1
%     %-------------------------------------------------------------------------%
%     %                      Direct method via SEM                              %  
%     %-------------------------------------------------------------------------%
%     % use this section to define the initial values of function fsolve
%     % angular frequency vector
%     wv=2*pi*(10^-4:100:1000e3);
%     % engeneering frequency vector
%     fv = wv./(2*pi);
% 
%     [kL_sem_PB,kL_sem_SB] = function_SEM_PZT(Gamma_c,K_a,B1,L_s,L_a,E_p,A_s,A_a,rho_p,e_p,alpha_p,wv);
% 
%     %-------------------------------------------------------------------------%
%     %                  Inverse Method via Transfer Matrix                     % 
%     %-------------------------------------------------------------------------% 
% 
%     [m1,m2,m3,m4] = function_find_mean_frequncies(kL_sem_PB);
%     m=[m1 m2 m3 m4];
%     kLv=-4*pi:8*pi/400:4*pi;
%     N=length(kLv); 
% 
%     fun= @caracteristic_equation_T_PZT_feedback;
%     options = optimoptions('fsolve','Display','off');
%     options.StepTolerance =1e-20;
%     options.OptimalityTolerance = 1e-20;
%     options.MaxIterations = 1000;
%     w_TM=zeros(4,N);
% 
%     % This loop runs through the vector kLv.
%     % PS: remember to adjust the initial conditions depending on the PC. You
%     % will need to adjust wv range also.
%     Gamma_c_aux =Gamma_c;
%     for n =1:N
%         kL=kLv(n);
%         Gamma_c = Gamma_c_aux.*exp(1i*a*kL);
%         for i=1:pb
%             w_TM(i,n) = fsolve(fun,fv(m(i))*2*pi,options);
%         end
%         w_TM(1,n) = fsolve(fun,fv(m1)*2*pi,options);
%         w_TM(2,n) = fsolve(fun,fv(m2)*2*pi,options);
%         w_TM(3,n) = fsolve(fun,fv(m3)*2*pi,options);
%         w_TM(4,n) = fsolve(fun,fv(m4)*2*pi,options);
% 
%     end    
%     Gamma_c = Gamma_c_aux;
%     % 
% % Detach frequency and damping factor for each bulk band:
%     for i=1:pb
%         w_real{i} = w_real_TM(i,:)/2/pi/1000;
%         w_imag{i} = w_imag_TM(i,:)/2/pi/1000;
%         if plotpassive==1
%             w_real_passive{i} = w_real_TM_passive(i,:)/2/pi/1000;
%             w_imag_passive{i} = w_imag_TM_passive(i,:)/2/pi/1000;
%         end
%     end
    %%
    %-------------------------------------------------------------------------%
    %                     Inverse Method via FEM                              % 
    %-------------------------------------------------------------------------%
    ni = np - 2;
    kLv=-4*pi:8*pi/400:4*pi;
    N=length(kLv); 
    K_nlocal = K;
    w_real_FEM =zeros(pb,N); w_imag_FEM = w_real_FEM;
    if plotpassive==1
        w_real_FEM_passive=w_real_FEM;w_imag_FEM_passive=w_real_FEM;
    end
    for n = 1:N
        kL = kLv(n);
        Tq = [1              zeros(1,ni) 
              zeros(ni,1)    eye(ni)
              exp(-1i*kL)    zeros(1,ni) ];

        % for non-local control  
        % Dispersion Relation uses frequency domain equations, so we can use 
        % Bloch Theorem her like we did before on SEM
        K_nlocal(n1,1) = K(n1,1) + Gamma_c(1,1)*exp(1i*a*kL);
        K_nlocal(n1,n1) = K(n1,n1) + Gamma_c(1,2)*exp(1i*a*kL);
        K_nlocal(np,1) = K(np,1) + Gamma_c(2,1)*exp(1i*a*kL);
        K_nlocal(np,n1) = K(np,n1) + Gamma_c(2,2)*exp(1i*a*kL);

        M_bar = Tq'*M*Tq;
        C_bar = Tq'*C*Tq;
        K_bar = Tq'*K_nlocal*Tq;
        w_FEM = polyeig(K_bar,1i.*C_bar,-M_bar);
        % The Quadratic Eigenvalue problem has 2*np eigenvalues, but we shall take 
        % only those with positive real(omega)
        f=find(real(w_FEM) >= 0);
        w_FEM_aux = w_FEM(f);
        w_FEM = sort(w_FEM_aux,'ComparisonMethod','real');
        w_real_FEM(:,n) = real(w_FEM(1:pb));
        w_imag_FEM(:,n) = imag(w_FEM(1:pb));

        if plotpassive==1
            K_passive_bar = Tq'*K_passive*Tq;
            w_FEM_passive = polyeig(K_passive_bar,1i.*C_bar,-M_bar);
            w_real_FEM_passive(:,n) = real(w_FEM_passive(1:pb));
            w_imag_FEM_passive(:,n) = imag(w_FEM_passive(1:pb));
            f=find(real(w_FEM_passive) >= 0);
            w_FEM_aux = w_FEM_passive(f);
            w_FEM_passive = sort(w_FEM_aux,'ComparisonMethod','real');
        end
    end
    % Detach frequency and damping factor for each bulk band and converting to kHz:
    for i=1:pb
        w_real{i} = w_real_FEM(i,:)/2/pi/1000;
        w_imag{i} = w_imag_FEM(i,:)/2/pi/1000;
        if plotpassive==1
            w_real_passive{i} = w_real_FEM_passive(i,:)/2/pi/1000;
            w_imag_passive{i} = w_imag_FEM_passive(i,:)/2/pi/1000;
        end
    end


end

%% Plots
%%
%-------------------------------------------------------------------------%
%          Unfolded Dispersion Diagram    %
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
        options{7}='w_real_passive';options{8}={w_real_passive};
        options{9}='w_imag_passive';options{10}={w_imag_passive};
    end
    function_unfolded_dispersion_diagrams(pb,kLv,w_real,w_imag,options{:});
    zlim([-1 1])
end

%%
%-------------------------------------------------------------------------%
%                   Bulk-Boundary Correspondence                          %
%-------------------------------------------------------------------------%
if bulkboundary==1 || eigenmodes==1
    clear options;
    options{1}='setplot';
    options{2}=labelplot;
    options{3}='periodic';
    options{4}=periodic;
    options{5}='plotpassive';
    options{6}=plotpassive;
    if plotpassive==1
        options{7}='w_real_passive';options{8}={w_real_passive};
        options{9}='w_imag_passive';options{10}={w_imag_passive};
    end
    function_unfolded_dispersion_diagrams(pb,kLv,w_real,w_imag,options{:});
    zlim([-1 1])
    hold on
 %    figure
    scatter(imag(wn/2/pi/1000),real(wn/2/pi/1000),'k')
 %   ylim([0 1.35])
 %    xlim([-0.3 0.3]);
 %   xlim([-2e-2 2e-2]);
%        ylabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
%      xlabel('$\Im(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
%     set(gca,'TickLabelInterpreter','Latex','fontsize',15);
%      box on
%      grid on
%      set(gcf, 'Color', 'w');
end
if FRF==1
    %% FRF on displacement by force
    F=1;
    [u_L,u_M,u_R,u_passive,x,kDg] = function_FRF_PZT_non_local(Gamma_c,K_a,B_s,L_s,L_a,E_p,A_s,A_a,rho_p,e_p,alpha_p,Lc,wv,F,ncell,a);
    u_L_SEM = u_L.';
    u_R_SEM = u_R.';
    if plotpassive==1
        u_passive_SEM = u_passive.'; 
    end
    figure
    plot(fv/1000,20*log10(abs(u_L_SEM(:,1))),'r','MarkerSize',8,'LineWidth',1)
    hold on
    plot(fv/1000,20*log10(abs(u_R_SEM(:,end))),'k','MarkerSize',8,'LineWidth',1)
    if plotpassive==1
        P_passive = output.ps_passive.';
        plot(fv/1000,20*log10(abs(P_passive(:,1))),'r--','MarkerSize',8,'LineWidth',1)
        plot(fv/1000,20*log10(abs(P_passive(:,end))),'k--','MarkerSize',8,'LineWidth',1)
    end
    ylabel('$ u / F $ [dB]', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    %legend('L','R','L filtered','R Filtered','passive')
    %xlim([0 flim/1000])
    box on
    set(gcf, 'Color', 'w');
    grid on
end
if transient==1
    %% Transient Response: applied force
    % Nc is the number of circles from central frequency fc.
    % T2 define the envelope frequency f2, which has Nc circles within.
    % Nt is the number of entries on time vector t. The heavside function
    % defines a window on half a period of the envelope.
    % time discretization period dt depends on N, in such manner that
    % the total time of the simulation T can be modified if nedeed.
    Nt=2*4096;
    % T=0.01; % total analysis time
    Ap=1;
    fc=300000;
    wc = 2*pi*fc;
    Nc=15; 
    Tc=1/fc;
    T2=Nc*Tc*2; 
    f2=1/T2; w2=2*pi*f2;
    dt=T2/Nt;
    % dt = T/N;
    t=0:dt:(Nt-1)*dt;
    % t=0:dt:1000e-6;
    F=Ap*sin(wc*t).*sin(w2*t).^2.*heaviside(T2/2-t); 
     df=1/t(end); 
    fm=1/dt;
    fv=0:df:fm/2;
    % fN=N*df/2;
    % f=0:df:fN-1;
    Fw=1/Nt*fft(F);
    % U = zeros(ndof,N);
   if e == 'm'
        m = floor(ndof/2);
        if mod(ndof,2) == 0
             F = [zeros(ndof-m-1,1);1;zeros(ndof-m,1)];
        else
            F = [zeros(m,1);1;zeros(m,1)];
        end
    % excitation on the left end
    elseif e == 'l'
        F = [1; zeros(ndof-1,1)];
    % excitation on the right end    
    elseif e=='r'
        F = [zeros(ndof-1,1);1];
    else 
        disp('invalid force position')
        return
   end

    A=[zeros(ndof,ndof) eye(ndof,ndof)
    -inv(Ms)*Ks -inv(Ms)*Cs];
    B=[zeros(ndof,1)
         inv(Ms)*F];
    C1=[eye(ndof,ndof) zeros(ndof,ndof)];
    C2=[zeros(ndof,ndof) eye(ndof,ndof)];
    D=0;
    sys_x = ss(A,B,C1,D);
    t=0:dt:(Nt-1)*dt;
    w=Ap*sin(wc*t).*sin(w2*t).^2.*heaviside(T2/2-t); 
    %% zero state response of the closed loop for input w
    Y=lsim(sys_x,w',t,zeros(2*ndof,1));
    %% plots
    Nplots=100;
    tend=length(t)/5;
    figure
    plot(t*10^6,F)
    %title('external force dynamics')
    xlabel('t [\mu s]')
    ylabel('F (t)')
    figure
    plot(fv/1000,abs(Fw(1:Nt/2)))
    %title('external force spectrum')
    xlabel('f(kHz)')
    ylabel('|F(\omega)|')
    xlim([0 fc/1000 + 200])
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
end