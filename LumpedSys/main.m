%% clear 
clear all 
close all 
clc 
%echo on
%%
% boundary == 0 : free-free (default)
% bundary == 1: fixed-fixed
% boundary == 2 : WABC
symmetricdimmer=0;%requires ndof=3*ncell-ncell+1 and boundary=0
boundary=1;
eigenmode=1;
bulkboundary=1;
Bode=0;
dispersion=1;
FRF3D=0;
s=tf('s');
%% lumped-parameter system set up
% m = mass 
% k = string elastic coefficient
% r=range of non-reciprocal coupling, r=0 means next-neighborhoud couple
% e = excitation point for the external 
k=10; m1=1; m2=m1;gp=-10;
gd=0;gi=0;gdd=0;%-11;
r=0;
b=0;%0.01;
e='m';
ndof=9;
% angular frequency vector limit
flim=3.5e0;
wv=2*pi*(0.01:0.002:flim);
fv=wv/2/pi;
%%
% D=[s^2*m+2*c*s+2*k  -s*c-k  0  0
%            -s*c-k+g    s^2*m+2*c*s+2*k-g -s*c-k 0 
%              0               -s*c-k+g        s^2*m+2*c*s+2*k-g -s*c-k
%             0                    0                      -s*c-k-g             s^2*m+2*c*s+2*k-g];
% P=simplify(inv(D));

%% periodic (infinity) system with one mass on the unit cell
if dispersion==1 || bulkboundary==1 || eigenmode==1
    muvec=pi*linspace(-1,1,701);
    omega_a1=zeros(size(muvec));omega_a2=omega_a1;
    for i=1:701
        cmu=1-cos(muvec(i));
        emu=1-exp(1i*muvec(i));
 %       a1=(-1i/m1)*(2*b*cmu-gd*emu);
 %      a2=(-1/m1)*(2*k*cmu-gp*emu);
  %     omega_a1(i)=-a1/2+sqrt(a1^2-4*a2)/2;
  %      omega_a2(i)=-a1/2-sqrt(a1^2-4*a2)/2;
        omega_a1(i)=sqrt(2*k*(1-cos(muvec(i))/m1 ) - (gp/m1)*(1-exp(1i*muvec(i)) )*exp(1i*muvec(i)*r ) );
       omega_a2(i)=-sqrt(2*k*(1-cos(muvec(i))/m1 ) - (gp/m1)*(1-exp(1i*muvec(i)) )*exp(1i*muvec(i)*r ) );
    a1=(1-gdd/m1*exp(1i*muvec(i)*r)*(1+exp(1i*muvec(i))))*1i;
    a2=2*b/m1*cmu-gd/m1*exp(1i*muvec(i)*r)*emu;
    a3=(gp/m1*exp(1i*muvec(i)*r)*emu-2*k/m1*cmu)*1i;
    a4=gi/m1*exp(1i*muvec(i)*r)*emu;
    omega_r(:,i)=roots([a1 a2 a3 a4]);
    end
     wn=sqrt(k/m1);
     Omega_r=omega_r/wn;
     omega_a1=omega_a1/wn;
     omega_a2=omega_a2/wn;
%%
    %Omega=omega_a1./wn;
    
    figure
   % plot3(real(omega_r(1,:)),imag(omega_r(1,:)),muvec,'r','LineWidth',1.5)
    plot3(real(Omega_r(2,:)),imag(Omega_r(2,:)),muvec/pi,'r','LineWidth',1.5)
       hold on
   plot3(real(Omega_r(3,:)),imag(Omega_r(3,:)),muvec/pi,'r','LineWidth',1.5)
    ylabel('$\Omega_I$','interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Omega_R$','interpreter', 'latex', 'fontsize', 15)
    zlabel('$\mu/\pi$','interpreter', 'latex', 'fontsize', 15)
    box on
    grid on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',20);
%%     
    figure
    plot(muvec/pi,real(Omega_r(2,:)),'r','LineWidth',1.5)
    %hold on
   %plot(real(omega_a2),muvec,'r','LineWidth',1.5)
    ylabel('$\Omega_R$','interpreter', 'latex', 'fontsize', 20)
    xlabel('$\mu/\pi$','interpreter', 'latex', 'fontsize', 20)
    box on
    grid on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',20;

    figure
    plot(muvec/pi,imag(Omega_r(2,:)),'r','LineWidth',1.5)
    %hold on
    %plot(imag(omega_a2),muvec,'r','LineWidth',1.5)
    ylabel('$\Omega_I$','interpreter', 'latex', 'fontsize', 15)
    xlabel('$\mu/\pi$','interpreter', 'latex', 'fontsize', 15)
    box on
    grid on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',20);
    %%
     figure
    plot(real(Omega_r(1,:)),imag(Omega_r(1,:)),'r','LineWidth',1.5)
    hold on
    plot(real(Omega_r(2,:)),imag(Omega_r(2,:)),'r','LineWidth',1.5)
    plot(real(Omega_r(3,:)),imag(Omega_r(3,:)),'r','LineWidth',1.5)
    xlabel('$\omega_R$','interpreter', 'latex', 'fontsize', 15)
    ylabel('$\omega_I$','interpreter', 'latex', 'fontsize', 15)
    box on
    grid on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',20);

end
if bulkboundary==1 || eigenmode==1 || Bode==1
    %% gains
    gamma_c=[gp gi gd];
    % quasiperiodicity parameters
    theta=0;phi=0;
    %[x,~,so_,Vn,sys_wz,nx_cl] = function_LumpedEigenfreq(ndof,m1,m2,k,eta,flim,boundary,e,gamma_c,theta,phi);
    %% Solving the eigenproblem under PBC
    output = function_LumpedEigenfreq(100,m1,m2,k,b,flim,2,e,gamma_c,theta,phi);
    omegan_wabc=output.wn;Vn_wabc=output.Vn;sys_wx_wabc=output.sys_wz;nx_cl_wabc=output.nx_cl;ndof_wabc=output.ndof;x_wabc=output.x;
    [Omegan_wabc,ind] = sort(omegan_wabc/wn,'ComparisonMethod','real');
    Vn_wabc = Vn_wabc(:,ind); 
       % solving the eigenproblem for the passive structure
%     gamma_c=[0*gp 0*gi 0*gd];
%     output= function_LumpedEigenfreq(ndof,m1,m2,k,b,flim,boundary,e,gamma_c,theta,phi);
%     wn_passive=output.wn;Vp=output.Vn;sys_passive=output.sys_wz;
%     [wn_passive,ind] = sort(wn_passive,'ComparisonMethod','real');
%     Vn = Vp(:,ind); 
    %% Computing dynamic stiffness matrix
    % number of elements
    if ndof>1
        mvec=zeros(ndof,1);
        mvec(1:2:ndof,1)=m1; mvec(2:2:ndof,1)=m2;
        Ms=diag(mvec);
         if symmetricdimmer==1
        % ndof=3*ncell-ncell+1 and boundary==1
              Ms(1,1)=m1/2;Ms(end,end)=Ms(1,1);
         end
        D2_=toeplitz([2 -1 zeros(1,ndof-2)]);
        D2_(1,1)=1; D2_(ndof,ndof)=1;
        Ks=k.*D2_;
       %viscous damping model
        R=b.*D2_;
        % Feedback Matrix
        Gp=gp.*(diag(-1*ones(1,ndof))+diag(ones(1,ndof-1),-1));
        Gp(1,:)=zeros(1,ndof);
        Gd=gd.*(diag(-1*ones(1,ndof))+diag(ones(1,ndof-1),-1));
        Gd(1,:)=zeros(1,ndof);
    else
    disp('Warning: only one mass isnt enough to have feedback' )
    end
    % structural/histerisis damping model
    %fm=flim/2;
    %Cs=eta.*Ks/fm/2/pi;
    %free-free
    ndof_n=ndof;
    %fixed-fixed
    if boundary==1
        Ms(1,:)=[];Ms(:,1)=[];Ms(end,:)=[];Ms(:,end)=[];
        Ks(1,:)=[];Ks(:,1)=[];Ks(end,:)=[];Ks(:,end)=[];
        R(1,:)=[];R(:,1)=[];R(end,:)=[];R(:,end)=[];
        Gp(1,:)=[];Gp(:,1)=[];Gp(end,:)=[];Gp(:,end)=[]; 
        Gd(1,:)=[];Gd(:,1)=[];Gd(end,:)=[];Gd(:,end)=[]; 
        ndof_n=ndof-2;
        Gp(1,:)=zeros(1,ndof_n);
        Gd(1,:)=zeros(1,ndof_n);
    % WABC
    elseif boundary==2
        Ms(:,1)=Ms(:,1)+Ms(:,end);
        Ms(1,:)=Ms(1,:)+Ms(end,:);
        Ms(:,end)=[];Ms(end,:)=[]; 
        Ks(:,1)=Ks(:,1)+Ks(:,end);  
        Ks(1,:)=Ks(1,:)+Ks(end,:);
        Ks(:,end)=[];Ks(end,:)=[];
        R(:,1)=R(:,1)+R(:,end);  
        R(1,:)=R(1,:)+R(end,:);
        R(:,end)=[];R(end,:)=[];
        ndof_n=ndof-1;
    end
    %% Assemblyng dynamic stiffness matrix and solving u(\omega)  
    W=length(wv);
    Dg=zeros(ndof_n,ndof_n,W);
     %% define where the external load is applied
     % F=1 for FRF
     F=1;
    % force in the middle
    if e == 'm'
        m = floor(ndof_n/2);
        if mod(ndof_n,2) == 0
             F = [zeros(ndof_n-m-1,1);F;zeros(ndof_n-m,1)];
        else
            F = [zeros(m,1);F;zeros(m,1)];
        end
    % force on the left end
    elseif e == 'l'
        F = [F; zeros(ndof_n-1,1)];
    % force on the right end    
    elseif e=='r'
        F = [zeros(ndof_n-1,1);F];
    else 
        disp('invalid force position')
        return
    end
    q=zeros(ndof_n,W);
    % This loop runs through the frequency vector
    % Computing the FRF through the dynamic stiffness matrix
    %     for i=1:W
    %        w=wv(i);
    %        Dm=-w^2*Ms-1i*w*Cs+Ks+Gs;
    %       q(:,i)=Dm\F;
    %     end
    % q1m=abs(q(1,:));
    % qNm=abs(q(end,:));
    % q1p=atan2d(imag(q(1,:)),real(q(1,:)));
    %qNp=atan2d(imag(q(end,:)),real(q(end,:)));
    %end
     % or from the state space 4 matrix representation
     nw=1;nz=ndof_n;
    A=[zeros(ndof_n) eye(ndof_n)
          -Ms\(Ks+Gp)    -Ms\(R+Gd)];
    B=[zeros(ndof_n,1)
            Ms\F];
    C=[eye(ndof_n)  zeros(ndof_n)];
    D=zeros(nz,nw);
    sys_simo=ss(A,B,C,D);
    P=tf(sys_simo);
    P1=P(1,1);Pn=P(ndof_n,1);
    % BODE diagram of tf for output u(1)
    [mag,phase]=bode(P1,wv);
    q1m=vec(mag);
    q1p=vec(phase);
    % BODE diagram of tf for output u(ndof)
    [mag,phase]=bode(Pn,wv);
    qNm=vec(mag);
    qNp=vec(phase);

     nz2=1;
     C=[ 1 zeros(1,ndof_n-2) -1 zeros(1,ndof_n)];
     D=zeros(nz2,nw);
     sys_siso=ss(A,B,C,D);
     %output1= function_stability_c(A,eye(size(A,1)))
    % output2=function_H2Norm(A,B,C,D)
    %% s=-i \omega
%((-i\omega)^2 M -i\omega C +K+G)q(\omega) = 0
if b ~=0
    [V,so]=polyeig(Ks+Gp,R,Ms);
    [V_passive,so_passive]=polyeig(Ks,R,Ms);
else
    [V,L]=eig(Ks+Gp,-Ms);
    so=[sqrt(diag(L));-sqrt(diag(L))];
    Vn=[V V];
    [V_passive,L]=eig(Ks,-Ms);
    so_passive=[sqrt(diag(L));-sqrt(diag(L))];
    Vn_passive=[V_passive V_passive];
end
omega_n=so.*1i;
[omega_n,ind]=sort(omega_n/wn,'ComparisonMethod','real');
Vn = Vn(:,ind); 
omega_npassive=so_passive.*1i;
[omega_npassive,ind]=sort(omega_npassive/wn,'ComparisonMethod','real');
Vn_passive=Vn_passive(:,ind);
%[V,L]=eig(Ks,-Ms);
%so_h=[sqrt(diag(L));-sqrt(diag(L))];
%%
%     figure(6)
%     scatter(real(so),imag(so),'kx','LineWidth',1.5)
%     hold on
%     scatter(real(so_h),imag(so_h),'bx','LineWidth',1.5)
%     ylim([0 8])
%     grid on
%     title('Eigenvalues', 'interpreter', 'latex', 'fontsize', 15)
%     ylabel('$\Im ( s )$', 'interpreter', 'latex', 'fontsize', 15)
%     xlabel('$\Re ( s )$', 'interpreter', 'latex', 'fontsize', 15')
%     set(gcf, 'Color', 'w');
%     set(gca,'TickLabelInterpreter','Latex','fontsize',20);
%     box on
    %%
    if Bode==1
        figure
        plot(fv,20*log10(q1m),'r','MarkerSize',8,'LineWidth',1)
        hold on
        plot(fv,20*log10(qNm),'k','MarkerSize',8,'LineWidth',1)
        hold off
        ylabel('$ |\hat{q} / \hat{F}| (i \omega)$ [dB]', 'interpreter', 'latex', 'fontsize', 15)
        xlabel('$f$ [Hz]', 'interpreter', 'latex', 'fontsize', 15)
        set(gca,'TickLabelInterpreter','Latex','fontsize',15);
        legend('L','R')
        xlim([0 flim])
        box on
        set(gcf, 'Color', 'w');
        set(gca,'TickLabelInterpreter','Latex','fontsize',15);
        grid on
        % figure(4)
        % plot(fv,q1p,'r','MarkerSize',8,'LineWidth',1)
        % hold on
        % plot(fv,qNp,'k','MarkerSize',8,'LineWidth',1)
        % hold off
        % ylabel('$ \hat{q} / \hat{F} (i \omega) $', 'interpreter', 'latex', 'fontsize', 15)
        % xlabel('$f$ [Hz]', 'interpreter', 'latex', 'fontsize', 15)
        % set(gca,'TickLabelInterpreter','Latex','fontsize',15);
        % legend('L','R')
        % xlim([0 flim])
        % box on
        % set(gcf, 'Color', 'w');
        % set(gca,'TickLabelInterpreter','Latex','fontsize',15);
        % grid on
    end
end
%%
if FRF3D==1
    % Computing the FRF through the dynamic stiffness matrix
    for i=1:W
       w=wv(i);
       Dm=-w^2*Ms-1i*w*R+Ks+Gp;
      q(:,i)=Dm\F;
    end
    figure(5)
    xs=1:1:ndof;
    if boundary==1
        q=[zeros(1,W); q; zeros(1,W)];
    elseif boundary==2
        xs(:,end)=[];
    end
    mesh(fv,xs,20*log10(abs(q)))
    %contour(fv(10:end)/1000,xs,20*log10(abs(u_M_SEM)))
    %view([0 -90])
    shading INTERP
    colormap jet
    h=colorbar;
    ylabel(h,'$| \hat{u}(x) |$ [dB]','interpreter', 'latex','FontSize', 15)
    %w_real_TM_4PB(end)/2/pi/1000 =  919.9540
    ylabel('$x$ [m]', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
    %zlabel('$| \hat{u}(x) | [dB]$', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    box on
    set(gcf, 'Color', 'w');
end
if dispersion==1 && eigenmode==1
    figure(7)
    scatter(imag(omega_n),real(omega_n),'ko','LineWidth',1.5)
    hold on
   plot(imag(omega_a1),real(omega_a1),'r','LineWidth',1.5)
    plot(imag(omega_a2),real(omega_a2),'r','LineWidth',1.5)

  %  scatter(real(so_h),imag(so_h),'bx','LineWidth',1.5)
   %scatter(real(so_wabc),imag(so_wabc),'ro')
    scatter(imag(omega_npassive),real(omega_npassive),'bo','LineWidth',1.5)
    ylim([0 3])
    grid on
    %title('Eigenvalues', 'interpreter', 'latex', 'fontsize', 15)
    ylabel('$\Omega_R$', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Omega_I$', 'interpreter', 'latex', 'fontsize', 15')
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',20);
    box on
    %% Eigenmodes animation
x=1:1:ndof_n;
if eigenmode==1
   % ylim([-20 0])
        [fi_real,fi_imag]=ginput(1);
       % fi_imag=0;
       % fi_real=0;
        dfi=1/2;
        j=1;
        % searching the eigenfrequency by its real and imaginary parts
        wnv=omega_n;
        ind_real=find(abs(real(wnv)-fi_real)<=dfi);
        ind_imag=find(abs(imag(wnv)-fi_imag)<=dfi);
        ind=intersect(ind_real,ind_imag);
        dfi=dfi/10;
        while length(ind) ~= 1 && j <=100
            ind_real=find(abs(real(wnv)-fi_real)<=dfi);
            ind_imag=find(abs(imag(wnv)-fi_imag)<=dfi);
            ind=intersect(ind_real,ind_imag);
            % updates accuracy until the only eigenfrequency find is the selected
            if isempty(ind)
                dfi=2*dfi;
            else
                dfi=dfi/10;
            end
    
            j=j+1;
        end
        if isempty(ind)
            j=1;
            % searching the eigenfrequency by its real and imaginary parts
            wnv=omega_npassive;
            ind_real=find(abs(real(wnv)-fi_real)<=dfi);
            ind_imag=find(abs(imag(wnv)-fi_imag)<=dfi);
            ind=intersect(ind_real,ind_imag);
            dfi=dfi/10;
            while length(ind) ~= 1 && j <=100
                ind_real=find(abs(real(wnv)-fi_real)<=dfi);
                ind_imag=find(abs(imag(wnv)-fi_imag)<=dfi);
                ind=intersect(ind_real,ind_imag);
                % updates accuracy until the only eigenfrequency find is the selected
                if isempty(ind)
                    dfi=2*dfi;
                else
                    dfi=dfi/10;
                end

                j=j+1;
            end
        end
        
        % norm1 = max(abs(real(Vn(:,mode1))));
        disp('ind [rad/s]')
        wo=wnv(ind);
        display(wo)
        %modo=sin(x)+1j*cos(x)/2;
        mode=Vn(:,ind);
        %wn=-1i*so(end);
        normm = max(abs(mode));
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
            scatter(x(1:ndof_n),y'/normm)
        %    y=[y;y];
        %    plot(x(1:end-1),y'/norm)
            axis([min(x),max(x), -1,1]);
            xlabel('$n$', 'interpreter', 'latex', 'fontsize', 15)
            ylabel('$ \hat{g}(\omega_n) $', 'interpreter', 'latex', 'fontsize', 15)     
            set(gca,'TickLabelInterpreter','Latex','fontsize',15);
            box on
            set(gcf, 'Color', 'w');
            legend(['\omega_n = ' num2str(wo,'%0.3f') 'rad/s']);
            drawnow
           
        end
end
end
%% zero state response of the closed loop for input w
% wc=0.5;
% fc=wc/2/pi;
% Nc=30; 
% Tc=1/fc;
% T2=Nc*Tc*2; 
% f2=1/T2; w2=2*pi*f2;
% dt=Tc/100;
% Nt=8*9046;
% t=0:dt:(Nt-1)*dt;
% tend=round(length(t)/10);
% 
% w=double(sin(wc*t).*sin(w2*t).^2.*heaviside(T2/2-t)); 
% df=1/t(end); fm=1/dt;
% fvec=0:df:fm/2;
% 
% Fw=1/Nt*fft(w);
% % 
% figure
% subplot(1,2,1)
% plot(t(1:tend),w(1:tend))
% xlabel('$t$ [s]', 'interpreter', 'latex','fontsize', 15)
% %xaxis([0 100])
% box on
% grid on
% title('External Force Dynamics', 'interpreter', 'latex','fontsize', 15)
% subplot(1,2,2)
% plot(2*pi*fvec,abs(Fw(1:Nt/2)))
% xlabel('$\omega $ [rad/s]', 'interpreter', 'latex','fontsize', 15)
% xlim([0 2])
% title('External Force Spectrum', 'interpreter', 'latex','fontsize', 15)
% box on
% grid on
% set(gcf, 'Color', 'w');
% %% zero state response of the closed loop for input w
% 
% ni=size(A,1);
% z=lsim(sys_wz,w,t,zeros(ni,1));
% %%
% figure
% Nplots=100;
% DI=round(tend/Nplots);
% %i=tend;
% for i=1:10:tend
%     scatter3(x(1:size(z,2)),(t(i)*ones(1,size(z,2))),z(i,:),'b');
%      hold on
% end
% xlabel('$n$', 'interpreter', 'latex', 'fontsize', 15)
% ylabel('$t$ [s]', 'interpreter', 'latex', 'fontsize', 15)
% zlabel('$q_n(t)$ [m]', 'interpreter', 'latex', 'fontsize', 15)
% hold off
% box on
% grid on
% set(gcf, 'Color', 'w');
% %%
% figure
% %tend=length(t);
% plot(t(1:tend),z(1:tend,1),'b')
% hold on
% plot(t(1:tend),z(1:tend,end),'k')
% set(gcf, 'Color', 'w');
% xlim([0,t(tend)])
% xlabel('$t [s]$','interpreter', 'latex','FontSize',15');
% legend('$q_1(t) [m]$','$q_N(t) [m]$','interpreter', 'latex','FontSize',15');
% grid on

%% Stabilizable dynamic output controller synthetis
%options.hinf=1e+20;
%output = function_Hinf_Dynamic_OutputFeedback(Ass,B1ss,B2ss,C1ss,C2ss,D11ss,D12ss,D21ss)
%% Internal stalibility test
%  output.feas=0 = > internally unstable
% output.feas=1 = > internally stable
% output = function_stability_c(Ass)
% output = function_Hinf_Optimal_Controller_c(Ass,B1ss,B1ss,C1ss,D11ss,D12ss)
% %%
% function_stabilizability_c(Ass,B2ss)
%  %%
%  % sample frequency for x domain
% fxsamp=1;% The distance between each mass is 0.5 for the dimer and 1 for monoatomic model 
% %(unitary distance is assumed as the length of the cell)
% % total length of x domain
% Nx = length(x);
% % resolution of spacial frequency (kx) mapped from x
% dkx = fxsamp/Nx;
% % Nysquist wavenumber
% kxN = fxsamp/2;
% % if Nx is odd
% if mod(Nx,2) ~= 0
% %    kx=-kxN+dkx/2:dkx:kxN-dkx/2;
%     kx = (1/dkx)*(-(Nx-1)/2:(Nx-1)/2)/Nx;
% else
%     kx = -kxN:dkx:kxN-dkx;
% end
% kx=2*pi*kx;
% p=0:.1:1;
% for j=1:size(Vn,2)
%     Vn_til(:,j)=1/Nx*fft(Vn(:,j));
% end
% 
% scatter(imag(wn),real(wn),'k')
% hold on
% normm = max(max(abs(Vn_til)));
% contour(wn,kx/pi,abs(Vn_til)./normm,p)
% set(gca,'TickLabelInterpreter','Latex','fontsize',15);
% box on
% grid on
% set(gcf, 'Color', 'w');
% colormap('jet')
% shading interp
% %ylim([0 2])
%echo off