%% clear 
clear all 
close all 
clc 
%echo on
%%
% boundary == 0 : free-free (default)
% bundary == 1: fixed-fixed
% boundary == 2 : WABC
boundary=1;
s=tf('s');
%% lumped-parameter system set up
% m = mass 
% k = string elastic coefficient
% range of non-reciprocal coupling, r=0 means next-neighborhoud couple
% e = excitation point for the external 
k=10; m1=1; m2=1;
%gp=-35:1:8;
gp=-35:1:35;
gd=gp;
r=0;
b=.01;
e='m';
ndof=9;
% angular frequency vector limit
flim=3.5e0;
wv=2*pi*(1e-2:2e-4:flim);
fv=wv/2/pi;
% number of elements
%     if ndof>1
%         mvec=zeros(ndof,1);
%         mvec(1:2:ndof,1)=m1; mvec(2:2:ndof,1)=m2;
%         Ms=diag(mvec);
%         D2_=toeplitz([2 -1 zeros(1,ndof-2)]);
%         D2_(1,1)=1; D2_(ndof,ndof)=1;
%         Ks=k.*D2_;
%        %viscous damping model
%         Cs=b.*D2_;
%     else
%     disp('Warning: only one mass isnt enough to have feedback' )
%     end
%     %fixed-fixed
    if boundary==1
        ndof_n=ndof-2;
%         Ms(1,:)=[];Ms(:,1)=[];Ms(end,:)=[];Ms(:,end)=[];
%         Ks(1,:)=[];Ks(:,1)=[];Ks(end,:)=[];Ks(:,end)=[];
%         Cs(1,:)=[];Cs(:,1)=[];Cs(end,:)=[];Cs(:,end)=[];
     % WABC
   elseif boundary==2
       ndof_n=ndof-1;
%         Ms(:,1)=Ms(:,1)+Ms(:,end);
%         Ms(1,:)=Ms(1,:)+Ms(end,:);
%         Ms(:,end)=[];Ms(end,:)=[]; 
%         Ks(:,1)=Ks(:,1)+Ks(:,end);  
%         Ks(1,:)=Ks(1,:)+Ks(end,:);
%         Ks(:,end)=[];Ks(end,:)=[];
%         Cs(:,1)=Cs(:,1)+Cs(:,end);  
%         Cs(1,:)=Cs(1,:)+Cs(end,:);
%         Cs(:,end)=[];Cs(end,:)=[];
    elseif boundary==0
        ndof_n=ndof;
    end
% Dynamic controllers add aditional poles
%so=zeros(length(gp),2*ndof_n);
area=zeros(size(gp,2),size(gd,2));lmi_feas=area;

q1m=zeros(size(gp,2),size(gd,2),size(wv,2));q1p=q1m;qnm=q1m;qnp=q1m;
for j=1:length(gd)
    for i=1:length(gp)
        gamma_c=[gp(i) 0  gd(j)];
        output=function_CL4MSS(ndof,m1,m2,k,b,flim,boundary,e,gamma_c);
        %% Computing FRF for each gain
         % Feedback Matrix
    %      Gs=gp(i).*(diag(-1*ones(1,ndof))+diag(ones(1,ndof-1),-1));
    %      Gs(1,:)=zeros(1,ndof);
    %          %fixed-fixed
    %     if boundary==1
    %         Gs(1,:)=[];Gs(:,1)=[];Gs(end,:)=[];Gs(:,end)=[]; 
    %         ndof_n=ndof-2;
    %         Gs(1,:)=zeros(1,ndof_n);
    %     % WABC
    %     elseif boundary==2
    %         ndof_n=ndof-1;
    %     end
             %% define where the external load is applied
    %  % F=1 for FRF
    %  F=1;
    % % force in the middle
    % if e == 'm'
    %     m = floor(ndof_n/2);
    %     if mod(ndof_n,2) == 0
    %          F = [zeros(ndof_n-m-1,1);F;zeros(ndof_n-m,1)];
    %     else
    %         F = [zeros(m,1);F;zeros(m,1)];
    %     end
    % % force on the left end
    % elseif e == 'l'
    %     F = [F; zeros(ndof_n-1,1)];
    % % force on the right end    
    % elseif e=='r'
    %     F = [zeros(ndof_n-1,1);F];
    % else 
    %     disp('invalid force position')
    %     return
    % end
        %% closed-loop state-space 4-matrix representation SIMO for FRFs
    %    nw=1;nz=ndof_n;
    %      A=[zeros(ndof_n) eye(ndof_n)
    %           -Ms\(Ks+Gs)    -Ms\Cs];
    %     B=[zeros(ndof_n,1)
    %             Ms\F];
    %     C=[eye(ndof_n)  zeros(ndof_n)];
    %     D=zeros(nz,nw);
    %     sys_simo=ss(A,B,C,D);
    %     
        sys_simo=output.CLsys;
        [A,B,C,D]=ssdata(sys_simo);
        P=tf(sys_simo);
        % Taking BODE diagrams of first and last DOFs
        P1=P(1,1);Pn=P(ndof_n,1);
        % BODE diagram of tf for output q(1)
        [mag,phase]=bode(P1,wv);
        q1m(i,j,:)=vec(mag);
        q1p(i,j,:)=vec(phase);
        % BODE diagram of tf for output q(n)
        [mag,phase]=bode(Pn,wv);
        qnm(i,j,:)=vec(mag);
        qnp(i,j,:)=vec(phase);
        dw=wv(2)-wv(1);
        area(i,j)=(trapz(qnm(i,j,:))-trapz(q1m(i,j,:)))*dw;
        area_rms(i,j)=(trapz(qnm(i,j,:).^2)-trapz(q1m(i,j,:).^2))*dw;
        sys_siso{i,j}.arearms=area_rms(i,j);
        sys_siso{i,j}.area=area(i,j);
        %% state-space 4-matrix representation SISO for H2 norm
        nw=1;nz=1;
       C=[1 zeros(1,ndof_n-2) -1 zeros(1,size(A,2)-ndof_n)];
        D=zeros(nz,nw);
        sys_siso{i,j}.statespace=ss(A,B,C,D);
        sys_siso{i,j}.norm2_matlab=norm(sys_siso{i,j}.statespace,2);
       % output=stability_c(A,eye(size(A)),'solver','mosek');
       % allowing marginal stability
       %     output=stability_c(A,zeros(size(A)),'solver','mosek');
    
      %  sys_siso{i}.stable=output.feas;
      %  sys_siso{i}.lyap=[];
       % if output.feas==1
      %      sys_siso{i}.lyap=output.P;
    %    end
    %output=h2_norm_c(A,B,C,'solver','mosek');
     %f output.feas==1
      %      lmi_feas(i)=1;
       %     sys_siso{i}.h2bound=output.h2;
    %end
        %% Computing Eigenspectum for each gain (free-vibration problem)
        % s=-i \omega
        %((-i\omega)^2 M -i\omega C +K+G)q(\omega) = 0
    %     if b ~=0
    %         [V,so(i,:)]=polyeig(Ks+Gs,Cs,Ms);
    %     else
    %         [V,L]=eig(Ks+Gs,-Ms);
    %         so(i,:)=[sqrt(diag(L));-sqrt(diag(L))];
    %     end
    so{i,j}=output.eig;
    
    end
end
%% Bode plot
% f0= 20;
% qN0m=qnm(f0,f0,:);
% q10m=q1m(f0,f0,:);
% q10p=q1p(f0,f0,:);
% qN0p=qnp(f0,f0,:);
% figure
% plot(fv,20*log10(vec(q10m).'),'r','MarkerSize',8,'LineWidth',1)
% hold on
% plot(fv,20*log10(vec(qN0m).'),'k','MarkerSize',8,'LineWidth',1)
% hold off
% ylabel('$ |\hat{q} / \hat{F}| (i \omega)$ [dB]', 'interpreter', 'latex', 'fontsize', 15)
% xlabel('$f$ [Hz]', 'interpreter', 'latex', 'fontsize', 15)
% set(gca,'TickLabelInterpreter','Latex','fontsize',15);
% legend('L','R')
% xlim([0 flim])
% box on
% set(gcf, 'Color', 'w');
% set(gca,'TickLabelInterpreter','Latex','fontsize',15);
% grid on
% figure
% plot(fv,vec(q10p)','r','MarkerSize',8,'LineWidth',1)
% hold on
% plot(fv,vec(qN0p)','k','MarkerSize',8,'LineWidth',1)
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
%% Spectrum
%     figure
%     scatter(real(so),imag(so),'kx','LineWidth',1.5)
%     hold on
%     plot(real(-1i.*omega_a),imag(-1i.*omega_a),'r','LineWidth',1.5)
%     scatter(real(so_h),imag(so_h),'bx','LineWidth',1.5)
%     scatter(real(so_wabc),imag(so_wabc),'r')
%     %ylim([0 8])
%     grid on
%     title('Eigenvalues', 'interpreter', 'latex', 'fontsize', 15)
%     ylabel('$\Im ( s )$', 'interpreter', 'latex', 'fontsize', 15)
%     xlabel('$\Re ( s )$', 'interpreter', 'latex', 'fontsize', 15')
%     set(gcf, 'Color', 'w');
%     set(gca,'TickLabelInterpreter','Latex','fontsize',20);
%     box on
figure
%wo=so.*i;
for i=1:length(gp)
 %   for j=1:length(gd)
    %  scatter3(gp(i).*ones(1,size(so{i,j},1)),gd(j).*ones(1,size(so{i,j},1)),real(so{i,j}),2,'k')
      scatter(gp(i).*ones(1,size(so{i,i},1)),real(so{i,i}),'.k')
       hold on
    %end
end
%plot(real(wo(1:nx,:)),phi/pi)
grid on
title('Modulated spectrum', 'interpreter', 'latex', 'fontsize', 15)
zlabel('$\Re(s)$', 'interpreter', 'latex', 'fontsize', 15)
xlabel('$g_p$', 'interpreter', 'latex', 'fontsize', 15')
ylabel('$g_i$', 'interpreter', 'latex', 'fontsize', 15')
set(gcf, 'Color', 'w');
box on
set(gca,'TickLabelInterpreter','Latex','fontsize',15);

%  figure
% for i=1:length(gp)
%     plot3(gp(i).*ones(1,size(A,1)),imag(so(i,:)),real(so(i,:)),'.k')
%     hold on
% end
% grid on
% title('Modulated spectrum', 'interpreter', 'latex', 'fontsize', 15)
% xlabel('$g_p$', 'interpreter', 'latex', 'fontsize', 15')
% ylabel('$\Im(s)$', 'interpreter', 'latex', 'fontsize', 15)
% zlabel('$\Re(s)$', 'interpreter', 'latex', 'fontsize', 15')
% %zlim([0 2.5])
% %ylim([-2 2])
% set(gcf, 'Color', 'w');
% box on
%%
stable_vec=zeros(size(gp,2),size(gd,2));h2bound_vec=stable_vec;stable_vec=stable_vec;
for i=1:length(gp)
    for j=1:length(gd)
   % stable_vec(i)=sys_siso{i}.stable;
    norm2_matlab(i,j)=sys_siso{i,j}.norm2_matlab;
    if lmi_feas(i,j)==1
        h2bound_vec(i,j)=sys_siso{i,j}.h2bound;
    end
    end
end
figure
aread=zeros(size(gp));
norm2d=aread;
for i=1:length(gp)
    aread(i)=area(i,i);
    norm2d(i)=norm2_matlab(i,i);
end
    scatter(gp,abs(aread),'k','LineWidth',1.5)
    hold on
    scatter(gp,aread,'b','LineWidth',1.5)
%scatter(gp,sqrt(abs(area_rms)),'c','LineWidth',1.5)
scatter(gp,norm2d,'r','LineWidth',1.5);
set(gcf, 'Color', 'w');
box on
grid on
xlabel('$g_p$=$g_d$', 'interpreter', 'latex', 'fontsize', 15')
legend('$|area|$','$area$','$H_2$ norm','interpreter', 'latex')
set(gca,'TickLabelInterpreter','Latex','fontsize',15);
%%
% stable_ind=f0(stable_vec);
% disp('values of gain for which the structure is stable according to the LMI')
% disp(gp(stable_ind))
% stableind=find(norm2_matlab~=Inf);
% gp(stableind(1))
% gp(stableind(end))

%echo off