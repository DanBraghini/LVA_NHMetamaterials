%% clear 
clear all 
close all 
clc 
%echo on
%% functionality flags
performance=0;
bodeplot=0;
rootlocus=0;
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
 %gp=[-37:0.01:-35+0.01 -35:1:10];
 %gd=[-100:1:-21 -20:0.1:-1 -0.99:0.01:0];
 gp=[-40:0.1:10];
 %gp=-gp;
%gd=-1000:1:-900;
%with gd=-100:1:0;
gd=0*gp;
%gp=0.*gd;
%gp=-10.*ones(size(gd));
r=0;
b=0.01;
e='m';
ndof=9;
% angular frequency vector limit
flim=3.5e0;
df=2e-4;
wv=2*pi*(1e-4:df:flim);
fv=wv/2/pi;
% number of elements
     if ndof>1
         mvec=zeros(ndof,1);
         mvec(1:2:ndof,1)=m1; mvec(2:2:ndof,1)=m2;
         Ms=diag(mvec);
        % Ms=diag([1 2 3 4 5 6 7 8 9]);
         D2_=toeplitz([2 -1 zeros(1,ndof-2)]);
         D2_(1,1)=1; D2_(ndof,ndof)=1;
         Ks=k.*D2_;
        %viscous damping model
         Cs=b.*D2_;
     else
          disp('Warning: only one mass isnt enough to have feedback' )
     end
%     %fixed-fixed
    if boundary==1
        ndof_n=ndof-2;
         Ms(1,:)=[];Ms(:,1)=[];Ms(end,:)=[];Ms(:,end)=[];
         Ks(1,:)=[];Ks(:,1)=[];Ks(end,:)=[];Ks(:,end)=[];
         Cs(1,:)=[];Cs(:,1)=[];Cs(end,:)=[];Cs(:,end)=[];
     % WABC
   elseif boundary==2
       ndof_n=ndof-1;
         Ms(:,1)=Ms(:,1)+Ms(:,end);
         Ms(1,:)=Ms(1,:)+Ms(end,:);
         Ms(:,end)=[];Ms(end,:)=[]; 
         Ks(:,1)=Ks(:,1)+Ks(:,end);  
         Ks(1,:)=Ks(1,:)+Ks(end,:);
         Ks(:,end)=[];Ks(end,:)=[];
         Cs(:,1)=Cs(:,1)+Cs(:,end);  
         Cs(1,:)=Cs(1,:)+Cs(end,:);
         Cs(:,end)=[];Cs(end,:)=[];
    elseif boundary==0
        ndof_n=ndof;
    end
so=zeros(2*ndof_n,length(gp));
area=zeros(size(gp));lmi_feas=area;
q1m=zeros(length(gp),length(wv));q1p=q1m;qnm=q1m;qnp=q1m;
for i=1:length(gp)
    %gamma_c=[gp(i) 0 0];
   % output=function_CL4MSS(ndof,m1,m2,k,b,flim,boundary,e,gamma_c);
    %% Computing FRF for each gain
     % Feedback Matrix
      D1_=(diag(-1*ones(1,ndof))+diag(ones(1,ndof-1),-1));
      D1_(1,:)=zeros(1,ndof);
      Gsp=gp(i).*D1_;
      Gsd=gd(i).*D1_;
      %fixed-fixed
     if boundary==1
         Gsp(1,:)=[];Gsp(:,1)=[];Gsp(end,:)=[];Gsp(:,end)=[]; 
         Gsp(1,:)=zeros(1,ndof_n);
         Gsd(1,:)=[];Gsd(:,1)=[];Gsd(end,:)=[];Gsd(:,end)=[]; 
         Gsd(1,:)=zeros(1,ndof_n);
    %WABC
     elseif boundary==2
         Gsp(:,1)=Gsp(:,1)+Gsp(:,end);
         Gsp(1,:)=Gsp(1,:)+Gsp(end,:);
         Gsp(:,end)=[];Gsp(end,:)=[]; 
         Gsd(:,1)=Gsd(:,1)+Gsd(:,end);
         Gsd(1,:)=Gsd(1,:)+Gsd(end,:);
         Gsd(:,end)=[];Gsd(end,:)=[]; 
     end
 %% define where the external load is applied
%  % F=1 for FRF
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
    %% closed-loop state-space 4-matrix representation SIMO for FRFs
    nw=1;nz=ndof_n;
    A=[zeros(ndof_n) eye(ndof_n)
           -Ms\(Ks+Gsp)    -Ms\(Cs+Gsd)];
    B=[zeros(ndof_n,1)
             Ms\F];
    C=[eye(ndof_n)  zeros(ndof_n)];
    D=zeros(nz,nw);
    sys_simo=ss(A,B,C,D); 
    %sys_simo=output.CLsys;
    %[A,B,C,D]=ssdata(sys_simo);
    P=tf(sys_simo);
    % Taking BODE diagrams of first and last DOFs
    P1=P(1,1);Pn=P(ndof_n,1);
    % BODE diagram of tf for output q(1)
    if bodeplot==1 || performance==1
        [mag,phase]=bode(P1,wv);
        q1m(i,:)=vec(mag);
        q1p(i,:)=vec(phase);
        % BODE diagram of tf for output q(n)
        [mag,phase]=bode(Pn,wv);
        qnm(i,:)=vec(mag);
        qnp(i,:)=vec(phase);
        %dw=wv(2)-wv(1);
        area(i)=2*(trapz(qnm(i,:))-trapz(q1m(i,:)))*df;
        sys_siso{i}.area=area(i);
        %% state-space 4-matrix representation SISO for H2 norm
        nw=1;nz=1;
        C=[1 zeros(1,ndof_n-2) -1 zeros(1,size(A,2)-ndof_n)];
        D=zeros(nz,nw);
        sys_siso{i}.statespace=ss(A,B,C,D);
        sys_siso{i}.norm2_matlab=norm(sys_siso{i}.statespace,2);
    end

    %% Computing Eigenspectum for each gain (free-vibration problem)
    % s=-i \omega
    %((-i\omega)^2 M -i\omega C +K+G)q(\omega) = 0
%     if b ~=0
%         [V,so(i,:)]=polyeig(Ks+Gs,Cs,Ms);
%     else
%         [V,L]=eig(Ks+Gs,-Ms);
%         so(i,:)=[sqrt(diag(L));-sqrt(diag(L))];
%     end
%so{i}=output.eig;
[V,D]=eig(A);
E{i}.V=V;
so(:,i)=diag(D);
E{i}.r=so(:,i);
end

%% Bode plot
if bodeplot==1
    if gp(2)~=gp(1)
        g=gp;
    else
        g=gd;
    end
    f0= find(g==-11);
    qN0m=qnm(f0,:);
    q10m=q1m(f0,:);
    q10p=q1p(f0,:);
    qN0p=qnp(f0,:);
    figure
    plot(fv,20*log10(q10m),'r','MarkerSize',8,'LineWidth',1)
    hold on
    plot(fv,20*log10(qN0m),'k','MarkerSize',8,'LineWidth',1)
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
    figure
    plot(fv,q10p,'r','MarkerSize',8,'LineWidth',1)
    hold on
    plot(fv,qN0p,'k','MarkerSize',8,'LineWidth',1)
    hold off
    ylabel('$\angle \hat{q} / \hat{F} (i \omega) $', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$f$ [Hz]', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    legend('L','R')
    xlim([0 flim])
    box on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    grid on
end
%%
if performance==1
for i=1:length(gp)
    norm2_matlab(i)=sys_siso{i}.norm2_matlab;
end
if gp(2)~=gp(1)
    g=gp;
else
    g=gd;
end
    figure
    scatter(g,abs(area),'k','LineWidth',1.5)
    hold on
    scatter(g,norm2_matlab,'r','LineWidth',1.5);
    scatter(g,max(real(so)),'b','LineWidth',1.5);
    set(gcf, 'Color', 'w');
    box on
    grid on
    xlabel('$g$', 'interpreter', 'latex', 'fontsize', 15')
    %legend('$|area|$','$area$','$H_2$ norm','interpreter', 'latex')
    legend('$Aclo$','$H_2$ norm','$\Re (\lambda_{max})$','interpreter', 'latex')
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
end
%% Root locus
if rootlocus==1
    if gp(2)~=gp(1)
        g=gp;
    else
        g=gd;
    end
    figure
    ind=find(g==0);
    scatter(real(so(:,ind)),imag(so(:,ind)),'mx','LineWidth',2,'SizeData',100)
    hold on
    son = function_SortMAC(so);
   %[son,~] = sort(so,1,'ComparisonMethod','abs');
    for i=1:size(so,1)
        scatter(real(so(i,:)),imag(so(i,:)),'.k','LineWidth',2);
        %plot(real(so(i,:)),imag(so(i,:)),'LineWidth',2)
    end
    grid on
    title('Root Locus','interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Re(s)$', 'interpreter', 'latex', 'fontsize', 15)
    ylabel('$\Im(s)$', 'interpreter', 'latex', 'fontsize', 15')
    set(gcf, 'Color', 'w');
    box on
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    hold off
end
%echo off
%%

