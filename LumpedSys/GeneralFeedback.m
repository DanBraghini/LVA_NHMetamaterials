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
% Laplace variable to define feedback law
s=tf('s');
% flag used for feedback in terms of velocity rather than displacement
%derivative=-1: integral feedback (position)
%derivative=0: position static feedback
%derivative=1: velocity static feedback
%derivative=2: acceleration static feedback
%derivative=1.1: position and velocity static feedback
derivative=2;
eigenmodes=1;
performance=0;
stability=0;
%% lumped-parameter system set up
% m = mass 
% k = string elastic coefficient
% range of non-reciprocal coupling, r=0 means next-neighborhoud couple
% e = excitation point for the external 
k=10; m1=1; m2=1;wn=sqrt(k/m1);
%g =[-100:1:-10 -9.9:0.1:0];
%g=[-20:1:-1 -0.9:0.1:0.9 1:20];
%g=-100:.1:100;gplot=-1*k
g=.1;
gplot=g;gn=g/k;
%g=-11;
%g=[-1:.001:1];
%g=[-0.2:1e-5:0];
% g=[-1000:10:-110 -100:1:-10 -9.99:0.01:0];
%g =-1.*[-1:1e-4:0];
%g=-100:1:10;
r=0;
b=0;%0.01;
e='m';
ndof=99;%9;
% angular frequency vector limit
flim=3.5e0;
%flim=10*flim;
%flim=500;
% frequency step, used to compute the performance metric
df=2e-4;
% angular frequency vector limit
wv=2*pi*(df:df:flim);
so=cell(size(g));
q1m=zeros(length(g),length(wv));q1p=q1m;qnm=q1m;qnp=q1m;
S=zeros(size(g));Norm2=S;Norm1=S;Norminf=S;
for i=1:length(g)
    %% Define the feedback law
    if derivative==-1 % integral case
        gamma_c=[0 0 g(i)];
    elseif derivative==1.1
        gamma_c=[g(i) g(i) 0];
    else
         gamma_c=[g(i) 0 0];
    end
    H=gamma_c(1)+gamma_c(2)+gamma_c(3)*s;
    %%
    output=function_CL4MSS(ndof,m1,m2,k,b,flim,boundary,e,r,derivative,gamma_c,H);
   %output=function_CL4MSS(ndof,m1,m2,k,b,flim,boundary,'l',derivative,gamma_c,H);
   %% State Space 4-matriz representation of closed loop system
    [A,B,C,D]=ssdata(output.CLsys);
    [V, Lambda] = eig(A);
    Lambdac{i} = diag(Lambda);
    Vn{i} = C*V;
    so{i}=Lambdac{i};
    ndof_n=output.ndof_n;
    xd=1:1:ndof_n;    
    nx_cl=size(A,1);
    nxc=nx_cl-output.nx;
    if performance==1 && max(real(Lambdac{i}))<=0
        %% state-space 4-matrix representation SIMO system for Performance Metric
%         for j=1:ndof_n
%             C2(j,:)=[zeros(1,ndof_n-(ndof_n-j)-1) 1 zeros(1,ndof_n-j)  zeros(1,ndof_n)];
%             sys_siso{j}=ss(A,B,C2(j,:),0);
%             [mag(:,j),~]=bode(sys_siso{j},wv);
%         end
        C1=[1 zeros(1,ndof_n-1)  zeros(1,ndof_n)];
        CN=[zeros(1,ndof_n-1) 1  zeros(1,ndof_n)];
                sys_1=ss(A,B,C1,0);sys_N=ss(A,B,CN,0);
        % Taking BODE diagrams of first and last DOFs
        [mag,~]=bode(sys_N,wv);qnm(i,:)=mag;
                [mag,~]=bode(sys_1,wv);q1m(i,:)=mag;
    %    qnp(i,:)=vec(phase);
      %   output=function_CL4MSS(ndof,m1,m2,k,b,flim,boundary,'r',derivative,gamma_c,H);
    %  [A,B,~,~]=ssdata(output.CLsys);
    %  sys_simo=ss(A,B,C,D); 
    %  P=tf(sys_simo);
    %    P1=P(1,1);
        % BODE diagram of tf for output q(n)
      %  q1p(i,:)=vec(phase);
        S(i)=2*(trapz(qnm(i,:))-trapz(q1m(i,:)))*2*pi*df;
        %% state-space 4-matrix representation SISO system for H2 norm
         nw=1;nz=1;
         C2=[1 zeros(1,ndof_n-2) -1 zeros(1,size(A,2)-ndof_n)];
         D=zeros(nz,nw);
         sys_siso=ss(A,B,C2,D);
         Norm2(i)=norm(sys_siso,2);
         Norminf(i)=norm(sys_siso,inf);
          [mag,~]=bode(sys_siso,wv);
         Q=vec(mag);
         Norm1(i)=abs(2*trapz(Q)*2*pi*df);
    end
end
%%
if stability==1
    ind=find(g==0);
    figure
    scatter(real(so{ind}),imag(so{ind}),'mx','LineWidth',2,'SizeData',100)
    hold on
    s1=so;
    s1{ind}=[];
    g1=gn;g1(ind)=[];
    R=cell2mat(so)';
    %[son,~] = sort(so,'ComparisonMethod','abs');
for i=1:size(R,1)
        if i~=ind
            scatter(real(R(i,:)),imag(R(i,:)),'.k','LineWidth',2);  
        end
end
    grid on
     xlim([-0.03 0.01])
    title('Root Locus','interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Re(s)$', 'interpreter', 'latex', 'fontsize', 15)
    ylabel('$\Im(s)$', 'interpreter', 'latex', 'fontsize', 15')
    set(gcf, 'Color', 'w');
    box on
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    hold off
    figure
    plot(gn,max(real(R')),'b','LineWidth',1.5)
    hold on;
    scatter(0,max(real(so{ind})),'mx','LineWidth',2,'SizeData',100)
    set(gcf, 'Color', 'w');
    box on
    grid on
    xlabel('$g / k$', 'interpreter', 'latex', 'fontsize', 15')
    ylabel('$\Re( \lambda)_{max}$','interpreter', 'latex')
    set(gca,'TickLabelInterpreter','Latex','fontsize',15)
end
%%
if performance==1
%    S_=zeros(size(S));Norm2_=S_;
%     for i=1:length(g)-1
%         dg=g(i+1)-g(i);
%           S_(i)=(S(i+1)-S(i))/dg;
%           Norm2_(i)=(Norm2(i+1)-Norm2(i))/dg;
%     end
%    S_(end)=(S(end)-S(length(g)-1))/dg;
 %   Norm2_(end)=(Norm2(end)-Norm2(length(g)-1))/dg;
     figure
     scatter(gn,abs(S),'k','LineWidth',1.5)
     hold on
     scatter(gn,Norm2,'r','LineWidth',1.5);
     scatter(gn,Norm1,'b','LineWidth',1.5);
     scatter(gn,Norminf/10,'m','LineWidth',1.5);
     ind=find(Norminf==0);
     aux1=find(g(ind)<0);
     aux2=find(g(ind)>0);
     if isempty(aux1)
         g0=g(1)/k;
     else
         g0=g(ind(aux1(end))+1)/k;
     end
     if isempty(aux2)
         g1=g(end)/k;
     else
         g1=g(ind(aux2(1))-1)/k;
     end
     xlim([g0 g1])
     set(gcf, 'Color', 'w');
     box on
     grid on
     xlabel('$g / k$', 'interpreter', 'latex', 'fontsize', 15')
     legend('$S$','$H_2$','$H_1$','$H_{\infty}/10$','interpreter', 'latex')
     set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    ylim([0 20])
    %xlim([-30 7])
    %ylim([-.001 .001])
    %%
%    gplot=-0.01*k;
      g0= find(abs(g-gplot)<=1e-4);
     if max(real(Lambdac{g0}))<=0
            qN0m=qnm(g0,:);
            q10m=q1m(g0,:);
            q10p=q1p(g0,:);
            qN0p=qnp(g0,:);
            figure
            fv=wv/2/pi;
            Omegav=wv/wn;
           % plot(fv(1:100:end),20*log10(q10m(1:100:end)),'r*','MarkerSize',8,'LineWidth',1)
            plot(Omegav,20*log10(q10m),'r','MarkerSize',8,'LineWidth',1)
                hold on
            plot(Omegav,20*log10(qN0m),'k','MarkerSize',8,'LineWidth',1)
            hold off
            ylabel('$ |\hat{q} / \hat{F}| (i \omega)$ [dB]', 'interpreter', 'latex', 'fontsize', 15)
            xlabel('$\Omega$', 'interpreter', 'latex', 'fontsize', 15)
            set(gca,'TickLabelInterpreter','Latex','fontsize',15);
            legend('L','R')
           % xlim([0, max(real(Omega_n))])
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
            xlabel('$\Omega$', 'interpreter', 'latex', 'fontsize', 15)
            set(gca,'TickLabelInterpreter','Latex','fontsize',15);
            legend('L','R')
            xlim([0 flim])
            box on
            set(gcf, 'Color', 'w');
            set(gca,'TickLabelInterpreter','Latex','fontsize',15);
            grid on
    end
%echo off
end
%%
 g0= find(abs(g-gplot)<=1e-4);
omega_n=-so{g0}.*1i;
[omega_n,ord]=sort(omega_n,'ComparisonMethod','real');Omega_n=omega_n/wn;
Vnv=Vn{g0}(:,ord);
% Dispersion
muvec=pi*linspace(-1,1,1000);
if derivative==-1
    omega_r=function_DispersionLumped(k,m1,b,r,0,0,g(g0),0,muvec);
elseif derivative==0
     omega_r=function_DispersionLumped(k,m1,b,r,g(g0),0,0,0,muvec);
elseif derivative==1
         omega_r=function_DispersionLumped(k,m1,b,r,0,g(g0),0,0,muvec);
elseif derivative==2
         omega_r=function_DispersionLumped(k,m1,b,r,0,0,0,g(g0),muvec);
elseif derivative==01
         omega_r=function_DispersionLumped(k,m1,b,r,g(g0),g(g0),0,0,muvec);
end
Omega_r=omega_r/wn;
%%
figure
plot(real(Omega_r(1,:)),imag(Omega_r(1,:)),'k',real(Omega_r(2,:)),imag(Omega_r(2,:)),'r',real(Omega_r(3,:)), imag(Omega_r(3,:)),'--r','LineWidth',2);hold on
scatter(real(Omega_n),imag(Omega_n),'ko','LineWidth',1)
xlabel('$\Omega_R$','interpreter', 'latex','FontSize',15)
ylabel('$\Omega_I$','interpreter', 'latex','FontSize',15)
box on
grid on
set(gcf, 'Color', 'w');
%% Eigenmodes Animation
if eigenmodes==1
%   eigenfrequency search
    fnv=Omega_n;
    [fi_real,fi_imag]=ginput(1)
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
    disp('Omega_n')
   % ind=105;
    fn=fnv(ind);
    display(fn)
    %modo=sin(x)+1j*cos(x)/2;
    mode=Vnv(:,ind);
    norm = max(abs(mode));
    % numero de frames
    Nf=70;
    % numero de ciclos
    nc=1; 
    teta=0:2*nc*pi/(Nf-1):2*nc*pi; 
    figure
    % l=length(x);
    % x = [x x(l).*ones(1,l-1)+x(2:l)];
    for i=1:Nf
        y=abs(mode).*cos(angle(mode)+teta(i));
        plot(xd,y'/norm)
    %    y=[y;y];
    %    plot(x(1:end-1),y'/norm)
        axis([min(xd),max(xd), -1,1]);
        xlabel('$x$ [m]', 'interpreter', 'latex', 'fontsize', 15)
      %  ylabel('$ \hat{P}(x,\omega_n)/ ||\hat{P}(x,\omega_n)||_{\infty}$', 'interpreter', 'latex', 'fontsize', 15)     
        set(gca,'TickLabelInterpreter','Latex','fontsize',15);
        box on
        set(gcf, 'Color', 'w');
        legend(['\Omega_n = ' num2str(fn,'%0.3f')]);
        drawnow
    end
    %%
        figure
        y=abs(mode).*sin(angle(mode));
        %y=imag(mode);
        plot(xd,y'/norm,'.-k','LineWidth',1)
    %    y=[y;y];
    %    plot(x(1:end-1),y'/norm)
        axis([min(xd),max(xd), -1,1]);
        xlabel('$n$', 'interpreter', 'latex')
        ylabel('$ \hat{u}$', 'interpreter', 'latex')     
        set(gca,'TickLabelInterpreter','Latex','fontsize',20);
        box on
        set(gcf, 'Color', 'w');
        legend(['\Omega_n = ' num2str(fn,'%0.3f') ]);
end
%%
figure
plot3(real(Omega_r(1,:)),imag(Omega_r(1,:)), muvec/pi,'k',real(Omega_r(2,:)),imag(Omega_r(2,:)), muvec/pi,'r',real(Omega_r(3,:)), imag(Omega_r(3,:)), muvec/pi,'--r','LineWidth',2);hold on
xlabel('$\Omega_R$','interpreter', 'latex','FontSize',15')
ylabel('$\Omega_I$ ','interpreter', 'latex','FontSize',15')
zlabel('$\mu/\pi$ ','interpreter', 'latex','FontSize',15')
box on
grid on
set(gcf, 'Color', 'w');

%%
figure
plot(muvec/pi,real(Omega_r(1,:)),'k',muvec/pi,real(Omega_r(2,:)),'r',muvec/pi,real(Omega_r(3,:)),'--r','LineWidth',2);hold on
ylabel('$\Omega_R$','interpreter', 'latex','FontSize',15)
xlabel('$\mu/\pi$ ','interpreter', 'latex','FontSize',15)
box on
grid on
set(gcf, 'Color', 'w');
%%
figure
plot(muvec/pi,imag(Omega_r(1,:)),'k',muvec/pi,imag(Omega_r(2,:)),'r',muvec/pi,imag(Omega_r(3,:)),'--r','LineWidth',2);hold on
ylabel('$\Omega_I$','interpreter', 'latex','FontSize',15)
xlabel('$\mu/\pi$ ','interpreter', 'latex','FontSize',15)
box on
grid on
set(gcf, 'Color', 'w');
