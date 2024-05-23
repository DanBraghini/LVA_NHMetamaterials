%% clear 
clear all 
%close all 
%clc 
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
 %gp=[-37:0.01:-35+0.01 -35:1:10]
 %gd=[-100:1:-21 -20:0.1:0];
 g=[-1000:10:-110 -100:1:-10 -9.99:0.01:0];
 %gp=-gp;
%gd=-1000:1:-900;
%gp=-100:1:0;
%gd=0*gp;
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
so=cell(size(g));
q1m=zeros(length(g),length(wv));q1p=q1m;qnm=q1m;qnp=q1m;
for i=1:length(g)
    gamma_c=[0  g(i) 0];
    output=function_CL4MSS(ndof,m1,m2,k,b,flim,boundary,e,gamma_c);
    so{i}=output.eig;
    [A,B,~,~]=ssdata(output.CLsys);
     C=[eye(ndof-2)  zeros(ndof-2,size(A,1)-ndof+2)];
    D=zeros(size(C,1),size(B,2));
    sys_simo=ss(A,B,C,D); 
    P=tf(sys_simo);
    % Taking BODE diagrams of first and last DOFs
    P1=P(1,1);Pn=P(ndof-2,1);
     [mag,phase]=bode(P1,wv);
    q1m(i,:)=vec(mag);
    q1p(i,:)=vec(phase);
    % BODE diagram of tf for output q(n)
    [mag,phase]=bode(Pn,wv);
    qnm(i,:)=vec(mag);
    qnp(i,:)=vec(phase);
end
%%
ind=find(g==0);
figure
scatter(real(so{ind}),imag(so{ind}),'mx','LineWidth',2,'SizeData',100)
hold on
so{ind}=[];
R=cell2mat(so);
[son,~] = sort(R,'ComparisonMethod','abs');
for i=1:size(R,1)
    scatter(real(R(i,:)),imag(R(i,:)),'.k','LineWidth',2);hold on;
    %plot(real(R(i,:)),imag(R(i,:)),'LineWidth',2)
end
grid on
title('Root Locus','interpreter', 'latex', 'fontsize', 15)
xlabel('$\Re(s)$', 'interpreter', 'latex', 'fontsize', 15)
ylabel('$\Im(s)$', 'interpreter', 'latex', 'fontsize', 15')
set(gcf, 'Color', 'w');
box on
set(gca,'TickLabelInterpreter','Latex','fontsize',15);
hold off
%%
  f0= find(g==-50);
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
%echo off


