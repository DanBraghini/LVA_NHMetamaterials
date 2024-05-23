%% clear 
clear all 
%close all 
clc 
%% lumped-parameter system set up
k=10; m1=1; m2=m1;
%k2=1; m2=2;
%r1=0; r2=0; %range of non-reciprocal coupling, r=0 means next-neighborhoud
%couple
r=0;
eta=0;
e='m';
% boundary == 0 : free-free (default)
% bundary == 1: fixed-fixed
% boundary == 2 : PBC
boundary=1;
% angular frequency vector limit
flim=3.5e3;
wv=2*pi*(2:2:flim);
fv=wv/2/pi;
Nd=100*2*pi*flim;
% gains
kp=0.6;
ki=kp;
kd=kp;
s=tf('s');
modulation = 'alpha';
ncell=10;
% qusiperiodicity parameters
if modulation == 'phiii'
    theta=1/3;phi=(0:1:ncell*2*pi)./ncell;
    parm=phi;
elseif modulation == 'theta'
     phi=0;%theta=0:.1:2*pi;
     theta=(0:1:ncell*2*pi)./ncell;
     parm=theta;
elseif modulation == 'alpha'
     phi=0;theta=0;
%     kp=-1e2:10:1e2;
     kp=linspace(-100,100,100);
%      for i=1:length(kp)
%          if kp(i)==0
%              kp(i) = kp(i-1);
%          end
%      end
     parm=kp;
     kp=kp';
end

gamma_c=[1*kp (0*ki).*ones(size(kp)) (0*kd).*ones(size(kp))];
% Feedback law in pressure by volume velocity
H_v=gamma_c(:,1)+gamma_c(:,2)/s+gamma_c(:,3)*Nd/(1+Nd/s);
%% FEM
output = function_buildFEM_xl(m1,m2,k, eta ,ncell,flim,'boundary',boundary);
%M=output.M;C=output.C;K=output.K;
Ms=output.Ms;Cs=output.Cs;Ks=output.Ks;
x=output.x;ndof=output.ndof;%% State Space of interconnection (S,H) 
%% S
output=function_buildSS_xl(Ms,Cs,Ks,e);
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


%% setting analog filters
% RC filter
% wf=1e3*2*pi;
% Hf = s/wf/(1+s/wf);
% no filter
Hf=1;
%butterworth filter
% w1=1800*2*pi;
% w2=2000*2*pi;
% % lp= 'low pass', bp = 'band pass'
% type='hp';
% out=function_filter(-3,-10,w1,w2,type);
% Hf=out.Hf;
%Hf=Filtro_ver2_completo_v1;
% final feedback transfer function
H=H_v*Hf;
%%
np=length(parm);
%wo=zeros(nxcl,np);
%Vo=zeros(length(x),nxcl,np);
for i=1:np
    Hm=H;
   if modulation == 'phiii'
       phim=phi(i);thetam=theta;
   elseif modulation == 'theta'
       phim=phi;thetam=theta(i);
   else
       phim=0;thetam=0;Hm=H(i);
   end
%    if norm(Hm) == 0
%         continue
%    end
   output = function_eigenspectrum_varyingp(Hm, Ass,B2ss,C1ss,C2ss,D12ss,thetam,phim);
   wo(:,i)=output.wo;
   Vo(:,:,i)=output.Vo;
end
%%
nx_cl=output.nx_cl;
% figure
% for i=1:np
%     plot(parm(i)/pi.*ones(1,nx_cl),real(wo(:,i)),'.k')
%     hold on
% end
% %plot(real(wo(1:nx,:)),phi/pi)
% grid on
% title('Modulated spectrum', 'interpreter', 'latex', 'fontsize', 15)
% ylabel('$\Re \{ \Omega \}$', 'interpreter', 'latex', 'fontsize', 15)
% xlabel('$parameter / \pi$', 'interpreter', 'latex', 'fontsize', 15')
% set(gcf, 'Color', 'w');
% box on
% %ylim([0 2.5])
% %%
% figure
% for i=1:np
%     plot(parm(i)/pi.*ones(1,nx_cl),imag(wo(:,i)),'.k')
%     hold on
% end
% %plot(real(wo(1:nx,:)),phi/pi)
% grid on
% title('Modulated spectrum', 'interpreter', 'latex', 'fontsize', 15)
% ylabel('$\Im \{ \Omega \}$', 'interpreter', 'latex', 'fontsize', 15)
% xlabel('$parameter / \pi$', 'interpreter', 'latex', 'fontsize', 15')
% set(gcf, 'Color', 'w');
% box on
% %%
% 
% % figure
% % for i=1:np
% %     plot(imag(wo(:,i)),real(wo(:,i)),'.k')
% %     hold on
% % end
% % %plot(real(wo(1:nx,:)),phi/pi)
% % grid on
% % title('Modulated spectrum', 'interpreter', 'latex', 'fontsize', 15)
% % xlabel('$\Im \{ \Omega \}$', 'interpreter', 'latex', 'fontsize', 15)
% % ylabel('$\Re  \{ \Omega \}$', 'interpreter', 'latex', 'fontsize', 15')
% % set(gcf, 'Color', 'w');
% % box on
% %%
 figure
for i=1:np
    plot3(parm(i).*ones(1,nx_cl),imag(wo(:,i)),real(wo(:,i)),'.k')
    hold on
end
%%
% hold on
% for i=1:nx_cl
%     plot3(parm,imag(wo(i,:)),real(wo(i,:)),'r')
%     hold on
% end
%%
grid on
title('Modulated spectrum', 'interpreter', 'latex', 'fontsize', 15)
xlabel('$parameter $', 'interpreter', 'latex', 'fontsize', 15')
ylabel('$\Im \{ \Omega \}$', 'interpreter', 'latex', 'fontsize', 15)
zlabel('$\Re  \{ \Omega \}$', 'interpreter', 'latex', 'fontsize', 15')
zlim([0 2.5])
ylim([-2 2])
set(gcf, 'Color', 'w');
box on


