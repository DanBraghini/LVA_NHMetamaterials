%% clear 
clear all 
%close all 
clc 
%echo on
%% functionality flags
performance=0;
bodeplot=0;
rootlocus=0;
dispersion=1;
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
k=10; m1=1; m2=1;wn=sqrt(k/m1);
%gp=[-37:0.01:-35+0.01 -35:1:10];
%gd=[-40:1:0];
% gd=[-10000:100:-500 -500:1:-50 -50:1e-2:-1 -0.99:1e-3:0];
%gp=[0:0.1:5 5.1:0.01:8 8:1e-3:10];
%gp=[-100:.1:-2 -2:1e-1:0];
gp=-100:1:100;
gdd=0.*gp;
gd=gdd;gi=gdd;
r=1;
b=0;
e='m';
ndof=91;
for i=1:length(gp)
     muvec=pi*linspace(-1,1,1000);
        %omega_a1=zeros(length(gp),length(muvec));omega_a2=omega_r;
     for j=1:length(muvec)
             cmu=1-cos(muvec(j));
             emu=1-exp(1i*muvec(j));
            a1=(1-gdd(i)/m1*exp(1i*muvec(j)*r)*emu)*1i;
            a2=2*b/m1*cmu-gd(i)/m1*exp(1i*muvec(j)*r)*emu;
            a3=(gp(i)/m1*exp(1i*muvec(j)*r)*emu-2*k/m1*cmu)*1i;
            a4=gi(i)/m1*exp(1i*muvec(j)*r)*emu;
            omega_r{i}(:,j)=roots([a1 a2 a3 a4])./wn;
    end

end
%%
   figure
   box on
    set(gcf, 'Color', 'w');
    for i=1:length(gp)
        plot(real(omega_r{i}(1,:)),imag(omega_r{i}(1,:)),'.r',real(omega_r{i}(2,:)),imag(omega_r{i}(2,:)),'.k',real(omega_r{i}(3,:)),imag(omega_r{i}(3,:)),'.b','LineWidth',1.5)
        %hold on
       % scatter(real(omega_r{i}(2,:)),imag(omega_r{i}(2,:)),'.k','LineWidth',1.5)
       % scatter(real(omega_r{i}(3,:)),imag(omega_r{i}(3,:)),'.b','LineWidth',1.5)
       mu=zeros(2,1);
       for k=2:3
           k1=1;k2=k1;
           for j=2:length(muvec)
                       if angle(omega_r{i}(k,j))-angle(omega_r{i}(k,j-1)) >= 0
                            k1=k1+1;
                       elseif angle(omega_r{i}(k,j))-angle(omega_r{i}(k,j-1)) <= 0
                           k2=k2+1;             
                       end                 
           end
           if abs(abs(angle(omega_r{i}(k,j)))-abs(angle(omega_r{i}(k,1))))<=1e-4
               if abs(k1-length(muvec))<=2
                   mu(k-1)=1;
               elseif  abs(k2-length(muvec))<=2
                   mu(k-1)=-1;
               end
           end
       end
       txt = ['\nu =' num2str(mu(1))];
       text(-20,10,txt,'FontSize',14)
        txt = [' \nu =' num2str(mu(2))];
       text(20,10,txt,'Color','blue','FontSize',14)
       txt = ['g =' num2str(gd(i))];
       text(0,10,txt,'FontSize',14)
        omegam=cell2mat(omega_r);
        xlim([-max(max(real(omegam))) max(max(real(omegam)))]);
        ylim([-max(max(imag(omegam))) max(max(imag(omegam)))])
       % ylim([-15 15]);xlim([-15 15]);
        grid on
        %hold on
        ylabel('$\Omega_I$ ','interpreter', 'latex', 'fontsize', 15)
        xlabel('$\Omega_R$','interpreter', 'latex', 'fontsize', 15)
        zlabel('$\mu$','interpreter', 'latex', 'fontsize', 15)
        set(gca,'TickLabelInterpreter','Latex','fontsize',20);
        drawnow
        %pause(0.05)
    end
    %%
    ind=90;
figure

plot3(real(omega_r{ind}(2,:)),imag(omega_r{ind}(2,:)), muvec/pi,'r',real(omega_r{ind}(3,:)), imag(omega_r{ind}(3,:)), muvec/pi,'--r','LineWidth',2);hold on
xlabel('$\Omega_r$','interpreter', 'latex','FontSize',15')
ylabel('$\Omega_i$ ','interpreter', 'latex','FontSize',15')
zlabel('$\mu/\pi$ ','interpreter', 'latex','FontSize',15')
box on
grid on
set(gcf, 'Color', 'w');
