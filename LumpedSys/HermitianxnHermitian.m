%%
k=10; m1=1;gp=-0.3;
gd=0;gi=0;gdd=0;%-11;
b=0.01;
% angular frequency vector limit
flim=3.5e0;
wv=2*pi*(0.01:0.002:flim);
fv=wv/2/pi;
    %Omega=omega_a1./wn;
       muvec=pi*linspace(-1,1,701);
    omega_a1=zeros(size(muvec));omega_a2=omega_a1;
    for i=1:701
        cmu=1-cos(muvec(i));
        emu=1-exp(1i*muvec(i));
        a1=1;
        a2=0;
        a3=gp/m1*emu-2*k/m1*cmu;
        omega_r_1(:,i)=roots([a1 a2 a3]);
        a2=-1i*2*b/m1*cmu;
        a3=-2*k/m1*cmu;
        omega_r_2(:,i)=roots([a1 a2 a3]);
        a2=0;
        omega_r_3(:,i)=roots([a1 a2 a3]);
    end
     wn=sqrt(k/m1);
     Omega_r_passive=omega_r_3/wn;
      Omega_r_damp=omega_r_2/wn;
       Omega_r=omega_r_1/wn;
       %%
figure
    plot(muvec(1:10:end)/pi,real(Omega_r_passive(2,1:10:end)),'ok','LineWidth',2)
    hold on
 %   plot(muvec(1:10:end)/pi,real(Omega_r_passive(2,1:10:end)),'ok','LineWidth',1.5)
        plot(muvec/pi,real(Omega_r_damp(1,:)),'k','LineWidth',2)
   %     plot(muvec/pi,real(Omega_r_damp(2,:)),'k','LineWidth',1.5)
       plot(muvec/pi,real(Omega_r(1,:)),'r','LineWidth',2)
 %      plot(muvec/pi,real(Omega_r(2,:)),'r','LineWidth',1.5)
    ylabel('$\Omega_R$','interpreter', 'latex', 'fontsize', 20)
    xlabel('$\mu/\pi$','interpreter', 'latex', 'fontsize', 20)
    box on
    grid on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',25);
    legend('$b=0, \gamma^d=0$','$b=0.01, \gamma^d=0$','$b=0, \gamma^d=-0.03$','interpreter', 'latex', 'fontsize', 20)
    export_fig Disp_Real.pdf
    %%
     figure
     plot(muvec(1:10:end)/pi,imag(Omega_r_passive(1,1:10:end)),'ok','LineWidth',2)
     hold on
     plot(muvec/pi,imag(Omega_r_damp(1,:)),'k','LineWidth',2)
    plot(muvec/pi,imag(Omega_r(1,:)),'r','LineWidth',2)
    
    %hold on
    %plot(imag(omega_a2),muvec,'r','LineWidth',1.5)
    ylabel('$\Omega_I$','interpreter', 'latex', 'fontsize', 20)
    xlabel('$\mu/\pi$','interpreter', 'latex', 'fontsize', 20)
    box on
    grid on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',20);
        legend('$b=0, \gamma^d=0$','$b=0.01, \gamma^d=0$','$b=0, \gamma^d=-0.03$','interpreter', 'latex', 'fontsize', 20)
 export_fig Disp_Imag_1.pdf
         figure
     plot(muvec(1:10:end)/pi,imag(Omega_r_passive(2,1:10:end)),'ok','LineWidth',2)
     hold on
     plot(muvec/pi,imag(Omega_r_damp(2,:)),'k','LineWidth',2)
    plot(muvec/pi,imag(Omega_r(2,:)),'r','LineWidth',2)
    %hold on
    %plot(imag(omega_a2),muvec,'r','LineWidth',1.5)
    ylabel('$\Omega_I$','interpreter', 'latex', 'fontsize', 20)
    xlabel('$\mu/\pi$','interpreter', 'latex', 'fontsize', 20)
    box on
    grid on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',20);
        legend('$b=0, \gamma^d=0$','$b=0.01, \gamma^d=0$','$b=0, \gamma^d=-0.03$','interpreter', 'latex', 'fontsize', 20)
 export_fig Disp_Imag_2.pdf
    %%
         figure
    plot(real(Omega_r(1,:)),imag(Omega_r(1,:)),'r','LineWidth',2)
    hold on
    plot(real(Omega_r(2,:)),imag(Omega_r(2,:)),'r','LineWidth',2)
    plot(real(Omega_r_damp(1,:)),imag(Omega_r_damp(1,:)),'k','LineWidth',2)
    plot(real(Omega_r_damp(2,:)),imag(Omega_r_damp(2,:)),'k','LineWidth',2)
    plot(real(Omega_r_passive(1,:)),imag(Omega_r_passive(1,:)),'ok','LineWidth',2)
    plot(real(Omega_r_passive(2,:)),imag(Omega_r_passive(2,:)),'ok','LineWidth',2)
    xlabel('$\Omega_R$','interpreter', 'latex', 'fontsize', 20)
    ylabel('$\Omega_I$','interpreter', 'latex', 'fontsize', 20)
    box on
    grid on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',20);
        legend('$b=0, \gamma^d=0$','$b=0.01, \gamma^d=0$','$b=0, \gamma^d=-0.03$','interpreter', 'latex', 'fontsize', 20)
 export_fig Disp_Plane.pdf