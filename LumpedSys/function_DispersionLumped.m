%% Dispersion of lumped mass system with feedback
% inputs:
% k=stiffness factor, m= mass. b=viscous damping factor such that m \ddot x
% + b \dot x+ k x = f
% r=range, in cells, between actuator and sensor (r=0 for local feedback)
% gp=proportional gain (feedback in terms of relative displacement);
% gd=derivative gain (feedback in terms of relative velocite);
% gdd=double derivative gain (feedback in terms of relative acceleration);
% g=integral gain (dynamic feedback);
% muvec-wavenumber vector Ex: muvec=pi*linspace(-1,1,1000)
%output:
% omega_r= dispersion matrtiz \omega(\mu) \in \C^{3 x length(muvec)}
function omega_r=function_DispersionLumped(k,m,b,r,gp,gd,gi,gdd,muvec)
        %omega_a1=zeros(length(gp),length(muvec));omega_a2=omega_r;
     for j=1:length(muvec)
            cmu=1-cos(muvec(j));
            emu=1-exp(1i*muvec(j));
            a1=(1-gdd/m*exp(1i*muvec(j)*r)*emu)*1i;
            a2=2*b/m*cmu-gd/m*exp(1i*muvec(j)*r)*emu;
            a3=(gp/m*exp(1i*muvec(j)*r)*emu-2*k/m*cmu)*1i;
            a4=gi/m*exp(1i*muvec(j)*r)*emu;
            omega_r(:,j)=roots([a1 a2 a3 a4]);
    end
end    