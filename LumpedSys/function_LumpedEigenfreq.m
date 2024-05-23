function output= function_LumpedEigenfreq(ncell,m1,m2,k,eta,flim,boundary,e,gamma_c,theta,phi)
s=tf('s');   
%% FEM
    output0 = function_buildFEM_xl(m1,m2,k, eta ,ncell,flim,'boundary',boundary);
    %M=output.M;C=output.C;K=output.K;
    Ms=output0.Ms;Cs=output0.Cs;Ks=output0.Ks;
    output.x=output0.x;output.ndof=output0.ndof;%% State Space of interconnection (S,H) 
    %% S
    % number of inputs=outputs
    output0=function_buildSS_xl(Ms,Cs,Ks,e);
    Ass=output0.Ass;nx=size(Ass,1);
    % PS: (Ms,Cs,Ks) were previously multiplied by B*rho,
    % so we need to multiply matrixes F (inside B1ss) and T(inside B2ss) also 
    B1ss=output0.B1ss;nw=size(B1ss,2);
    B2ss=output0.B2ss;nu=size(B2ss,2);
    C1ss=output0.C1ss;nz=size(C1ss,1);
    C2ss=output0.C2ss;ny=size(C2ss,1);
    D11ss=output0.D11ss;
    D12ss=output0.D12ss;
    D21ss=output0.D21ss;

    %% H 
    Nd=100*2*pi*flim;
    % Feedback law in pressure by volume velocity
    H=gamma_c(1)+gamma_c(2)/s+gamma_c(3)*Nd/(1+Nd/s);
    %%  closed loop (S,H)
    if  norm(gamma_c,1)~=0
    %     figure
    %     bode(H)
    %     title('bode diagram of the feedback law')
        [num,dem]=tfdata(H,'v');
        if size(num,1) > size(dem,1)
            disp('improper feedback law\n');
        end
        %% controlable state space realization of H
        [Ac1,Bc1,Cc1,Dc1]=tf2ss(num,dem);
        nxc=size(Ac1,1);
        %nu=ny
        Acss=zeros(ny*nxc,ny*nxc);
        Bcss=zeros(ny*nxc,ny);
        Ccss=zeros(ny,ny*nxc);
        Dcss=zeros(ny,ny);
        %% building state space realization of H (Acss,Bcss,Ccss,Dcss)
        % by blocks, since each cell has a decoupled feedback law
        for k=1:ny
            lines=(k-1)*nxc+1:k*nxc;
            Acss(lines,lines)=Ac1;
            Bcss(lines,k)=Bc1;
            Ccss(k,lines)=Cc1.*cos(theta*k+phi);
            Dcss(k,k)=Dc1.*cos(theta*k+phi);
        end
        %Acss=Ac1*eye(ncell);   Bcss=Bc1*eye(ncell);
        %Ccss=Cc1*eye(ncell);  
        %Dcss=Dc1*eye(ny);
        %% closed loop (S,H) interconnection
        Ass_cl=[Ass+B2ss*Dcss*C2ss B2ss*Ccss
                Bcss*C2ss           Acss];
        Bss_cl=[B1ss+B2ss*Dcss*D21ss
                    Bcss*D21ss];
        Css_cl=[C1ss+D12ss*Dcss*C2ss     D12ss*Ccss];
        Dss_cl=D11ss+D12ss*Dcss*D21ss;
        output.sys_wz = ss(Ass_cl,Bss_cl,Css_cl,Dss_cl);
        Css=Css_cl;
        output.nx_cl=size(Ass_cl,1);
    else 
        % only logical case other then the previous is gamma_c==0, meaning open
        % loop system, u =0
        output.sys_wz = ss(Ass,B1ss,C1ss,D11ss);
        Css=C1ss;
    end
    %% verify internal stability of closed loop system
    if  norm(gamma_c,1)~=0
        [V, Lambda] = eig(Ass_cl);
    else
        [V, Lambda] = eig(Ass);
    end
    output.Lambda = diag(Lambda);
    output.Vn = Css*V;
    output.wn = -1i*Lambda;
 end