function [x,ndof,wn,Vn,sys_wz,nx_cl] = function_LumpedEigenfreq(ncell,m1,m2,k,c,flim,boundary,e,gamma_c,theta,phi)
s=tf('s');   
%% FEM
    output = function_buildFEM_xl(m1,m2,k, c ,ncell,flim,'boundary',boundary);
    %M=output.M;C=output.C;K=output.K;
    Ms=output.Ms;Cs=output.Cs;Ks=output.Ks;
    x=output.x;ndof=output.ndof;%% State Space of interconnection (S,H) 
    %% S
    % number of inputs=outputs
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

    %% H 
    Nd=100*2*pi*flim;
    % PID Feedback law in pressure by volume velocity
    H_v=gamma_c(1)+gamma_c(2)/s+gamma_c(3)*Nd/(1+Nd/s);

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
        sys_wz = ss(Ass_cl,Bss_cl,Css_cl,Dss_cl);
        Css=Css_cl;
        nx_cl=size(Ass_cl,1);
    else 
        % only logical case other then the previous is gamma_c==0, meaning open
        % loop system, u =0
        sys_wz = ss(Ass,B1ss,C1ss,D11ss);
        Css=C1ss;
        nx_cl=nx;
    end
    %% verify internal stability of closed loop system
    if  norm(gamma_c,1)~=0
        [V, Lambda] = eig(Ass_cl);
    else
        [V, Lambda] = eig(Ass);
    end
    Lambda = diag(Lambda);
    Vn = Css*V;
    wn=Lambda;
    %wn = -1i*Lambda;
    [wn,ind] = sort(wn,'ComparisonMethod','real');
    Vn = Vn(:,ind); 

 end