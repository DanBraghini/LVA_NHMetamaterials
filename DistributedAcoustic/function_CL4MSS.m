function output = function_CL4MSS(ncell,m1,m2,k,c,flim,boundary,e,gamma_c)
s=tf('s');   
%% FEM
    output1 = function_buildFEM_xl(m1,m2,k, c ,ncell,flim,'boundary',boundary);
    %M=output.M;C=output.C;K=output.K;
    Ms=output1.Ms;Cs=output1.Cs;Ks=output1.Ks;
    x=output1.x;ndof=output1.ndof;%% State Space of interconnection (S,H) 
    %% 4-Matriz State Space representation of passive system S
    output1=function_buildSS_xl(Ms,Cs,Ks,e);
    Ass=output1.Ass;nx=size(Ass,1);
    % PS: (Ms,Cs,Ks) were previously multiplied by B*rho,
    % so we need to multiply matrixes F (inside B1ss) and T(inside B2ss) also 
    B1ss=output1.B1ss;nw=size(B1ss,2);
    B2ss=output1.B2ss;nu=size(B2ss,2);
    C1ss=output1.C1ss;nz=size(C1ss,1);
    C2ss=output1.C2ss;ny=size(C2ss,1);
    D11ss=output1.D11ss;
    D12ss=output1.D12ss;
    D21ss=output1.D21ss;

    %% 4-Matriz State Space representation of Feedback system H
    Nd=100*2*pi*flim;
    % PID Feedback law in pressure by volume velocity
    H=gamma_c(1)+gamma_c(2)/s+gamma_c(3)*Nd/(1+Nd/s);
    %H=gamma_c(1)/(s+100)^2;
    if  norm(gamma_c,1)~=0
    %     figure
    %     bode(H)
    %     title('bode diagram of the feedback law')
        [num,dem]=tfdata(H,'v');
        if size(num,1) > size(dem,1)
            disp('improper feedback law\n');
        end
        % controlable state space realization of H
        [Ac1,Bc1,Cc1,Dc1]=tf2ss(num,dem);
        nxc=size(Ac1,1);
        %nu=ny
        Acss=zeros(ny*nxc,ny*nxc);
        Bcss=zeros(ny*nxc,ny);
        Ccss=zeros(ny,ny*nxc);
        Dcss=zeros(ny,ny);
        % building state space realization of H (Acss,Bcss,Ccss,Dcss)
        % by blocks, since each cell has a decoupled feedback law
        for k=1:ny
            lines=(k-1)*nxc+1:k*nxc;
            Acss(lines,lines)=Ac1;
            Bcss(lines,k)=Bc1;
            Ccss(k,lines)=Cc1;%.*cos(theta*k+phi);
            Dcss(k,k)=Dc1;%.*cos(theta*k+phi);
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
    else 
        % only logical case other then the previous one is gamma_c==0, meaning open
        % loop system, u =0
        sys_wz = ss(Ass,B1ss,C1ss,D11ss);
        Css=C1ss;
    end
    %% verify internal stability of closed loop system
    if  norm(gamma_c,1)~=0
        [V, Lambda] = eig(Ass_cl);
    else
        [V, Lambda] = eig(Ass);
    end
    Lambda = diag(Lambda);
    Vn = Css*V;
    so=Lambda;
    %wn = -1i*Lambda;
    %[so,ind] = sort(so,'ComparisonMethod','real');
   % Vn = Vn(:,ind); 
   %% Evaluating outputs
   %Closed loop 4-matrix state-space representation
   output.CLsys=sys_wz;
   output.eig=so;
   output.eigvec=Vn;
   output.FEMmesh=x;
 end