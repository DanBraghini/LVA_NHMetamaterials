function output = function_CL4MSS(ncell,m1,m2,k,c,flim,boundary,e,r,derivative,gamma_c,H)
%% FEM
    outputFEM = function_buildFEM_xl(m1,m2,k, c ,ncell,flim,'boundary',boundary);
    %M=output.M;C=output.C;K=output.K;
    Ms=outputFEM.Ms;Cs=outputFEM.Cs;Ks=outputFEM.Ks;
    x=outputFEM.x;N=size(Ms,1);
%% State Space of interconnection (S,H) 
    %% 4-Matrix State Space representation of passive system S
    outputSS=function_buildSS_xl(Ms,Cs,Ks,e,r,derivative);
    Ass=outputSS.Ass;nx=size(Ass,1);
    % PS: (Ms,Cs,Ks) were previously multiplied by B*rho,
    % so we need to multiply matrixes F (inside B1ss) and T(inside B2ss) also 
    B1ss=outputSS.B1ss;nw=size(B1ss,2);
    B2ss=outputSS.B2ss;nu=size(B2ss,2);
    C1ss=outputSS.C1ss;nz=size(C1ss,1);
    C2ss=outputSS.C2ss;ny=size(C2ss,1);
    D11ss=outputSS.D11ss;
    D12ss=outputSS.D12ss;
    D21ss=outputSS.D21ss;
    D22ss=outputSS.D22ss;
    %% 4-Matrix reresentation Passive System
    sys_open=ss(Ass,B1ss,C1ss,D11ss);
    %% 4-Matriz State Space representation of Feedback system H
    %H=gamma_c(1)/(s+100)^2;
    if  norm(gamma_c,1)~=0
            if derivative==1.1
                %nu=2ny
                    Acss=zeros(0,0);
                    Bcss=zeros(0,ny);
                    Ccss=zeros(nu,0);
                    Dcss=[gamma_c(1)*eye(N-1) gamma_c(2)*eye(N-1)];
             else

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
            end
            %% closed loop (S,H) interconnection
            Ass_cl=[Ass+B2ss*Dcss*C2ss B2ss*Ccss
                    Bcss*C2ss           Acss];
            Bss_cl=[B1ss+B2ss*Dcss*D21ss
                        Bcss*D21ss];
            Css_cl=[C1ss+D12ss*Dcss*C2ss     D12ss*Ccss];
            Dss_cl=D11ss+D12ss*Dcss*D21ss;
            sys_wz = ss(Ass_cl,Bss_cl,Css_cl,Dss_cl);
    else 
        % only logical case other then the previous one is gamma_c==0, meaning open
        % loop system, u =0
        sys_wz = ss(Ass,B1ss,C1ss,D11ss);
    end
   %% Evaluating outputs
   %Closed loop 4-matrix state-space representation
   output.CLsys=sys_wz;
   output.OLsys=sys_open;
   output.FEMmesh=x;
   output.ndof_n=outputFEM.ndof;
   output.nx=nx;
 end