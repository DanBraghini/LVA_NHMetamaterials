function output = function_eigenspectrum_varyingp(H, Ass,B2ss,C1ss,C2ss,D12ss,theta,phi)
    %% H 
    % Feedback law in pressure by volume acceleration
    %%  closed loop (S,H)
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
    ny=size(C2ss,1);
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
    %Bss_cl=[B1ss+B2ss*Dcss*D21ss
   %             Bcss*D21ss];
    Css_cl=[C1ss+D12ss*Dcss*C2ss     D12ss*Ccss];
   % Dss_cl=D11ss+D12ss*Dcss*D21ss;
    nx_cl=size(Ass_cl,1);

    %% verify internal stability of closed loop system
    [V, Lambda] = eig(Ass_cl);
    Lambda = diag(Lambda);
    Vn = Css_cl*V;
    wn = -1i*Lambda;
    [wn,ind] = sort(wn,'ComparisonMethod','real');
    Vn = Vn(:,ind); 
    %% outputs
    output.wo=wn;
    output.Vo=Vn;
    output.nx_cl=nx_cl;
end