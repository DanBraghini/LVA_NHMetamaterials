function u_FEM= function_FRF_FEM(F,wv,M,C,K,e,n1,n2,gamma_c_ph,ndofo,ne,boundary)
%  u_FEM= function_FRF_FEM(F,wv,M,C,K,e,n1,n2,gamma_c_ph,ndof,ne,boundary);
% Danilo Braghini 08, January, 2022
% Given an external force F at the excitation pont e, and structure's
% FEM discretization (Ms,Cs,Ks), this algorithm computs the response 
% on the dofs of the FEM model, for each frequency of the frequency vector wv.
% Control is introduced on the physical model (M,C,K) as
% gamma_c_ph, on the unit cell. The FRF is computed for an
% finite system (metasctructure), thus we need to assemble 
% the matrices for ncells cells, each one with ne finite elements, which 
% results on ndof degrees of freedom.
% output: complex response of the whole structure  in pressure by volume velocity (acoustic)
% on the range of frequencies given by wv, to impulsive excitation. 
    %% for PID local control  
    % Dispersion Relation uses frequency domain equations, so we can use 
    % Bloch Theorem here like we did before on SEM
    M(n2,n1) = M(n2,n1)-gamma_c_ph(3);
    K(n2,n1) = K(n2,n1)-gamma_c_ph(2);
    C(n2,n1) = C(n2,n1)-gamma_c_ph(1);
    %% joining cells to make a structure
    Ms=zeros(ndofo,ndofo);Cs=Ms;Ks=Ms;
    for i=1:ne:ndofo-ne
        Ms(i:i+ne,i:i+ne) = Ms(i:i+ne,i:i+ne) + M;
        Cs(i:i+ne,i:i+ne) = Cs(i:i+ne,i:i+ne) + C;
        Ks(i:i+ne,i:i+ne) = Ks(i:i+ne,i:i+ne) + K;
    end
    %% correcting dimensions depending on boundary conditions
    if boundary==1
       ndof=ndofo-2;
    elseif boundary==2
      ndof=ndofo-1;
    elseif boundary==0
      ndof=ndofo;
    end
    nf=length(wv);
    u_FEM = zeros(ndof,nf);
    %% frequency loop
    for i = 1:nf
        w = wv(i);
        % Dinamic stiffness matrix via FEM
        Dg = -w^2*Ms + 1i*w*Cs + Ks;
        % External force in constant volume velocity 
        % PS: In this discretization, source is -B/A times volume velocity.
        % Also, the matrixes (Ms,Cs,Ks) were previously multiplied by A*rho
        %% External force in volume velocity 
        Fe=F*1i*w;
        %% force in the middle
        if e == 'm'
            m = floor(ndof/2);
            if mod(ndof,2) == 0
                 Fv = [zeros(ndof-m-1,1);Fe;zeros(ndof-m,1)];
            else
                Fv = [zeros(m,1);Fe;zeros(m,1)];
            end
        % force on the left end
        elseif e == 'l'
            Fv = [Fe; zeros(ndof-1,1)];
        % force on the right end    
        elseif e=='r'
            Fv = [zeros(ndof-1,1);Fe];
        else 
            disp('invalid force position')
            return
        end
        %% applying boundary conditions
    if boundary==1
        Dg(1,:)=[];Dg(:,1)=[];Dg(end,:)=[];Dg(:,end)=[];
    elseif boundary==2
        Dg(:,1)=Dg(:,1)+Dg(:,end); Dg(1,:)=Dg(1,:)+Dg(end,:);Dg(1,1)=Dg(1,1)+Dg(end,end);
        Dg(:,end)=[];Dg(end,:)=[];
    end
        %% solving linear system
        u_FEM(:,i) = Dg\Fv;
    end
 end