function [w_real_FEM,w_imag_FEM]= function_dipersion_FEM(ni,n1,n2,neig,M,C,K,gamma_c_ph,a,dmodel,kLv)
% [w_real_FEM,w_imag_FEM]= function_dipersion_FEM(ni,n1,n2,neig,M,C,K,gamma_c,a,dmodel,kLv);
% Danilo Braghini 08, January, 2022    
% Dispersion Diagram Aproximation via FEM
% FEM matrices given as inputs discribre the passive system. On this
% algorithm, control is introduced on the physical model (M,C,K) as
% gamma_c_ph, on the unit cell. The dispersion relation is computed for an
% infinity system (metamaterial) an due to periodicity we only need the
% cell's matrices.
% output: dispersion relation (real and imaginary parts of frequencies solutions 
% corresponding to each wavenumber) in rad/s
%% running trhough the wavenumber vector, created to calculate the dispersion relation
    N=length(kLv);    
    w_real_FEM =zeros(neig,N); w_imag_FEM = w_real_FEM;
    for n = 1:N
        kL = kLv(n);
        Tq = [1              zeros(1,ni) 
              zeros(ni,1)    eye(ni)
              exp(-1i*kL)    zeros(1,ni) ];

        %% for PID generally non-local control    
        % Dispersion Relation uses frequency domain equations, so we can use 
        % Bloch Theorem here like we did before on SEM
        M(n2,n1) = M(n2,n1)-gamma_c_ph(3)*exp(1i*a*kL);
        K(n2,n1) = K(n2,n1)-gamma_c_ph(2)*exp(1i*a*kL);
        C(n2,n1) = C(n2,n1)-gamma_c_ph(1)*exp(1i*a*kL);

        M_bar = Tq'*M*Tq;
        C_bar = Tq'*C*Tq;
        K_bar = Tq'*K*Tq;

        w_FEM = polyeig(K_bar,1i.*C_bar,-M_bar);

        %% The Quadratic Eigenvalue problem has 2*np eigenvalues, but we can take 
        % only those with positive real(omega)
        f=find(real(w_FEM) >= 0);
        w_FEM_aux = w_FEM(f);
        w_FEM = sort(w_FEM_aux,'ComparisonMethod','real');
        w_real_FEM(:,n) = real(w_FEM(1:neig));
        w_imag_FEM(:,n) = imag(w_FEM(1:neig));

        %% If C is the null matriz, polyeig shall be replaced for eig (better conditioned algorithm)    
         if dmodel == 'n' && gamma_c_ph(1) == 0  
             w_FEM_square = eig(K_bar,M_bar);
             w_FEM_square = sort(w_FEM_square,'ComparisonMethod','real');
             w_FEM = sqrt(w_FEM_square(1:neig));
             w_real_FEM(:,n) = real(w_FEM);
             w_imag_FEM(:,n) = imag(w_FEM);
         end
    end
end