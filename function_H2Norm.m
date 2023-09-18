function output = function_H2Norm(A,B,C,varargin)
%function output = hinf_norm_c(A,B,C,D,varargin)
%
% Compute the H-infininty norm of a continuous-time linear system using linear
% matrix inequalities. The LMIs are  programmed using YALMIP and can be solved 
% by any LMI solver supported by YALMIP (SeDuMi was used with the default options).
% inputs:  A,B,C,D          -> matrices of the system.
%          solver (optional)-> solver used to solve the LMIs. The default
%                              is SeDuMi.
%          h2 (optional)  -> if this parameter is different of zero, the routine will
%                              check the feasibility of the LMIs for gamma = h2.
%          tol (optional)   -> tolerance of the feasbility of the LMIs. 
%
% outputs: output.h2     -> value of the H-2 norm (0 if unfeasible)
%          output.cpusec   -> cpu time to solve the LMIs (seconds)
%          output.Pk       -> Solution variables P
%          output.K        -> number of scalar variables used in the optmization problem
%          output.L        -> number of LMI rows used in the optmization problem
%          output.delta    -> minimal primal residual returned by the LMI solver (SeDuMi is the default).
%

if nargin > 4
    if nargin == 5
        options = varargin{1};
    else
        options = struct(varargin{:})
    end
else
    options = [];
end

if ~isfield(options,'solver')
    options.solver = 'sedumi';
end
if ~isfield(options,'h2')
    options.h2 = 0;
end
if ~isfield(options,'tol')
    options.tol = 1e-7;
end

%% Dimensions of the system

%determine the number of vertices of polytope
order = size(A,1);

%% Decision variables

% Lyapunov matrix
P = sdpvar(order,order,'symmetric');
X = sdpvar(order,order,'full');
% F =sdpvar(order,order,'full');
%G =sdpvar(order,order,'full');

% if the user wants to calculate the norm, define objective function
% otherwise, the norm is input to the function.
if options.h2 == 0
    gamma = sdpvar();
    obj = gamma;
else
    gamma = options.h2;
    obj = [];
    options.tol = 0;
end

%% Constraints (LMIs)

% Controlability gramian
%LMIs = [P>=0, gamma>=trace(X), X>=C*P*C', A*P+P*A'+B*B'<=0];
% Observability gramian
LMIs = [P>=0, gamma>=trace(B'*P*B), A'*P+P*A+C'*C<=0];

%% Finsler lemma
% inputs = size(B,2);
% M11=F*A';
% M12=P+A*G-F;
% M13=B;
% M22=-G-G';
% M33=-eye(inputs);
% M23=zeros(order,inputs);
% 
% M = [M11  M12  M13
%      M12' M22  M23
%      M13' M23' M33];
%  LMIs = [P>=0, gamma>=trace(C*P*C'), M<=0]; 
%% Projection lemma
% M11=C'*C+P-F-F';
% M12=A'*P+F';
% M22=-P;
% M21=M12';
% M = [M11 M12
%      M21 M22];
% LMIs = [P>=0, gamma>=trace(B'*P*B), M<=0]; 
 
%% Numerical complexity of the problem

% determine the number of LMI rows
output.L = 0;
for i=1:size(LMIs,1)
    output.L = output.L + size(sdpvar(LMIs(i)),1);
end
% determine the number of scalar variables
output.V = size(getvariables(LMIs),2);

%% Solving the LMIs

% solve the LMIs
sol = optimize(LMIs,obj,sdpsettings('verbose',0,'solver',options.solver));
% evaluate the elapsed time to solve the LMIs set
output.cpusec = sol.solvertime;

%% Checking the feasibility and saving the solution (if exists)
% retrieving the minimal primal residual
output.pr=min(checkset(LMIs));

output.feas = 0;
if sol.problem == 1
    return;
end

% capturing the solutions (if ones exist)
if output.pr > -options.tol
    output.W = double(P);
    output.h2 = sqrt(double(gamma));
    output.feas = 1;
end