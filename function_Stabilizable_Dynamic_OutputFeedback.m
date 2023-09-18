function output = function_Stabilizable_Dynamic_OutputFeedback(A,B1,B2,C1,C2,D11,D12,D21,varargin)
%function output =  function_Stabilizable_Dynamic_OutputFeedback(A,B1,C1,C2,D21,varargin)
% inputs:  (A,B1,C1,C2,D11,D21   -> state space matrixes for discrite
%                             system with perturbations w, measured output
%                             y and estimated output z
%          solver (optional)-> solver used to solve the LMIs. The default
%                               is SeDuMi.
%
% outputs: output.feas     -> stable (1) unstable (0)
%          output.cpusec   -> cpu time to solve the LMIs (seconds)
%          output.controller   -> state space realization of the resulted
%          output.V        -> number of scalar variables used in the optmization problem
%          output.L        -> number of LMI rows used in the optmization problem
%
%%
if nargin > 5
    if nargin == 6
        options = varargin{1};
    else
        options = struct(varargin{:});
    end
else
    options = [];
end

if ~isfield(options,'solver')
    options.solver = 'sedumi';
end

%% Dimensions of the system

% number of states
nx = size(A,1);
% number of exogenous inputs
nw = size(B1,2);
% number of applied feedback inputs
nu = size(B2,2);
% number of estimated outputs
nz = size(C1,1);
% number of measured outputs
ny = size(C2,1);

%% Decision variables

% Variables of the problem
% Blocks of Lyapunov matrix
X1 = sdpvar(nx,nx,'symmetric');
Y1 = sdpvar(nx,nx,'symmetric');
% other variables 
An=sdpvar(nx,nx,'full');
Bn=sdpvar(nx,ny,'full');
Cn=sdpvar(nu,nx,'full');
Dn=sdpvar(nu,ny,'full');

% Define if we have an minimization or feasibility problem
if options.hinf == 0
    % Minimization of objective function gamma
    gamma = sdpvar();
    obj = gamma;
    options.tol =-10^-7;
else
    % Feasibility problem
    gamma = options.hinf^2;
    obj = [];
    options.tol = 0;
end

H= A*Y+Y*A'+B2*L+L'*B2';
G= A'*X+X*A+F*C2+C2'*F';
T11=Y;
T12=eye(nx);
T22=X;
T= [ T11 T12
     T12' T22];
%% Constraints (LMIs)
% Stability condition
LMIs = [-H >=0, -G>=0, T >=0];

%% Numerical complexity of the problem

% determine the number of LMI rows
output.L = 0;
for i=1:size(LMIs,1)
    output.L = output.L + size(sdpvar(LMIs(i)),1);
end
% determine the number of scalar variables
output.V = size(getvariables(LMIs),2);

%% Solving the LMIs
sol = optimize(LMIs,[],sdpsettings('verbose',0,'solver',options.solver));
% evaluate the elapsed time to solve the LMIs set
output.cpusec = sol.solvertime;

%% Checking the feasibility and saving the solution (if exists)

output.feas = 0;
% retriving calculated value of mu
%output.mu = double(mu);
if sol.problem == 1
    disp('sol.problem == 1, the problem is asured to be unfeasable')
    return;% the solver asserted unfeasibility
end
% retrieving the minimal primal residual
pr=min(checkset(LMIs));
% capturing the solutions (if ones exist)
if pr > 0
% recovering controller variables(except Df, which is a variable of the problem)
    Ud =double(U);
    Vd = double(V);
    Yd = double(Y);
    Xd = double(X);
    Ld =double(L);
    Fd = double(F);
    Z = A + Yd*A'*Xd+Ld'*B2'*Xd+Yd*C2'*Fd';
    M=-Z;
    Ac= (inv(Vd)*M*inv(Ud'))';
    Bc =inv(Ud)*Fd;
    Cc = Ld*inv(Vd');
   
    Dc = 0;
    controller = ss(Ac,Bc,Cc,Dc);
    output.H = controller;
    output.feas = 1;
end