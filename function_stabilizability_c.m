function output = function_stabilizability_c(A,B,varargin)
%function output = function_stabilizability_c(A,B,varargin)
%
% Verify the stability of a continuous-time linear system using linear
% matrix inequalities. The LMIs are  programmed using YALMIP and can be solved 
% by any LMI solver supported by YALMIP (SeDuMi was used with the default options).
% inputs:  A   -> dynamic matrix (cell array)
%          solver (optional)-> solver used to solve the LMIs. The default
%                              is SeDuMi.
%
% outputs: output.feas     -> stable (1) unstable (0)
%          output.cpusec   -> cpu time to solve the LMIs (seconds)
%          output.P        -> Solution variables P
%          output.K        -> static stabilizable state feedback controller
%          output.V        -> number of scalar variables used in the optmization problem
%          output.L        -> number of LMI rows used in the optmization problem
%
% Example: 2 states
% A = [-0.9  0.2;
%       -0.5  -1.9 ];
% output = stability_c(A)
%
% Date: 02/10/2020
% Author: ricfow@dt.fee.unicamp.br

if nargin > 2
    if nargin == 3
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
order = size(A,1);
inputs = size(B,2);
%% Decision variables

% Variable of the problem
P = sdpvar(order,order,'symmetric');

% eps = 1*10^-8; % very small value to enforce positive definite
% epsI = eps*eye(order); % matrix needed to enforce positive definite 
%% Constraints (LMIs)

% Stability condition
LMIs = [P>=0, A'*P + P*A + B*B' <= 0];%-epsI];

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
if sol.problem == 1
    return;% the solver asserted unfeasibility
end
% retrieving the minimal primal residual
pr=min(checkset(LMIs));
% capturing the solutions (if ones exist)
if pr > 0
    output.P = double(P);
    output.K = (1/2)*B'*inv(double(P));
    output.feas = 1;
end