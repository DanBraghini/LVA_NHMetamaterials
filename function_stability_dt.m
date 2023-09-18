function output = function_stability_dt(Ad,varargin)
%function output = function_stability_dt(A,varargin)
%
% Verify the stability of a discrite-time linear system using linear
% matrix inequalities. The LMIs are  programmed using YALMIP and can be solved 
% by any LMI solver supported by YALMIP (SeDuMi was used with the default options).
% inputs:  A   -> dynamic matrix (cell array)
%          solver (optional)-> solver used to solve the LMIs. The default
%                              is SeDuMi.
%
% outputs: output.feas     -> stable (1) unstable (0)
%          output.cpusec   -> cpu time to solve the LMIs (seconds)
%          output.P        -> Solution variables P
%          output.V        -> number of scalar variables used in the optmization problem
%          output.L        -> number of LMI rows used in the optmization problem
%          output.mu       -> Solution variable mu         
%

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
order = size(Ad,1);

%% Decision variables

% Variable of the problem
P = sdpvar(order,order,'symmetric');

%% Constraints (LMIs)

% Stability condition
LMIs = [P>=0, Ad'*P*Ad - P + eye(order) <= 0];

%% Numerical complexity of the problem

% determine the number of LMI rows
output.L = 0;
for i=1:size(LMIs,1)
    output.L = output.L + size(sdpvar(LMIs(i)),1);
end
% determine the number of scalar variables
output.V = size(getvariables(LMIs),2);

%% Solving the LMIs

sol = optimize(LMIs,trace(P),sdpsettings('verbose',0,'solver',options.solver));
% evaluate the elapsed time to solve the LMIs set
output.cpusec = sol.solvertime;

%% Checking the feasibility and saving the solution (if exists)

output.feas = 0;
if sol.problem == 1
    disp('sol.problem == 1, the problem is asured to be unfeasable')
    return;% the solver asserted unfeasibility
end
% retrieving the minimal primal residual
pr=min(checkset(LMIs));
% capturing the solutions (if ones exist)
if pr > -10^-7
    output.P = double(P);
    output.feas = 1;
end
