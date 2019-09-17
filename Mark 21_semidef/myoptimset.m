function myoptions = myoptimset
% MYOPTIMSET Sets default options for user-developed nonlinear programming
% routine, both unconstrained and constrained.
%
%   OUTPUTS:
%           myoptions       =       structure containing all options
%                                   required by the NLP solvers
% 
% %% General options
% myoptions.display    	=	'Iter';     % Display iteration output
% myoptions.xsequence 	=	'off';      % Store sequence of points {xk}
% myoptions.tolgrad    	=	1e-6;       % Termination tolerance on the norm
%                                       % of the directional derivative
% myoptions.tolX       	=	1e-12;      % Termination tolerance on the relative
%                                       % change of optimization variables
% myoptions.tolfun    	=	1e-12;      % Termination tolerance on the relative
%                                       % improvement of the cost function
% myoptions.tolconstr   =	1e-6;       % Constraint tolerance satisfaction
%                                         tolerance
% myoptions.nitermax   	=	50;         % Termination tolerance on the number of
%                                       % iterations
%                                         
% %% Differentiation options
% myoptions.gradmethod 	=	'FD';       % Method for gradient computation
%                                       % FD    =   Forward Differences
%                                       % CD    =   Central Differences
%                                       % IM    =   Imaginary-part trick
%                                       % UP    =   User Provided
% myoptions.graddx     	=	2^-26;      % Perturbation for gradient computation
%                                       % use 2^-26 for FD
%                                       % use 2^-17 for CD
%                                         
% %% Line search options
% myoptions.ls_tkmax   	=	1;          % Maximum step size
% myoptions.ls_beta    	=	0.8;        % Beta scaling factor for 
%                                       % back-tracking line search
% myoptions.ls_c       	=	0.1;        % c coefficient factor for 
%                                       % back-tracking line search
% myoptions.ls_nitermax	=	20;         % max. number of back-tracking 
%                                       % line search iterations
% 
% %% Quasi-Newton method options
% myoptions.Hessmethod 	=	'BFGS';     % Method for Hessian computation
%                                       % Exact =   Exact Newton. In this
%                                                   case, function f(x) shall 
%                                                   provide as second
%                                                   output the gradient of
%                                                   f, and as third output 
%                                                   the Hessian of f
%                                       % BFGS  =   BFGS algorithm
%                                       % GN    =   Gauss-Newton algorithm
%                                       % GD    =   Gradient descent
% myoptions.BFGS_gamma 	=	0.2;        % gamma factor for Powell's trick in
%                                       % BFGS algorithm

%% General options
myoptions.display    	=	'Iter';     % Display iteration output
myoptions.xsequence    	=	'off';      % Store sequence of points {xk}
myoptions.tolgrad    	=	1e-6;       % Termination tolerance on the norm
                                        % of the directional derivative
myoptions.tolconstr    	=	1e-6;       % Constraint satisfaction tolerance
myoptions.tolx          =	1e-12;      % Termination tolerance on the relative
                                        % change of optimization variables
myoptions.tolfun        =	1e-12;      % Termination tolerance on the relative
                                        % improvement of the cost function
myoptions.nitermax      =	50;         % Termination tolerance on the number of
                                        % iterations
myoptions.outputfcn     =	[];         % Handle for output function

%% Differentiation options
myoptions.gradmethod  	=	'FD';       % Method for gradient computation
                                        % FD    =   Forward Differences
                                        % CD    =   Central Differences
                                        % IM    =   Imaginary-part trick
                                        % UP    =   User Provided
myoptions.graddx        =	2^-26;      % Perturbation for gradient computation
                                        % use 2^-26 for FD
                                        % use 2^-17 for CD
                                        
%% Line search options
myoptions.ls_tkmax      =	1;          % Maximum step size
myoptions.ls_beta       =	0.8;        % Beta scaling factor for 
                                        % back-tracking line search
myoptions.ls_c          =	0.1;        % c coefficient factor for 
                                        % back-tracking line search
myoptions.ls_nitermax   =	20;         % max. number of back-tracking 
                                        % line search iterations

%% Quasi-Newton method options
myoptions.Hessmethod  	=	'BFGS';     % Method for Hessian computation
                                        % Exact =   Exact Newton. In this
                                        %           case, it is expected that the
                                        %           cost function provides as second
                                        %           output the gradient, and as third
                                        %           output the Hessian matrix.
                                        % BFGS  =   BFGS algorithm
                                        % GN    =   Gauss-Newton algorithm
                                        % SD    =   Steepest Descent
myoptions.BFGS_gamma  	=	1e-1;       % gamma factor for Powell's trick in
                                        % BFGS algorithm
myoptions.GN_funF       =	[];         % function providing the value of F
                                        % in Gauss-Newton method, assuming
                                        % f(x)=F(x)'*F(x)
myoptions.GN_sigma      =	0;          % coefficient to ensure Hessian is 
                                        % positive definite in Gauss-Newton
                                        % method (H=GradF*GradF'+eye(n)*sigma)
myoptions.QPoptions     =   ...
    optimset('Display','none','Algorithm','interior-point-convex');         % options for QP solver in SQP methods
end

