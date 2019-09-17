function [fxk,gradient] = mygradient(fun,xk,method,dx)
% MYGRADIENT Computes the gradient (i.e. the Jacobian transpose) of a 
% given function, using one of several possible methods.
%   INPUTS:
%           fun         =   function whose gradient shall be evaluated (in
%                           general a vector of N functions)
%           xk          =   value of the n-dimensional input argument at which to evaluate
%                           the gradient (dfun(x)/dx)'
%           method      =   string indicating the differentiation method: 'FD' 
%                           (Forward Finite Differences), 'CD' (Central
%                           Finite Differences), 'UP' (User provided), 
%                           'IM' (Imaginary-part trick).
%                           If 'UP' is used, function fun(x) shall return as
%                           second output the Jacobian dfun(x)/dx evaluated
%                           at x = xk.
%           dx          =   perturbation step used in Finite Difference
%                           approximation
%
%   OUTPUTS:
%           fxk         =   Value of fun(xk)
%           gradient    =   n-by-N matrix containing in each column the
%                           partial derivatives of the corresponding element of function
%                           fun(x) with respect to each component of x,
%                           evaluated at xk
 
if strcmp(method,'UP')
    [fxk,gradient]  =   fun(xk);
elseif strcmp(method,'FD')
    n               =   length(xk);
    p               =   zeros(n,1);
    fxk             =   fun(xk);
    N               =   length(fxk);
    gradient        =   zeros(n,N);
    p(1,1)          =   dx;
    gradient(1,:)   =   (fun(xk+p)-fun(xk))'/dx;
    for ind = 2:n
        p(ind-1,1)      =   0;
        p(ind,1)        =   dx;
        gradient(ind,:) =   (fun(xk+p)-fun(xk))'/dx;
    end
elseif strcmp(method,'CD')
    n               =   length(xk);
    p               =   zeros(n,1);
    fxk             =   fun(xk);
    N               =   length(fxk);
    gradient        =   zeros(n,N);
    p(1,1)          =   dx;
    gradient(1,:)   =   (fun(xk+p)-fun(xk-p))'/(2*dx);
    for ind = 2:n
        p(ind-1,1)      =   0;
        p(ind,1)        =   dx;
        gradient(ind,:) =   (fun(xk+p)-fun(xk-p))'/(2*dx);
    end
elseif strcmp(method,'IM')
    n               =   length(xk);
    p               =   zeros(n,1);
    fxk             =   fun(xk);
    N               =   length(fxk);
    gradient        =   zeros(n,N);
    tol             =   1e-100;
    p(1,1)          =   1j*tol;
    gradient(1,:)   =   imag(fun(xk+p))'/tol;
    for ind = 2:n
        p(ind-1,1)      =   0;
        p(ind,1)        =   1j*tol;
        gradient(ind,:) =   imag(fun(xk+p))'/tol;
    end
end