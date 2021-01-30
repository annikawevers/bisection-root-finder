function x = newtond(g, jacfd, h, x, tol)

% f:     Function which implements the nonlinear system of equations.
%        Function is of the form f(x) where x is a length-d vector, and
%        returns length-d column vector.
% jacfd: Function which is of the form jacfd(f, x, h) where f is the above
%        function, x is a length-d vector, and h is the finite difference
%        parameter. jacfd returns the d x d matrix of approximate Jacobian
%        matrix elements.
% h:     Finite differencing parameter.
% x0:    Initial estimate for iteration (length-d column vector).
% tol:   Convergence criterion: routine returns when relative magnitude
%        of update from iteration to iteration is <= tol.
% x:     Estimate of root (length-d column vector)

x0 = [-1.00 ; 0.75 ; 1.50];
h = 10^(-5); 
tol = 10^(-12);  
jacfd = @jacfd; 
g = @g; 

dx = 1; 
normdx = 1;
count = 0;
%res = [0 ; 0; 0];  
%J = [[0, 0, 0] ; [0, 0, 0]; [0, 0, 0]]; 

x = x0;
while normdx > tol 
    if count < 1000
        count = count + 1;
        res = g(x);
        J = caller(@jacfd, @g, x, h);
        dx = J\res;
        x = x - dx
        normdx = sqrt((dx(1)^2 + dx(2)^2 + dx(3)^2)/3);
    else 
        break
    end
end
end

function ty = caller(jacfd, g, x, h)
    jacfd = @jacfd;
    g = @g;
    ty = jacfd(g, x, h);
end




   
  



