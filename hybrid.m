function x = hybrid(f, dfdx, xmin, xmax, tol1, tol2)

% f: Function whose root is sought.
% dfdx: Derivative function.
% xmin: Initial bracket minimum.
% xmax: Initial bracket maximum.
% tol1: Relative convergence criterion for bisection.
% tol2: Relative convergence criterion for Newton iteration.
% x: Estimate of root.


f = @f;
dfdx= @dfdx;
xmin = 0;
xmax = 0.1;
tol1 = 0.000000000001;

converged = false;
	
% Bisection 
while converged == false
    xmid = (xmin + xmax) / 2;
    fmid = f(xmid);
    fmax = f(xmax);
    if fmid == 0 
        break
    elseif (fmid * fmax) < 0
        xmin = xmid;
    else 
        xmax = xmid;
    end
    if (xmax - xmin) / abs(xmid) < tol1
        converged = true;
    end
end

xn = xmid;

% Newton's method 
tol2 = 0.000000000001;
count = 0;
change = 0;
converge2 = false;

while converge2 == false
    if count == 100
        break
    elseif count < 100
        count = count + 1;
        xnew = xn - f(xn)/dfdx(xn);
        if (abs(xnew-xn) < tol2)
            converge2 = true
        end
    end  
end

x = xnew;
end





