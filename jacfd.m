% Compute jacobian matrix for system of d = 3 equations
% approx jacobian element using difference quotient 
% x is a length d vector 
% h is finite difference parameter 
function [J] = jac(g,x,h)
g = @g;
gx=feval(g,x);
xp=x;

for i=1:3
    xp(i)=xp(i)+h;
    J(:,i)=(feval(g,xp)-gx)/h;
    xp(i)=x(i);
end
end



