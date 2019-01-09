function [xmin,ymin] = ournewton(f,x0,tol)
% This function finds the minimum by Newton's method.
f=sym(f);
syms x
df(x)=diff(f);
ddf(x)=diff(f,2);
f=matlabFunction(f);
df=matlabFunction(df);
ddf=matlabFunction(ddf)
x1=x0-df(x0)/ddf(x0);
while abs(x1-x0)>tol
    x0=x1;
    x1=x0-df(x0)/ddf(x0);
end
xmin=x1;
ymin=f(xmin);
end

