function Tn = tr(f,a,b,h)
% This MATLAB function appromates the integral 
% of f from a to be by using the trapezoidal rule
% with n subintervals.
n=(b-a)/h;

Tn=(f(a)+f(b))/2;
for i=1:n-1
    Tn=Tn+f(a+i*h);
end
Tn=h*Tn;
end

