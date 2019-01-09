function ot2(f,t0,y0,b,n)
% This MATLAB function finds both numerical and
% exact solutions to y'=f(t,y), y(t0)=y0 in [t0,b]
% using Taylor's method of order 2 with subintervals.

% First, compute first order partial derivatives ft and fy.
syms t y ft(t,y) fy(t,y)
f=sym(f);
ft(t,y)=diff(f,t);
fy(t,y)=diff(f,y);
% Then, transform f, ft, and fy from symbolic to MATLAB functions.
f=matlabFunction(f);
ft=matlabFunction(ft);
fy=matlabFunction(fy);
h=(b-t0)/n;
t(1)=t0;
y(1)=y0;
fprintf('\n')
disp('      Taylors method of order 2    ')
fprintf('\n')
disp('  i              ti           yi')
disp('-----------------------------------')

fprintf('%6.0f %12.4f %12.4f\n',0,t(1),y(1))
for i=1:n
    t(i+1)=t(i)+h;
    y(i+1)=y(i)+f(t(i),y(i))*h+(ft(t(i),y(i))+...
        fy(t(i),y(i))*f(t(i),y(i)))*(h^2/2);
    fprintf('%6.0f %12.4f %12.4f\n',i,t(i+1),y(i+1))
end
plot(t,y,'r*')
hold on
% Solve for the exact solution (if possible).
syms t y(t)
y=dsolve(diff(y)==f(t,y),y(t0)==y0);
y=simplify(y)
fplot(y,[t0 b])
hold off
title('Graphs of the numerical and exact solutions')
grid on
xlabel('t-axis')
ylabel('y-axis')
legend('Numerical solution','Exact solution')




end

