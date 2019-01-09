function oe3x3(f,t0,y0v,y1v,y2v,b,N)
% This MATLAB function solves the initial value problem
% x'=f(t,x,y), y'=g(t,x,y), x(t0)=x0, y(t0)=y0, 
% on the interval [t0,b] using Euler's method with n subintervals.
fprintf('\n')
disp('                The (vectorized) Euler Method             ')
disp('------------------------------------------------------')
disp('    n       tn        y0n        y1n        y2n    ')
fprintf('\n')
t(1)=t0;
y0(1)=y0v;
y1(1)=y1v;
y2(1)=y2v;

x1



fprintf('%6.0f %6.2f %12.2f %12.2f %12.2f \n',0,t(1),y0(1),y1(1),y2(1))
h=(b-t0)/N;
for n=1:N
    t(n+1)=t(n)+h;
    y0(n+1)=y0(n)+h*g(t(n),y0(n),y1(n),y2(n));
    y1(n+1)=y1(n)+h*g(t(n),y0(n),y1(n),y2(n));
    y2(n+1)=y2(n)+h*g(t(n),y0(n),y1(n),y2(n));
    fprintf('%6.0f %6.2f %12.2f %12.2f %12.2f \n',n,t(n+1),y0(n+1),y1(n+1),y2(n+1))
end
plot(t,y0,'b*',t,y1,'r*',t,y2,'g*')
hold on
syms t y0(t) y1(t) y2(t)
disp('The exact solution is:')
[y0,y1,y2]=dsolve(diff(y0)==f(t,y0,y1,y2),diff(y1)==f(t,y0,y1,y2),diff(y1)==f(t,y0,y1,y2),y0(t0)==y0v,y1(t0)==y1v,y2(t0)==y2v);
y0=simplify(y0)
y1=simplify(y1)
y2=simplify(y2)
fplot(y0,[t0 b])
hold on
fplot(y1,[t0 b])
hold on
fplot(y2,[t0 b])
hold off
legend('Numerical solution y0','Numerical solution y1','Numerical solution y2','Exact solution y0','Exact solution y1','Exact solution y2')
grid on
xlabel('t-axis')
ylabel('y-axis')
end



Example problem: The angle y of an undamped pendulum with a driving force sin(5 t) satisfies the differential equation

y'' = -sin(y) + sin(5 t)

and the initial conditions

y(0) = 1
y'(0) = 0.

If your problem is of order 2 or higher: rewrite your problem as a first order system. Let y1=y and y2=y', this gives the first order system

y1' = y2,
y2' = -sin(y1) + sin(5 t)

with the initial conditions

y1(0) = 1
y2(0) = 0.
