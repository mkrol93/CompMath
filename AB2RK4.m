function AB2RK4(f,t0,y0,b,n)
% This MATLAB function solves the initial value problem 
% y'=f(t,y), y(t0)=y0, on the interval [t0,b],
% using the 2-step Adams-Bashforth method with n grid points (or n steps). 
% If possible, it also provides the exact solution. 
% For example, type eulern(@(t,y) t*y+1,0,1,2,10) to solve
% y'=ty+1, y(0)=1 on [0,2] (with 10 knots, i.e, step size h=0.2).
fprintf('\n')
disp('               The A-B 2 step Method')
disp('_____________________________________________')
disp('    i     ti       yi   ')
disp('_____________________________________________')
fprintf('\n')
y(1)=y0;
t(1)=t0;
fprintf('%6.0f %6.2f %12.6f \n',0,t(1),y(1))
h=(b-t0)/n;

t(2)=t(1)+h;
    k1=f(t(1),y(1));
    k2=f(t(1)+h/2,y(1)+h*k1/2);
    k3=f(t(1)+h/2,y(1)+h*k2/2);
    k4=f(t(1)+h,y(1)+h*k3);
y(2)=y(1)+h*(k1+2*k2+2*k3+k4)/6; % Runge-Kutta for the second step

fprintf('%6.0f %6.2f %12.6f \n',0,t(2),y(2))
for i=2:n
   t(i+1)=t0+i*h;
   y(i+1)=y(i)+(h/2)*(3*f(t(i),y(i))-f(t(i-1),y(i-1)));
   fprintf('%6.0f %6.2f %12.6f\n',i,t(i+1),y(i+1));
end
plot(t,y,'r*')
xlabel('x-axis')
ylabel('y-axis')
legend('Approximate Solution')
hold on
syms t y(t)
disp('The exact solution is:')
y=dsolve(diff(y)==f(t,y),y(t0)==y0);
y=simplify(y)
fplot(y,[t0,b],'b-')
hold off
title('Graphs of the numerical and exact solutions')
grid on
legend(' Numerical solution y',' Exact solution y')
xlabel('x-axis')
ylabel('y-axis')


 
