function AB4RK4(f,t0,y0,b,n)
% This MATLAB function solves the initial value problem 
% y'=f(t,y), y(t0)=y0, on the interval [t0,b],
% using the 4-step Adams-Bashforth method with n grid points (or n steps). 
% If possible, it also provides the exact solution. 
% For example, type AB4RK4(@(t,y) t*y+1,0,1,2,10) to solve
% y'=ty+1, y(0)=1 on [0,2] (with 10 knots, i.e, step size h=0.2).
fprintf('\n')
disp('               The A-B 4 step Method')
disp('_____________________________________________')
disp('    i     ti       yi   ')
disp('_____________________________________________')
fprintf('\n')

y(1)=y0;
t(1)=t0;

h=(b-t0)/n;

for z=1:4
    t(z+1)=t(z)+h;
    k1=f(t(z),y(z));
    k2=f(t(z)+h/2,y(z)+h*k1/2);
    k3=f(t(z)+h/2,y(z)+h*k2/2);
    k4=f(t(z)+h,y(z)+h*k3);
    y(z+1)=y(z)+h*(k1+2*k2+2*k3+k4)/6;
    fprintf('%6.0f %6.2f %12.6f \n',0,t(z),y(z))
end   


for i=4:n
   t(i+1)=t0+i*h;
   y(i+1)=y(i)+(h/24)*(55*f(t(i),y(i))-59*f(t(i-1),y(i-1))+37*f(t(i-2),y(i-2))-9*f(t(i-3),y(i-3)));
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


 
