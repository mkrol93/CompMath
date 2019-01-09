function rk4(f,t0,y0,b,N)
% This MATLAB function finds both numerical and
% exact solutions to y'=f(t,y), y(t0)=y0 in [t0,b]
% using Runge-Kutta's method of order 4 with N subintervals.
h=(b-t0)/N;
t(1)=t0;
y(1)=y0;
fprintf('\n')
disp('      RK method of order 4    ')
fprintf('\n')
disp('  n              tn           yn')
disp('-----------------------------------')

fprintf('%6.0f %12.4f %12.4f\n',0,t(1),y(1))
for n=1:N
    t(n+1)=t(n)+h;
    k1=f(t(n),y(n));
    k2=f(t(n)+h/2,y(n)+h*k1/2);
    k3=f(t(n)+h/2,y(n)+h*k2/2);
    k4=f(t(n)+h,y(n)+h*k3);
    y(n+1)=y(n)+h*(k1+2*k2+2*k3+k4)/6;
    fprintf('%6.0f %12.4f %12.4f\n',n,t(n+1),y(n+1))
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

