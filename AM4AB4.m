function AM4AB4(f,t0,y0,b,n)
% This MATLAB function solves the initial value problem 
% y'=f(t,y), y(t0)=y0, on the interval [t0,b],
% using the 2-step Adams-Moulton method with n grid points (or n steps). 
% If possible, it also provides the exact solution. 
% For example, type eulern(@(t,y) t*y+1,0,1,2,10) to solve
% y'=ty+1, y(0)=1 on [0,2] (with 10 knots, i.e, step size h=0.2).
fprintf('\n')
disp('               The A-M 2 step Method')
disp('_____________________________________________')
disp('    i     ti       yi   ')
disp('_____________________________________________')
fprintf('\n')
y(1)=y0;
t(1)=t0;
fprintf('%6.0f %6.2f %12.6f \n',0,t(1),y(1))
h=(b-t0)/n;
t(2)=t(1)+h;

%RK TO APROXIMATE FIRST 4
for i=1:3
    t(i+1)=t(i)+h;
    k1=f(t(i),y(i));
    k2=f(t(i)+h/2,y(i)+h*k1/2);
    k3=f(t(i)+h/2,y(i)+h*k2/2);
    k4=f(t(i)+h,y(i)+h*k3);
    y(i+1)=y(i)+h*(k1+2*k2+2*k3+k4)/6;
    fprintf('%6.0f %6.2f %12.6f \n',i,t(i+1),y(i+1))
end   

for i=4:n
   t(i+1)=t(i)+h;
   
   f1=f(t(i),y(i));
   f2=f(t(i-1),y(i-1));
   f3=f(t(i-2),y(i-2)); 
   
   w=y(i)+(h/24)*(55*f(t(i),y(i))-59*f(t(i-1),y(i-1))+37*f(t(i-2),y(i-2))-9*f(t(i-3),y(i-3))); % predictor of y(i+1) by using 4-step A-B method
   f5=f(t(i+1),w); % it approximates f(t(i+1),y(i+1))
   y(i+1)=y(i)+ h*(9*f5+19*f1-5*f2+f3)/24; % corrector stage
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


 
