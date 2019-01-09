function oeh(f,t0,y0,b,h)
%This matlab function solves the initial value problem 
%y'-f(t,y), y(t0)=y0, on the interval [t0,b] usig Euler's method with n
%subinervals. 

fprintf('\n')
disp('           The Euler Method              ')
disp('__________________________________________')
disp('    i     ti       yi       ')
fprintf('\n')

n=(b-t0)/h; 

t(1)=t0;
y(1)=y0;

%p.8 documentation 
fprintf('%6.0f, %6.2f, %12.2f \n', 0, t(1), y(1))


%shift unit from zero to one (matlab cant use 0)
for i=1:n
    t(i+1)=t(i)+h; 
    y(i+1)=y(i)+h*f(t(i),y(i)); 
    fprintf('%6.0f, %6.2f, %12.2f \n',i,t(i+1),y(i+1))
end

plot(t,y,'r-')
xlabel('t-axis')
ylabel('y-axis')

grid on 
hold on 

syms t y(t)
disp('The exact solution is:')
y=dsolve(diff(y)==f(t,y(t)), y(t0)==y0); 
y=simplify(y)
fplot(y,[t0,b])
hold off
legend('Numerical Solution', 'Exact Solution')
xlabel('t-axis')
ylabel('y-axis')

end

