function bisection(f,a,b,tol,n)
% Bisection method for solving f(x)=0 in [a,b].
% tol=error tolerance
% n= maximum number of iterrations 
a0=a;
b0=b;
iter=0;
u=f(a);
v=f(b);
c=(a+b)/2;
err=abs(b-a)/2;
disp('__________________________________________________')
disp('iter     a        c        b        f(c)    err')
disp('__________________________________________________')
fprintf('\n')
if (u*v<=0)
    while (err>tol)&&(iter<=n)
    w=f(c);
    if w==0
        fprintf('Exact solution is c=%.f\n',c)
        break
    else
        fprintf('%2.0f %10.4f %10.4f %10.4f %10.4f %10.4f\n',iter,a,c,b,f(c),err)
        if (w*u<0)
            b=c;v=w;
        else
            a=c;u=w;
        end
        iter=iter+1;
        c=(a+b)/2;
        err=err/2;
    end
    end
    if (iter>n)
        disp('Method failed to converge')
    end
else
    disp('The method may not be applicable.')
end
%Plot f(x) in [a,b].
ezplot(f,[a0 b0])
title('Graph of the function')
xlabel('x');ylabel('y')
grid on 
end

