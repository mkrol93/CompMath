function falsep(f,a,b,tol,n)
% False position method for solving the nonlinear
% equation f(x)=0.
a0=a;
b0=b;
iter=0;
u=f(a);
v=f(b);
p=(v*a-u*b)/(v-u);
w=f(p);
disp('_______________________________________________________')
disp(' iter     a             b          c            f(p)   ')
disp('_______________________________________________________')
fprintf('\n')
if (u*v<=0)
   while (abs(w)>tol)&&(abs(b-a)>tol)&&(iter<=n)&&((v-u)~=0)
      w=f(p);
      if w==0
          fprintf('Exact solution is p=%.f\n',p)
          break
      else
      fprintf('%2.0f %12.4f %12.4f  %12.6f  %10.6f\n',iter,a,b,p,w)
      if (w*u<0)
         b=p;v=w;
      else
         a=p;u=w;
      end
      iter=iter+1;
      p=(v*a-u*b)/(v-u);
      end
   end
   if (iter>n)
      disp('  Method failed to converge')
   end
   if (v-u==0)
      disp('  Division by zero')
   end
else
   disp('  The method cannot be applied f(a)f(b)>0')  
end 
% Plot f(x) in the interval [a,b].
fplot(f,[a0 b0])
xlabel('x');ylabel('f(x)'); grid

  
 
      
         
   