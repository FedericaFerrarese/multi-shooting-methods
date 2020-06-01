% BRACHISTOCHRONE
% This code solves the brachistocrone problem with a simple shooting
% method.
clc
clear all
close all

% function
f1=@(t,y) [y(2);-(1+(y(2))^2)/(2*y(1))];

% boundary conditions
xa=0;
ya=1;

xb=1;
yb=1.1;


f=@(y) [y(1)*(y(3)-sin(y(3)))/2+y(2)-xa;...
                y(1)*(1-cos(y(3)))/2-ya;...
                y(1)*(y(4)-sin(y(4)))/2+y(2)-xb;...
                y(1)*(1-cos(y(4)))/2-yb];
            
              
y0=[1;-1;2.5;3.6];

y=fsolve(f,y0);

r=y(1);
% time interval
a=y(3);
b=y(4);  
% boundary values
alfa=ya; 
beta=yb; 

% initial approximation
s0=1;


% Newton method
 
df1=@(y) (1+y(2)^2)/(2*y(1)^2);
df2=@(y) y(2)/y(1);
[sN,iter1] = shootingNewton(f1,df1,df2,s0,xa,xb,alfa,beta)
% [sN,iter] = shootingbisezione(f1,s0,s1,xa,xb,alfa,beta)
%[sN,iter] = shootingSecanti(f1,s0,s1,xa,xb,alfa,beta)
% soluzione analitica in forma parametrica 
xesatta=@(t) r*(t-sin(t))/2+y(2); 
yesatta=@(t) r*(1-cos(t))/2;


% analitical vs numerical solution

mrange=10:10:100;
iter=0;

% ordine 2 
A=diag(2/3,-1);
c=[0,2/3];
b1=[1/4,3/4];
for m=mrange+1
    iter=iter+1;
    t1=linspace(xa,xb,m)';
    t2=linspace(a,b,m)';
   
    y0=[alfa;sN];
    
    k=(xb-xa)/(m-1);
    Y=zeros(2,m);
    Y(:,1)=y0;
    for n=1:m-1
        Y(:,n+1)=Y(:,n)+k*f1(1,Y(:,n));
    end
    err(iter)=norm(yesatta(t2)-Y(1,:)', inf); 
end

% plots
figure
plot(t1,Y(1,:),'ro','linewidth',2)
hold on
plot(xesatta(t2),yesatta(t2),'b','linewidth',3)

axis equal
legend('numerical solution','analitical solution')
xlabel('x')
ylabel('y(x)') 

    
ord1 = err(1)*(mrange/mrange(1)).^(-1);
ord2 = err(1)*(mrange/mrange(1)).^(-2);
ord3 = err(1)*(mrange/mrange(1)).^(-3);
ord4 = err(1)*(mrange/mrange(1)).^(-4);

figure
loglog(mrange,err,'*r',mrange,ord1,mrange,ord2,mrange,ord3,mrange,ord4)
legend('order','h','h^2','h^3','h^4')
xlabel('m')
ylabel('error')
