% BRACHISTOCHRONE
% This code solves the brachistocrone problem with a multiple shooting
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

% analitical solution
xesatta=@(t) r*(t-sin(t))/2+y(2); 
yesatta=@(t) r*(1-cos(t))/2;

dy=@(t) r*sin(t)/2; 



A=diag(2/3,-1);
c=[0,2/3];
b1=[1/4,3/4];

counter=0;
m1=10;
s0=ones(1,m1);
tic
s=solveMulti(f1,xa,xb,alfa,beta,m1);
toc
x=linspace(a,b,m1);
sx=dy(a);
norm(s(2)-sx,inf)

s01=1; 
df1=@(y) (1+y(2)^2)/(2*y(1)^2);
df2=@(y) y(2)/y(1);
tic
 [sN,iter1] = shootingNewton(f1,df1,df2,s01,xa,xb,alfa,beta);
 toc
norm(sN-sx,inf) 
mrange=10:10:100;
iter=0; 
 
 for m=mrange+1
     counter=counter+1;
     iter=0;
     m1=25;
     s0=ones(1,m1);
    s=solveMulti(f1,xa,xb,alfa,beta,m1);

     x0=linspace(a,b,m1*(m1-1));

     x2=linspace(xa,xb,m1);
    x3=linspace(a,b,m1);   

    j=1;
    g=1;
    g1=m1;
 for k=1:m1-1
     iter=iter+1;
    y0=[s(j);s(j+1)];
    x1=linspace(x2(k),x2(k+1),m1);
    x11=linspace(x3(k),x3(k+1),m1);
    h=(x2(k+1)-x2(k))/(m1-1);
    y1=zeros(2,m1);
    y1(:,1)=y0;
    for n=1:m1-1
        y1(:,n+1)=y1(:,n)+h*f1(1,y1(:,n));
    end

    Y1(g:g1)=y1(1,:);
    X(g:g1)=x1;
    g=g1+1; 
    g1=g+m1-1;
    j=j+2;
     err1(iter)=norm(yesatta(x11)-y1(1,:),inf);
 end
    err(counter)=sum(err1)/m;
 end
  

figure 
plot(X,Y1,'ro','linewidth',2)
  hold on
    plot(xesatta(x0),yesatta(x0),'b','linewidth',3)
axis equal
legend('numerical solution','analitical solution')

% grafico errore 
ord1 = err(1)*(mrange/mrange(1)).^(-1);
ord2 = err(1)*(mrange/mrange(1)).^(-2);
ord3 = err(1)*(mrange/mrange(1)).^(-3);
ord4 = err(1)*(mrange/mrange(1)).^(-4);

figure
loglog(mrange,err,'*r',mrange,ord1,mrange,ord2,mrange,ord3,mrange,ord4)
legend('order','h','h^2','h^3','h^4')
xlabel('m')
ylabel('error')
