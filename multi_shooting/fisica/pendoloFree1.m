% PENDULUM
% This code finds the time and initial velocity 
clc
clear all 
close all

l=1;
g=9.8;

% linear approximation of the pendulum equation
f=@(t,y) [y(2); -g/l*y(1)];

s0=pi/2*sqrt(l/g); % approximation 1
s1=pi*sqrt(l/g); % approximation 2

a=2/3*pi; % initial time

alfa= pi; % solution at final time
beta=0; % velocity at final time
gamma=pi/100; % initial angol 


% BISECTION METHOD
tic
[sB,vB,iterB]= shootingBisFree(f,alfa,beta,gamma,s0,s1,a);
toc

% NEWTON METHOD
df1=@(y) -g/l;
df2=@(y) 1;
% dato iniziale 2*s0 vicino alla radice

% METODO DELLE SECANTI
tic
 [sN,vN,iterN]=shootingNewFree(f,df1,df2,alfa,beta,gamma,2*s0,a);
toc
tic 
 [sS,vS,iterS]= shootingSecFree(f,s0,s1,alfa,beta,gamma,a);
toc

% NUMERICAL SOLUTION

% tableau order 2 
A=diag(2/3,-1);
c=[0,2/3];
b1=[1/4,3/4];

mrange=20:10:50;
counter=0; 

for m=mrange+1
    counter=counter+1;
    t=linspace(a,sS,m); 
    
    % analitical solution
    omega=sqrt(g/l);

    f1=@(y) [y(1)*cos(omega*a+y(2))-gamma;
        -y(1)*omega*sin(omega*a+y(2))-vS];
    y0=[1;-1];
    y=fsolve(f1,y0);

    sol=@(t) y(1)*cos(omega*t+y(2));
    
    % numerical solution 
    y0=[gamma;vS];
    Y=rk(f,a,sS,A,b1,c,m,y0);
    err(counter)=norm(sol(t)-Y(1,:),inf);

end


plot(t,Y(1,:),'ro','linewidth',2)
hold on 
plot(t,sol(t),'b','linewidth',3)
legend('numerical solution','analitical solution')xlabel('x')
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

