function [s,iter]=multiNew(f,a,b,alfa,beta,df1,df2,s0,m) 
% DESCRIZIONE: funzione che risolve il metodo di shooting
% multiplo applicando il metodo di Newton.
% INPUT: 
% f=funzione del problema 
% a=intervallo temporale sx
% b=intervallo temporale dx
% alfa=valore della soluzione nell'estremo sinistro
% beta=valore della soluzione nell'estremo destro
% df1=derivata di f rispetto alla variabile y(1) 
% df2=derivata di f rispetto alla variabile y(2) 
% m=numero di nodi 
%
% OUTPUT: 
% s=vettore [alfa,s1,alfa2,s2,...,beta,sm] contenente il valore della
% soluzione e della sua derivata ad ogni nodo
% iter=numero di iterazioni necessarie 
options=odeset('AbsTol',1e-4);

maxiter=100;
% Uso un'approssimazione iniziale della soluzione data dalla retta passante
% per i punti (a,alfa) e (b,beta) 
 y1=@(x)alfa+(beta-alfa)*(x-a)/(b-a);

 g=@(t,y) [f(t,y);y(4); df1(y)*y(3)+df2(y)*y(4)]; % sistema ausiliario


x1=linspace(a,b,m);
% trovo un'approssimazione dei valori della soluzione nei nodi x_k
alfa1=y1(x1);  
% trovo un'approssimazione dei vaolori della derivata prima della soluzione
% nei nodi x_k

s1=zeros(1,2*m);
s1(1:2:2*m)=alfa1;
s1(2:2:2*m)=s0;

tol=1e-4;
iter=0;
sx=s1+1;
sx=sx';
s1=s1';
% METODO DI NEWTON 
while norm(sx-s1)>tol && iter<maxiter
    iter=iter+1;
     sx=s1; 
     F=zeros(2*m,1);
     JF=zeros(2*m,2*m); 
     j=2;
     K=10;
     
     for k=1:m-1
        y0=[sx(j-1);sx(j);0;1];
        x2=linspace(x1(k),x1(k+1),K);
        [x,y]=ode45(g,x2,y0,options); % risoluzione del sistema ausiliario
        F1=y(end,1)-sx(j+1); % F(s)
        F2=y(end,2)-sx(j+2);
        F(j+1:j+2)=[F1;F2];
        p=f(0,y(end,1:2));
        G=[y(end,2),y(end,3); p(2),y(end,4)];
        JF(j+1:j+2,j-1:j+2)=[G,-eye(2)];
        j=j+2;
     end 
     
     y0=[alfa;sx(2);0;1];
     [x,Y]=ode45(g,[a b],y0,options);
     F(1)=Y(1,1)-alfa;
     F(2)=Y(end,1)-beta;
     JF(1,1)=1;
     JF(2,2)=Y(end,3);
     s1=sx-inv(JF)*F;
     
end
s=s1;
