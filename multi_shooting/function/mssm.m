function [s,iter]=mssm(f,a,b,alfa,beta,df1,df2,m) 
% METODO DI SHOOTING MODIFICATO 
% DESCRIZIONE:  ad ogni iterazione k risolvo con il metodo di shooting semplice 
% sull'intervallo [a,x(k)] per k=2...m-1 considerando beta del metodo di
% shooting semplice pari al valore dell'approssimazione iniziale y1 in
% x(k). In questo modo, trovo un s0 iniziale con un approssimazione migliore
% per applicare il metodo di shooting semplice. 
% 
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
% s=valore di y'(a)
% iter=numero di iterazioni necessarie 
options=odeset('AbsTol',1e-4);

maxiter=100; % numero max di iterazioni metodo di Newton

 y1=@(x)alfa+(beta-alfa)*(x-a)/(b-a); % approssimazione iniziale 

 g=@(t,y) [f(t,y);y(4); df1(y)*y(3)+df2(y)*y(4)]; % funzione sistema ausiliario


x1=linspace(a,b,m);

s0=(beta-alfa)/(b-a); % approssimazione iniziale per s0


tol=1e-4; % tolleranza assoluta di ode45 

iter=0;

s1=s0+1;
F=1;

% per ogni k nodo applico il metodo di Newton per trovare lo zero di F
for k=2:m-1
    while abs(s0-s1)>tol && iter<maxiter
            iter=iter+1;
            s0=s1;

            y0=[alfa;s0;0;1];
            [x,y]=ode45(g,[a x1(k)],y0,options);
            F=y(end,1)-y1(x1(k));
            s1=s0-(F/(y(end,3)));
    end
       
        F=1;
end
% dopo aver trovato un'approssimazione iniziale di s0 applico il metodo di
% shooting semplice sull'intervallo [a,b]
% end

[s,iter] = shootingNewton(f,df1,df2,s0,a,b,alfa,beta);
 
    
    