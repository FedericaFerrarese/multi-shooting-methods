function [s,iter] = shootingNewton(f,df1,df2,s0,a,b,alfa,beta)
% DESCRIZIONE: Dato il problema di Cauchy
% y''(t)=f(t,y(t),y'(t))     con t in (a,b)
% y(a)=alfa
% y(b)=beta
%
% dove f,a,b,alfa, beta sono assegnati.
% Il metodo trasforma il problema ai limiti in uno di Cauchy del tipo:
%
% y1'(t)=y2(t)         con t in (a,b)
% y2'(t)=f(t,y1(t),y2(t))
% y1(a)=alfa
% y2(a)=s
%
% e determina il valore del parametro s che soddisfa la seguente relazione:
%             F(s)=|y2(b;s)-beta|<=tol
% mediante il metodo di Newton. Per il calcolo di F(s) e F'(s) e'
% necessario introdurre un sistema ausiliario rappresentato dalla funzione
% g. Quindi F(s)=y(1) e F'(s)=y(3) al tempo t=b.
%
% INPUT:
% f = funzione del problema di Cauchy
% df1 = derivata di f rispetto a a y1
% df2 = derivata di f rispetto a y2
% s0 = approssimazione iniziale di s
% a = estremo sinistro intervallo di integrazione
% b = estremo destro intervallo di integrazione 
% alfa = valore della soluzione nell'estremo sinistro
% beta = valore della soluzione nell'estremo destro
%
% OUTPUT:
% s = valore di y'(a) 
% iter = numero di iterazioni necessarie
options=odeset('AbsTol',1e-4);

g=@(t,y) [f(t,y);y(4); df1(y)*y(3)+df2(y)*y(4)]; % sistema ausiliario

y0=[alfa;s0;0;1]; % condizione iniziale sistema ausiliario

[x,y]=ode45(g,[a b],y0,options); %  risoluzione del sistema ausiliario
F=y(end,1)-beta; % F(s)

s1=s0-(F/(y(end,3))); % aggiornamento di s0 con il metodo di Newton

tol=1e-4; % tolleranza di ode45
maxiter=100; % numero massimo iterazioni
iter=0; % contatore iterazioni

% procedo risolvendo il sistema ausiliario aggiornando la condizione
% iniziale con il nuovo s0 trovato con il metodo di Newton fino a che il
% metodo converge alla soluzione cercata 
while abs(F) > tol && iter<=maxiter
    s0=s1;
    iter=iter+1;

    y0=[alfa;s0;0;1];
    
    [x,y]=ode45(g,[a b],y0,options);
     F=y(end,1)-beta;
    
    s1=s0-(F/(y(end,3)));
end
s=s1; % soluzione finale di s
