function [s,v,iter]=shootingNewFree(f,df1,df2,alfa,beta,gamma,s0,b)
% DESCRIZIONE: trovare il tempo s risolvendo il problema ai 
% limiti a frontiera libera
% y''(t)=f(t,y(t),y'(t))     con t in (s,b)
% y(s)=alfa
% y'(s)=beta
% y(b)=gamma
%
% dove f,b,alfa,beta,gamma sono assegnati.
% Il metodo determina il valore del parametro s che soddisfa la seguente relazione:
%             F(s)=|y(b;s)-alfa,y'(b;s)-beta|-gamma<=toll 
% mediante il metodo di Newton. 
% INPUT:
% f = funzione del problema
% df1= derivata di f rispetto a y
% df2= derivata di f rispetto a y'
% alfa = valore della soluzione al tempo s 
% beta = valore della derivata della soluzione al tempo s
% gamma= valore della soluzione al tempo b
% s0= prima approssimazione di s
% a= estremo sinistro intervallo di integrazione (finale)
%
% OUTPUT:
% s = tempo finale da trovare con il metodo di Newton
% v = velocita' al tempo iniziale a 
% iter = numero di iterazioni necessarie 
options=odeset('AbsTol',1e-4);

g=@(t,y) [f(t,y); y(4); df1(y)*y(3)+df2(y)*y(4)]; % sistema ausiliario

g1=f(s0,[alfa;beta]);
y0=[alfa;beta;-beta;-g1(2)]; % condizione iniziale sistema ausiliario
[x,y]=ode45(g,[s0;b],y0,options); % risoluzione del sistema ausiliario con tempo 
                          % iniziale s0

F=y(end,1)-gamma; % F(s)
s1=s0-(F/(y(end,3)));


tol=1e-4; % tolleranza di ode45 
maxiter=100; % numero massimo iterazioni
iter=0; % contatore iterazioni

% procedo risolvendo il sistema ausiliario aggiornando la condizione
% iniziale con il nuovo s0 trovato con il metodo di Newton fino a che il
% metodo converge alla soluzione cercata 
while abs(s0-s1) > tol && iter<=maxiter
    s0=s1;
    iter=iter+1;
    g1=f(s0,[alfa;beta]);
    y0=[alfa;beta;-beta;-g1(2)]; % condizione iniziale sistema ausiliario
    [x,y]=ode45(g,[s0;b],y0,options);
     F=y(end,1)-gamma;
    
    s1=s0-(F/(y(end,3)));

end
s=s1; % tempo finale s
v=y(end,2); % velocita' iniziale 
