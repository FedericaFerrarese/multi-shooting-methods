function [s,v,iter]= shootingSecFree(f,s0,s1,alfa,beta,gamma,b)
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
% mediante il metodo delle secanti. 
% INPUT:
% f = funzione del problema
% s0= approssimazione iniziale di s
% s1= approssimazione iniziale di s 
% alfa = valore della soluzione al tempo s 
% beta = valore della derivata della soluzione al tempo s 
% gamma= valore della soluzione al tempo b
% b= estremo destro intervallo di integrazione
%
% OUTPUT:
% s = tempo iniziale
% v = velocita' al tempo iniziale s 
% iter = numero di iterazioni necessarie 
options=odeset('AbsTol',1e-4);

y0=[alfa; beta]; % condizioni iniziali approssimazione s0
y1=[alfa; beta]; % condizioni iniziale approssimazione s1

% soluzioni problemi con tempo iniziale s0 e s1
[x,Ys0]=ode45(f,[s0;b],y0,options); 
[x1,Ys1]=ode45(f,[s1;b],y1,options);

f0=Ys0(end,1)-gamma; % F(s0)
f1=Ys1(end,1)-gamma; % F(s1)


tol=1e-4; % tolleranza assoluta di ode45
maxiter=100;
iter=0;

while abs(s0-s1) > tol && iter <= maxiter
    iter=iter+1;
    s2=s1-(f1*(s1-s0)/(f1-f0)); % calcolo il nuovo valore di s con il 
                                % metodo delle secanti
    
    y2=[alfa; beta]; % condizioni iniziale approssimazione s2
    
    [x2,Ys2]=ode45(f,[s2;b],y2,options);

    f0=f1; % F(s0)
    f1=Ys2(end,1)-gamma; % F(s1) f1=f2
    s0=s1;
    s1=s2;
    
end
s=s2; % tempo iniziale 
v=Ys2(end,2); % velocita' al tempo iniziale 