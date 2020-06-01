function [s,v,iter]= shootingBisFree(f,alfa,beta,gamma,s0,s1,b)
% DESCRIZIONE: trovare il tempo s (finale) risolvendo il problema ai 
% limiti a frontiera libera
% y''(t)=f(t,y(t),y'(t))     con t in (s,b)
% y(s)=alfa
% y'(s)=beta
% y(b)=gamma
%
% dove f,b,alfa,beta,gamma sono assegnati.
% Il metodo determina il valore del parametro s che soddisfa la seguente relazione:
%             F(s)=|y(b;s)-alfa,y'(b;s)-beta|-gamma<=toll 
% mediante il metodo di bisezione. 
% INPUT:
% f = funzione del problema
% alfa = valore della soluzione al tempo s 
% beta = valore della derivata della soluzione al tempo s 
% gamma= valore della soluzione al tempo a 
% s0= prima approssimazione di s
% s1= seconda approssimazione di s
% a= estremo sinistro intervallo di integrazione 
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

% se la funzione non cambia segno nell'intervallo considerato non e'
% possibile utilizzare il metodo di bisezione e viene inviato un segnale di
% errore 
if sign(f0)'*sign(f1)>0
    error('intervallo non valido, la funzione non cambia di segno')
end

iter=0;

tol=1e-4; % tolleranza di ode45
maxiter=101;

% la funzione cambia segno quindi posso aggiornare s0 e s1 e applicare il
% metodo di bisezione per calcolare F(s)
while abs(s0-s1)>tol && iter<=maxiter
    iter=iter+1;
    sm=(s0+s1)/2;
    
    ym=[alfa; beta]; 
    [xm,Ysm]=ode45(f,[sm;b],ym,options);

    fm=Ysm(end,1)-gamma;
    
   % calcolo il nuovo s con il metodo di bisezione  
    if sign(fm)*sign(f0)<=0
        s1=sm;
        f1=fm;
    else 
        s0=sm;
        f0=fm;
    end
    
end

s=(s0+s1)/2; % valore di s che soddisfa |F(s)|<tol (tempo finale) 
v=Ysm(end,2); % velocita' iniziale (tempo t=s)