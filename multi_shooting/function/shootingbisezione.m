function [s,iter] = shootingbisezione(f,s0,s1,a,b,alfa,beta)
% DESCRIZIONE: Dato il problema di Cauchy
% y''(t)=f(t,y(t),y'(t))     con t in (a,b)
% y(a)=alfa
% y(b)=beta
%
% dove f,a,b,alfa,beta sono assegnati.
% Il metodo trasforma il problema ai limiti in uno di Cauchy del tipo:
%
% y1'(t)=y2(t)         con t in (a,b)
% y2'(t)=f(t,y1(t),y2(t))
% y1(a)=alfa
% y2(a)=s
%
% e determina il valore del parametro s che soddisfa la seguente relazione:
%             F(s)=|y2(b;s)-beta|<=toll 
% mediante il metodo di bisezione. 
%
% INPUT
% f=funzione del problema
% s0=approssimazione iniziale di s
% s1=approssimazione iniziale di s 
% a=tempo iniziale
% b=tempo finale 
% alfa=valore di y(a)
% beta=valore di y(b)
% 
% OUTPUT
% s valore di y'(a) trovato come zero della funzione F(s) 
% iter numero di iterazioni necessarie 
options=odeset('AbsTol',1e-4);

y0=[alfa; s0]; % condizioni iniziali approssimazione s0
y1=[alfa; s1]; % condizioni iniziale approssimazione s1
% soluzioni problemi con condizione iniziale y'(a)=s0 e y'(a)=s1
[x,Ys0]=ode45(f,[a,b],y0,options); 
[x1,Ys1]=ode45(f,[a,b],y1,options);

f0=Ys0(end,1)-beta; % F(s0)


f1=Ys1(end,1)-beta; % F(s1)

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
    
    
    ym=[alfa; sm]; 
    [xm,Ysm]=ode45(f,[a,b],ym,options);

    fm=Ysm(end,1)-beta;
    
    
    if sign(fm)*sign(f0)<=0
        s1=sm;
        f1=fm;
    else 
        s0=sm;
        f0=fm;
    end
    
end

s=(s0+s1)/2; % valore di s che soddisfa |F(s)|<tol

 
