function [s,iter] = shootingSecanti(f,s0,s1,a,b,alfa,beta)
% DESCRIZIONE: Dato il problema di Cauchy
% y''(t)=f(t,y(t),y'(t))     con t in (a,b)
% y(a)=alfa
% y(b)=beta
%
% dove f,a,b,alfa, beta sono assegnati 
% Il metodo trasforma il problema ai limiti in uno di Cauchy del tipo:
%
% y1'(t)=y2(t)         con t in (a,b)
% y2'(t)=f(t,y1(t),y2(t))
% y1(a)=alfa
% y2(a)=s
%
% e determina il valore del parametro s che soddisfa la seguente relazione:
%             F(s)=|y2(b;s)-beta|<=tol
% mediante il metodo delle secanti.
%
% INPUT:
% f = funzione del problema di Cauchy
% s0 = approssimazione iniziale di s
% s1 = approssimazione iniziale di s
% a = estremo sinistro intervallo di integrazione
% b = estremo destro intervallo di integrazione 
% alfa = valore della soluzione nell'estremo sinistro
% beta = valore della soluzione nell'estremo destro
% 
%
% OUTPUT:
% s = valore di y'(a) 
% iter = numero di iterazioni necessarie
options=odeset('AbsTol',1e-4);

y0=[alfa; s0]; % condizioni iniziali approssimazione s0
y1=[alfa; s1]; % condizioni iniziale approssimazione s1

% soluzioni problemi con condizione iniziale y'(a)=s0 e y'(a)=s1
[x,Ys0]=ode45(f,[a,b],y0,options); 
[x1,Ys1]=ode45(f,[a,b],y1,options);

f0=Ys0(end,1)-beta; % F(s0)
f1=Ys1(end,1)-beta; % F(s1)



tol=1e-4; % tolleranza di ode45
maxiter=100;
iter=0;

while abs(s0-s1) > tol && iter <= maxiter
    iter=iter+1;
    s2=s1-(f1*(s1-s0)/(f1-f0)); % calcolo il nuovo valore di s con il 
                                % metodo delle secanti
    
    y2=[alfa; s2]; % condizioni iniziale approssimazione s2
    
    [x2,Ys2]=ode45(f,[a,b],y2,options);

    f0=f1; % F(s0)
    f1=Ys2(end,1)-beta; % F(s1) f1=f2
    s0=s1;
    s1=s2;
    
end
s=s2; % valore di s calcolato con il metodo delle secanti





















