function s=solveMulti(f,a,b,alfa,beta,m) 
% DESCRIZIONE: funzione che risolve il metodo di shooting o shooting
% multiplo con la funzione fsolve di MATLAB. 
% 
% INPUT: 
% f=funzione del problema
% a=estremo sx intervallo temporale 
% b=estremo dx intervallo temporale 
% alfa=valore di y(a) 
% beta=valore di y(b) 
% m=numero di nodi (1=metodo di shooting semplice, m>1 metodo di shooting
% multiplo) 
%
% OUTPUT: 
% s=vettore [alfa,s1,alfa2,s2,...,beta,sm] contenente il valore della
% soluzione e della sua derivata ad ogni nodo

% Uso un'approssimazione iniziale della soluzione data dalla retta passante
% per i punti (a,alfa) e (b,beta) 
y1=@(x)alfa+(beta-alfa)*(x-a)/(b-a);

x1=linspace(a,b,m);
% trovo un'approssimazione dei valori della soluzione nei nodi x_k
alfa1=y1(x1);
% trovo un'approssimazione dei vaolori della derivata prima della soluzione
% nei nodi x_k
s0=ones(1,m)*(beta-alfa)/(b-a);

% creo una prima approssimazione del vettore s
sx(1:2:2*m)=alfa1;
sx(2:2:2*m)=s0;

% trovo F
F=@(s)fsolveShootNew(f,s,a,b,alfa,beta,m); 
% trovo s attraverso la funzione fsolve di MATLAB
s=fsolve(F,sx);

end