function U=rk(f,a,b,A,b1,c,m,w0)
% DESCRIZIONE: implementa il metodo di Runge-Kutta di tableau A,b1,c
% INPUT: 
% f: funzione del problema
% a,b: estremi temporali
% A: matrice dei coefficienti del tableau
% b1: vettore dei coefficienti del tableau
% c: vettore dei coefficienti del tableau
% m: numero di nodi
% w0: soluzione iniziale al passo 0
%
% OUTPUT: 
% U: soluzione calcolata con il metodo di Runge-Kutta
nu=length(c); % numero di stadi
U = NaN(2,m); % vettore contenente la soluzione approssimata 
  U(:,1) = w0; 
  x0 = linspace(a,b,m);
  % applico il metodo definito dal tableau
  for n = 1:m-1
    h = (b-a)/(m-1); % ampiezza sottointervalli 
    K = zeros(2,nu);
   
    for i = 1:nu
      K(:,i) = feval(f,x0(n)+c(i)*h,U(:,n) +h* K*A(i,:)');
    end
    U(:,n+1) = U(:,n) + h*K*b1';
  end
  