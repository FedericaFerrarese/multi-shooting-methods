function F=fsolveShootNew(f,sx1,a,b,alfa,beta,m) 
% DESCRIZIONE: questa funzione crea la funzione
% F=@(s) [y(x(k+1);x(k);s(k))-s(k+1); y'(x(k+1);x(k);s(k))] valutata in sx1
%
% INPUT: 
% f=funzione del problema
% sx1= vettore in cui F dovrà essere valutata
% a= estremo sx intervallo temporale
% b=estremo dx intervallo temporale
% alfa=valore di y(a) 
% beta=valore di y(b)
% m=numero di nodi in cui si vuole dividere l'intervallo [a,b]
%
% OUTPUT: 
% F=funzione F=@(s) [y(x(k+1);x(k);s(k))-s(k+1); y'(x(k+1);x(k);s(k))] 
% valutata in sx1

 
 y0=[sx1(1);sx1(2)]; % condizione iniziale sistema ausiliario
 x1=linspace(a,b,m);
        
 [x,Y]=ode45(f,[a b],y0);
 F(1)=Y(1,1)-alfa;
 F(2)=Y(end,1)-beta;
        
        j=2;
 for k=1:m-1
        y0=[sx1(j-1);sx1(j)];

        [x,y]=ode45(f,[x1(k) x1(k+1)],y0); %  risoluzione del sistema ausiliario
        F1=y(end,1)-sx1(j+1); % F(s)
        F2=y(end,2)-sx1(j+2);
        F(j+1:j+2)=[F1;F2];
        j=j+2;
    end
  end
