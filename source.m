load data.mat
%Salvam datele din id si val si construim vectorii de timp:
u=id.InputData; y=id.OutputData; uv=val.InputData; yv=val.OutputData;
ts=id.Ts; t=(1:length(u))*ts; tv=(1:length(uv))*ts; 

nk= 2 ; %> alegem intarzierea = 1
M = 4 ; %> gradul maxim al polinomului y(k)
nn= 6 ; %> ordinul maxim ales pentru (na,nb)
N=length(y); NV=length(yv); %numarul de date din id si val
PHI=zeros(N,1); %> PHI-ul creat din datele de identificare,sub forma phi
phi=[];         %> un vector simbolic format din termenii polinomului y(k)
%                  de gradul m si cu toti coeficientii 1

%retinem valorile calculate in urmatoarele matrici:
MSE_pred=[]; Y_pred=[]; MSE_simu=[]; Y_simu=[]; NA=[]; NB=[]; MM=[];

%In continuare vom apela functia 'sym' pentru a crea un vector de variabile 
%simbolice de care ne vom folosi pentru crearea polinomului phi(k). 
%Deoarece numarul termenilor variaza in functie de na si nb, aceasta 
%functie se apeleaza de (nn!)*nn ori:

for nb=1:nn,for na=nb:nn %(na>=nb)
  y_k = 0; p=zeros(N,na+nb);
  v_sym=sym('v',[1 na+nb]);  
  for m=1:M
     %y(k) de gradul m si coeficientii obtinuti din ridicarea la putere:
      y_k = y_k+ expand(sum(v_sym(:))^m);  
      phi=[children(y_k) 1]; %> separam termenii, adaugam termenul liber
      PHI=zeros(N,length(phi));
      for i=1:length(phi)
          phi(i)=phi(i)/coeffs(phi(i)); %> coeficientii devin 1
      end
      
     %Transformam phi in functie anonima pentru a evita substitutia pe
     %simboluri, care determina un timp de executie mare:
      phimf=matlabFunction(phi);
     %Populam matricea PHI inlocuind in phimf cu datele din id 
      for k=2:N
          for i=1:min(na,k-1)
              p(k,i)=y(k-i);
          end
          for i=1:min(nb,k-nk)
              p(k,i+na)=u(k-nk-i+1);
          end
          args=num2cell(p(k,:));
          PHI(k,:)=feval(phimf, args{:});
      end
     %Obtinem vectorul de parametrii si initializam iesirile cu
     %coeficientul termenului liber:
      teta=PHI\y;
      y_pred=zeros(NV,1); y_pred(1)=teta(end); 
      y_simu=zeros(NV,1); y_simu(1)=teta(end);
      
     %Evaluam functia phimf pe datele de validare prin modul predictiei 
     %cu un pas inainte si al simularii: 
      for k=2:NV
          args=num2cell([yv(k-1:-1:max(1,k-na))' zeros(1,na-k+1)...
              uv(k-1:-1:max(1,k-nb))' zeros(1,nb-k+1)]);
          y_pred(k)=feval(phimf,args{:})*teta;
          args=num2cell([y_simu(k-1:-1:max(1,k-na))' zeros(1,na-k+1)...
              uv(k-1:-1:max(1,k-nb))' zeros(1,nb-k+1)]);
          y_simu(k)=feval(phimf,args{:})*teta;
      end
        
     %Calculam erorile medii patratice si memoram in matrici rezultatele:
      e=yv-y_pred; mse=mean(e.^2); MSE_pred=[MSE_pred mse];
      e=yv-y_simu; mse=mean(e.^2); MSE_simu=[MSE_simu mse];
      Y_pred=[Y_pred y_pred]; Y_simu=[Y_simu y_simu];
      NA=[NA na]; NB=[NB nb]; MM=[MM m];
    end
end,end
%Alegem cea mai mica eroare medie simulata si plotam iesirile aproximate:
[~,i]=min(MSE_simu);
figure(1); plot(tv,yv,'k',tv,Y_pred(:,i),'m',tv,Y_simu(:,i),'b'); 
title({['Eroare predictie ',num2str(MSE_pred(i))],...
       ['Eroare simulare ',num2str(MSE_simu(i))],...
       ['M=',num2str(MM(i)),' na=',num2str(NA(i)),' nb=',num2str(NB(i))]}); 
legend('y','y_{predictie}','y_{simulare}');

%Cream un tabel cu toate datele obtinute si il afisam:
TT=[];
for i=1:length(MM)
    col=[NA(i),NB(i),MM(i),MSE_pred(i)',MSE_simu(i)'];
    TT=[TT; col];
end
T=table(TT);
figure(2);
uitable('Data',T{:,:},'ColumnName',{'na','nb','M','MSE_pred','MSE_simu'},...
'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);