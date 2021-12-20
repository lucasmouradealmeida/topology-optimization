%nome do arquivo MEF_chapa_elemento_retangular_ger.m
%programa de EF para chapas planas com elementos retangulares
%Baseado em Argyris 1954
%Prof Reyolando Brasil - abril/2018
%IC - Otimizacao Topologica de Chapas
%Orientador: Marcelo Araujo / Aluno: Lucas Moura

%VIGA ENGASTADA E EM BALANÇO - CARGA CONCENTRADA - NONLCON
%Gerador de malha Elemento Retangular

clc
clear

global compr alt ndx ndy nglpn nnel EM nu V;

%dimensoes do retangulo
compr=8;
alt=2;

%divisoes
ndx=16; ndy=4; %numero de divisoes em x e y

nel=ndx*ndy; %numero de elementos
nno=(ndx+1)*(ndy+1); %numero de nos
dx=compr/ndx; dy=alt/ndy; %tamanho para cada elemento

%tensao admissível do aço e espessura minima
csadm = 225e6;
ctmin = 0.01;
sadm = ones(nel,1);
sadm = sadm*csadm;
tmin = ones(nel,1);
tmin = tmin*ctmin;

%geracao das coordenadas dos nos
gcoord=zeros(nno,2); %gerador de coordenadas de cada elemento
x=0;y=0;k=0;
for i=1:ndx+1
    for j=1:ndy+1
        k=k+1;
        gcoord(k,1)=x; %coordenada x
        gcoord(k,2)=y; %coordenada y
        y=y+dy;
    end
    y=0;
    x=x+dx;
end

%geracao da conectividade dos elementos
nodel=zeros(nel,4); %matriz elementos finitos
kel=0;kaux=0;
for i=1:ndx
    for j=1:ndy
        kaux=kaux+1;
        kel=kel+1;
        nodel(kel,1)=kaux+ndy+2;
        nodel(kel,2)=kaux+1;
        nodel(kel,3)=kaux;
        nodel(kel,4)=kaux+ndy+1;
    end
    kaux=kaux+1;
end

%dimensao do problema
nglpn=2;%numero de graus de liberdade por no
nds=nno*nglpn; %numero de deslocamentos do sistema
nnel=4; %numero de nos por elemento
ndpel=nnel*nglpn;%numero de deslocamentos por elemento

% matriz de numero de graus de liberdade por no
LN=zeros(nno,nglpn);

%condicoes de contorno (n. de no restrito e direções restritas ou livres)
%engaste de viga em balanco
for i= k-ndy:k
    LN(i,:)=[-1 -1];
end

%dados fisicos dos elementos
%(esp., modulo elasticidade e coef.de Poisson)
t0=0.2; %espessura inicial (2x-0.2 3x-4 4x-?)
EM=200e9; %modulo de elasticidade do aço
nu=0.3; %coeficiente de Poisson do aço
EL=EM/(1-nu*nu);
G=EM/2/(1+nu); %modulo de cisalhamento
E=EL*[1 nu 0;nu 1 0;0 0 (1-nu)/2]; %Matriz elástica

%determinação das matrizes
%matriz LN
ngl=0;
for i=1:nno
    for j=1:nglpn
        if  LN(i,j)== 0
            ngl=ngl+1;
            LN(i,j)=ngl;
        end
    end
end
ngr=ngl;
for i=1:nno
    for j=1:nglpn
        if LN(i,j)<0
            ngr=ngr+1;
            LN(i,j)=ngr;
        end
    end
end

%inicializacao de matrizes e vetores
K=zeros(nds,nds);
p=zeros(nds,1);
P=zeros(nds,1);
Tens=zeros(nel,3);
q=zeros(8,1);


%Vetor Espessura
vt = ones(nel,1);  %vetor de espessuras inicializado igual a 1
vt = vt*t0;


%Area e Volume
Area = ones(1,nel);
Area = Area*dx*dy;
Vol = Area*vt;

%carga vertical na extremidade livre de viga em balanco
V = -100e4;

%FMINCON
disp('-----------------------------');
disp('   Algoritmo de Otimização   ');
disp('-----------------------------');

%Initial Guess
x0 = vt;

tic;
options = optimset('MaxIter',10000,'MaxFunEvals',1000000,'Display','iter','PlotFcn',{  @optimplotfunccount @optimplotfval @optimplotconstrviolation });
vol = @(vt) Area*vt;
[xopt,fval,exitflag,output] = fmincon(vol,x0,[],[],[],[],[],[],@nonlcon,options);
disp('Valor final função objetivo:');
disp(vol(xopt));
vt = xopt;
toc;

 h1 = gcf;
 h2 = figure(2); %change this variable for each figure
 copyobj(get(h1,'children'),h2);

% vetor de carregamento P
%carga vertical na extremidade livre de viga em balanco
%V=-100e4;
for i=1:ndy+1
    P(LN(i,2))=V/(ndy+1);
end

%introduzir aqui recalques de apoio
%matrizes de rigidez dos elementos
for iel=1:nel
    for j=1:nnel
        nd(j)=nodel(iel,j);
    end
    xa=gcoord(nd(1),1);xb=gcoord(nd(2),1);
    yb=gcoord(nd(2),2);yc=gcoord(nd(3),2);
    t = vt(iel,1);

%dimensoes do retangulo
    a=(xa-xb)/2;
    b=(yb-yc)/2;


%constantes
    c1=EL*t*b/3/a;
    c2=c1/2;
    c3=EL*t*nu/4;
    c4=G*t*a/3/b;
    c5=c4/2;
    c6=G*t/4;
%
    kd(1,1)=c1;kd(1,2)=c3;kd(1,3)=-c1;kd(1,4)=c3;kd(1,5)=-c2;kd(1,6)=-c3;kd(1,7)=c2;kd(1,8)=-c3;
    kd(2,2)=c1;kd(2,3)=-c3;kd(2,4)=c2;kd(2,5)=-c3;kd(2,6)=-c2;kd(2,7)=c3;kd(2,8)=-c1;
    kd(3,3)=c1;kd(3,4)=-c3;kd(3,5)=c2;kd(3,6)=c3;kd(3,7)=-c2;kd(3,8)=c3;
    kd(4,4)=c1;kd(4,5)=-c3;kd(4,6)=-c1;kd(4,7)=c3;kd(4,8)=-c2;
    kd(5,5)=c1;kd(5,6)=c3;kd(5,7)=-c1;kd(5,8)=c3;
    kd(6,6)=c1;kd(6,7)=-c3;kd(6,8)=c2;
    kd(7,7)=c1;kd(7,8)=-c3;
    kd(8,8)=c1;
%
    ks(1,1)=c4;ks(1,2)=c6;ks(1,3)=c5;ks(1,4)=-c6;ks(1,5)=-c5;ks(1,6)=-c6;ks(1,7)=-c4;ks(1,8)=c6;
    ks(2,2)=c4;ks(2,3)=c6;ks(2,4)=-c4;ks(2,5)=-c6;ks(2,6)=-c5;ks(2,7)=-c6;ks(2,8)=c5;
    ks(3,3)=c4;ks(3,4)=-c6;ks(3,5)=-c4;ks(3,6)=-c6;ks(3,7)=-c5;ks(3,8)=c6;
    ks(4,4)=c4;ks(4,5)=c6;ks(4,6)=c5;ks(4,7)=c6;ks(4,8)=-c5;
    ks(5,5)=c4;ks(5,6)=c6;ks(5,7)=c5;ks(5,8)=-c6;
    ks(6,6)=c4;ks(6,7)=c6;ks(6,8)=-c4;
    ks(7,7)=c4;ks(7,8)=-c6;
    ks(8,8)=c4;
%
    k=kd+ks;

%simetria
    for i=2:8
        for j=1:i-1
            k(i,j)=k(j,i);
        end
    end

%soma na matriz de rigidez do sistema
    kl=0;
    for n=1:nnel
        kl=kl+1;
        d(kl)=LN(nd(n),1);
        kl=kl+1;
        d(kl)=LN(nd(n),2);
    end
    for i=1:ndpel
        for j=1:ndpel
            K(d(i),d(j))=K(d(i),d(j))+k(i,j);
        end
    end
end

%Solucao do sistema
%calculo dos deslocamentos
disp('Deslocamentos')
p(1:ngl)=K(1:ngl,1:ngl)\(P(1:ngl)-K(1:ngl,ngl+1:nds)*p(ngl+1:nds));
disp(p)

%calculo das reações de apoio
disp('Esforços Nodais inclusive reacoes de apoio')
P(ngl+1:nds)=K(ngl+1:nds,1:ngl)*p(1:ngl)+K(ngl+1:ngl+1,ngl+1:ngl+1)...
    *p(ngl+1:nds);
disp(P)

%Tensões
disp('Tensoes no pto central dos elementos')
disp('sigma_x,    sigma_y,    tau_xy')
for iel=1:nel
    for j=1:nnel
        nd(j)=nodel(iel,j);
    end
    xa=gcoord(nd(1),1);xb=gcoord(nd(2),1);
    yb=gcoord(nd(2),2);yc=gcoord(nd(3),2);

%  dimensões do retangulo
    a=(xa-xb)/2;b=(yb-yc)/2;

%  constantes
    ca=1/4/a;cb=1/4/b;
    
%  matriz B=L*N calculada no centro do elemento x=y=0
    B=[ca 0 -ca 0 -ca 0 ca 0;0 cb 0 cb 0 -cb 0 -cb;cb ca cb -ca -cb -ca -cb ca];
%
    kl=0;
    for n=1:nnel
        kl=kl+1;
        d(kl)=LN(nd(n),1);
        kl=kl+1;
        d(kl)=LN(nd(n),2);
    end
    for i=1:ndpel
        q(i)=p(d(i));
    end
    tau=E*B*q;
    Tens(iel,:)=tau';
end
disp(Tens)

%Tensao normal (MPa)
disp('Tensão normal em cada elemento [MPa]')
disp('sigma_1,    sigma_2,    sigma_max')
smax = zeros(3,3);
for i=1:nel
    smax(i,1)=(Tens(i,1)+Tens(i,2))/2+(((Tens(i,1)-Tens(i,2))/2)^2+Tens(i,3)^2)^0.5;
    smax(i,2)=(Tens(i,1)+Tens(i,2))/2-(((Tens(i,1)-Tens(i,2))/2)^2+Tens(i,3)^2)^0.5;
    if abs(smax(i,1))>abs(smax(i,2))
       smax(i,3)=smax(i,1);
    else 
       smax(i,3)=smax(i,2);
    end
end
for i=1:nel
    smax(i,3)=abs(smax(i,3));
end
disp(smax)

%Vetor Tensão Máxima
disp('Tensão normal maxima dos elementos [MPa]')
vmax = zeros(nel,1);
for i=1:nel
    vmax(i)= smax(i,3);
end
disp(vmax);

%Generate the X and Y grid arrays using the MESHGRID function.
clf;
%TENSAO MAXIMA
x = (1:ndx+1);
y = (1:ndy+1);
[X,Y] = meshgrid(x,y);

%ESPESSURA
x1 = (1:ndx+1);
y1 = (1:ndy+1);
[X1,Y1] = meshgrid(x1,y1);

%Gerador da Matrix Z
Z = ones(ndy+1,ndx+1);
Z1 = ones(ndy+1,ndx+1);
kx=0;
for j=1:ndx
    for i=1:ndy
        kx=kx+1;
        Z(i,j)= vmax(kx);
        Z(ndy+1,j)= 0;
        Z(i,ndx+1)= 0;
        Z(ndy+1,ndx+1)=0;
        
        Z1(i,j)= vt(kx);
        Z1(ndy+1,j)= 0;
        Z1(i,ndx+1)= 0;
        Z1(ndy+1,ndx+1)=0;
    end
end

%Plots
 
%Tensão Máxima
subplot(2,1,1)
s = pcolor(X,Y,Z);
title('Tensão Máxima')
colormap(flipud(gray));
colorbar;
axis image;

%Espessura
subplot(2,1,2)
s1 = pcolor(X1,Y1,Z1);
title('Espessura')
colormap(flipud(gray));
colorbar;
axis image;
