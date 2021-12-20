%IC - Otimizacao Topologica de Chapas
%Baseado no programa desenvolvido pelo Prof Reyolando Brasil
%Orientador: Marcelo Araujo / Aluno: Lucas Moura

%Algoritmo: KNITRO

clear
clc

disp('----------------Programa de Otimização de Chapas e Vigas----------------')
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%VARIAVEIS GLOBAIS%%%%%%%%%%%%%%

global Opt1 Mult EM nu csadm;

%%%%%%%%%%%%%%%%%%%%%%%%DEFAULT%%%%%%%%%%%%%%%%%%%%%%%%
%ORDEM: OT/OTVBL/OTVBU/OTVBC

Compr = [8;5;10;10];
Alt = [2;5;2;2];
Carreg = [-100e4;100e4;-100e4;-100e4];
Descrip = {'Chapa Engastada/Carregamento Concentrado na Extremidade';...
    'Chapa Biapoiada/Carregamento Concetrado Lateralmente';...
    'Viga Biapoiada/Carregamento Uniformemente Distribuido';...
    'Viga Biapoiada/Carregamento Concentrado no Centro'};


T1 = table(Compr,Alt,Carreg,Descrip,...
    'VariableNames',{'Comprimento','Altura','Carregamento','Descricao'},...
    'RowNames',{'(1) OT','(2) OTVBL','(3) OTVBU','(4) OTVBC'});
disp(T1)

Opt1 = input('Deseja qual processo de otimização ? ');

Mult = input('Qual multiplicador deseja ? ');

fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%MATERIAL%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------Material----------------')
fprintf('\n')

EM = 200e9; %Modulo de elasticidade 
nu = 0.3; %Coeficiente de Poisson

 T2 = table(EM,nu,...
     'VariableNames',{'ModuloElasticidade','CoeficientePoisson',},...
     'RowNames',{'Aço'});
 disp(T2)


fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%OTIMIZAÇÃO%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------Processo de Otimização----------------')
fprintf('\n')

otimizacao;

fprintf('\n\n')
%%%%%%%%%%%%%%%%%%%%%%%%CONFIABILIDADE%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------Processo de Confiabilidade----------------')
fprintf('\n')

%massa limite

porcent = input('Qual a porcentagem para massa limite ? ');

mel = (fval*7800)*(1+(porcent/100)); %Aumento em 10% do valor obtido [kg] %Modificar a porcentagem no teste de confiabilidade

Minter = input('Quantidade de iterações que deseja realizar ?  ');
me = ones(Minter,1);


tic;
for h=1:Minter
    csadm = norminv(rand(),260e6,26e6);
    fprintf('Rodando Confiabilidade %d ...',h);
    confiabilidade;
    me(h,1) = (fval*7800); %[kg]
end
time2=toc;

%M
M = mel - me;

%média
M_mean = mean(M);

%desvio padrao
M_std = std(M);

%beta
beta = (M_mean)/(M_std);

%função cumulativa normal padrão
fun2 = @(z) exp(-(((z).^2)/2));
int_fun2 = integral(fun2,-Inf,-beta);
P_phi = (1/(sqrt(2*pi)))*int_fun2;

T4 = table(beta,P_phi,time2,...
    'VariableNames',{'Beta','ProbabilidadeFalha','Tempo',},...
    'RowNames',{'Resultado (Confiabilidade)'});
disp(T4)

%Confiabilidade - Segundo Método

cont = 0;

for i = 1:Minter
    if(M(i)<0)
		cont = cont +1;
    end
end

falha2 = cont/Minter;
probabilidadefalha2 = falha2*100;

T5 = table(probabilidadefalha2,...
    'VariableNames',{'ProbabilidadeFalha2'},...
    'RowNames',{'Resultado (Confiabilidade2)'});
disp(T5)

