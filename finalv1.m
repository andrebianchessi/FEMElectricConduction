clc
close all
global h e1 e2 es ei nELinha emax nNosLinha noMax nNos sa sb sc TESTE1 TESTE2
%% parametros
h=0.005; %delta x e delta y

xmax=0.32/2; %dimensoes da malha
ymax=0.3;


vi=0.1; %potencial da borda inferior
vs=220; %potencial da borda superior

p=0.4; %profundidade da barra

sa=5.9*10^7*p; %sigma a
sb=6.3*10^7*p; %sigma b
sc=4.7*10^7*p; %sigma c

TESTE1=0; %1 para teste de barra homogenea de ferro
TESTE2=1; %1 para teste de buraco em malha homogenea
TESTE3=0; %1 para teste com diferentes materiais

if TESTE3
    sa=5.9*10^7*p;
    sc=2.38*10^6*p;
    sb=10^2;
end

%% definicoes

%tipos de elementos
e1=1;%Ec1 ou Ee
e2=2;%Ec2 ou Es
es=3;%Es
ei=4;%Ei

nELinha=(xmax/h)*2; %numero de elem por linhas
emax=(xmax/h)*(ymax/h)*2; %ultimo index de elemento
nNosLinha=(xmax/h)+1;
noMax=nNosLinha*((ymax/h)+1);
nNos=noMax;

%% construcao da malha

nelementos=(xmax/h)*(ymax/h)*2; %numero de elementos
elementos=1:1:nelementos; %e=[1,2,...,nelementos]

n=1:1:nNos; %n=[1,2,...,imax]

%% vetores para criacao matrizes do sistema

nElementosT=nELinha*(ymax/h);%total de elem
nElementosBSI=nELinha;%total de elem na borda superior e inf

totalVars=noMax;
totalEqs=totalVars+nNosLinha*2;

A=sparse(totalEqs,totalVars);

J=zeros(totalEqs,1);
%% main
%eqs dos elementos
for e=elementos
    [u1,u2,u3]=n123(e);
    sig=sigma(e);
    tipo=tipoE(e);
    
    sig1=1*sig;
    sig5=sig*0.5;
    
    a=sparse(totalEqs,totalVars);%matriz local
    
    if ((tipo==e1)||(tipo==ei))
        %linha1
        a(u1,u1)=sig1;
        a(u1,u2)=-sig5;
        a(u1,u3)=-sig5;
        
        %linha2
        a(u2,u1)=-sig5;
        a(u2,u2)=sig5;
        %a(u2,u3)=0;
        
        %linha3
        a(u3,u1)=-sig5;
        %a(u3,u2)=0;
        a(u3,u3)=sig5;
    end
    
    if (tipo==e2||(tipo==es))
        %linha1
        a(u1,u1)=sig5;
        a(u1,u2)=-sig5;
        %a(u1,u3)=0;
        
        %linha2
        a(u2,u1)=-sig5;
        a(u2,u2)=sig1;
        a(u2,u3)=-sig5;
        
        %linha3
        %a(u3,u1)=0;
        a(u3,u2)=-sig5;
        a(u3,u3)=sig5;
    end
    A=A+a;
    
end
%condicoes de contorno
lg=noMax;
for i=1:nNosLinha%baixo
    lg=lg+1;
    A(lg,i)=1;
    J(lg)=vi;
end
for i=(noMax-nNosLinha+1):noMax%cima
    lg=lg+1;
    A(lg,i)=1;
    J(lg)=vs;
end

%remocao de linhas extras
nosLinhaInferior=zeros(nNosLinha,1);
nosLinhaSuperior=zeros(nNosLinha,1);

for i=1:(nNosLinha)
    nosLinhaInferior(i,1)=i;
end
j=1;
for i=(noMax-nNosLinha+1):noMax
    nosLinhaSuperior(j,1)=i;
    j=j+1;
end
delete=[nosLinhaInferior,nosLinhaSuperior];
A(delete, :) = [];
J(delete)=[];

%resolucao do sistema
U=A\J;

%% pos processamento

figure(1)
x=zeros(nNosLinha,1);
y=zeros(round(noMax/nNosLinha),1);
m=zeros(round(noMax/nNosLinha),nNosLinha);

for i=1:nNosLinha
    x(i,1)=xn(i);
end

l=1-nNosLinha;
for i=1:round(noMax/nNosLinha)
    l=l+nNosLinha;
    y(i,1)=yn(l);
end

for i=1:round(noMax/nNosLinha)%y
    for j=1:nNosLinha%x
        m(i,j)=U(no(i,j));
    end
end
surf(x,y,m);
colorbar
title(['Distribuição de potenciais na barra para h=' num2str(h) 'm'])

figure(2)

x=[];
y=[];
i=[];
j=[];
for el=1:2:emax
    [cxj,cyj,xj,yj]=vJ(el,U);
    x=[x,cxj];
    y=[y,cyj];
    i=[i,xj];
    j=[j,yj];
end
quiver(x,y,i,j);
title(['Vetor densidade de corrente para h=' num2str(h) 'm']);



%% calculo da resistencia do bloco
j=[];
for el=1:2:nELinha-1%borda inferior
    [cxj,cyj,xj,yj]=vJ(el,U);
    j=[j,yj];
end

integral=0;
for i=1:length(j)-1
    integral=integral+(j(i)+j(i+1))*h/2;
end
I=0.4*integral


R=-(vs-vi)/I
if TESTE1
    rho=9.71*10^(-8);
    l=ymax;
    area=p*xmax;
    R2=rho*l/area;
    erroPercentual=(R-R2)/R2*100
end


%%funcoes
function s=sigmaPto(x,y)
%retorna sigma do ponto x y
global sa sb sc TESTE1 TESTE2

s=sb;
if (y>0.26) && (x>0.06)%reg1
    s=sc;
end

if (x>0.1) && (y>0.22)%reg2
    s=sc;
end

if (x<0.14) && (y>0.18) && (y<0.22) %reg3
    s=sa;
end

if (x>0.08) && (x<0.14) && (y>0.12) && (y<0.18) %reg4
    s=sa;
end

if (x>0.04) && (y>0.08) && (y<0.12) %reg5
    s=sc;
end

if (x>0.04) && (x<0.1) && (y<0.08) %reg6
    s=sc;
end

if TESTE1
    s=1.03*10^7;
end

if TESTE2
    s=sa;
    if (x>0.02) && (x<0.14) && (y>0.14) && (y<0.16) %reg4
        s=0;
    end
end

end
function s=sigma(elemento)
%retorna sigma do elemento
[x,y]=xyCentro(elemento);
s=sigmaPto(x,y);
end
function [x,y]=xyCentro(elemento)
%retorna x e y do centro do elemento
global e1 e2 ei es
t=tipoE(elemento);
[n1,n2,n3]=n123(elemento);
x1=xn(n1);
x2=xn(n2);
x3=xn(n3);
y1=yn(n1);
y2=yn(n2);
y3=yn(n3);

if (t==e1 || t==ei)
    x=(x3+x2)/2;
    y=(y3+y2)/2;
end

if (t==e2 || t==es)
    x=(x3+x1)/2;
    y=(y3+y1)/2;
end
end
function r=tipoE(elemento)
%retorna o tipo do elemento de acordo com seu numero
global e1 e2 es ei nELinha emax
r=0;
if (elemento>emax || elemento<1) %checa index elemento
    disp('Erro indexacao elemento')
end
    
if (mod(elemento,2)==1)
    r=e1;
end
if (mod(elemento,2)==0)
r=e2;
end

if ((elemento<=nELinha)&& r==e1)
    r=ei;
end
if ((elemento>(emax-nELinha)) && r==e2)
    r=es;
end
end

function l=linhaE(elemento)
%retorna linha do elemento
global nELinha
if mod(elemento,nELinha)==0
    elemento=elemento-1;
end
l=fix(elemento/nELinha)+1;
end

function [n1,n2,n3]=n123(elemento)
%retorna index global dos nos 1 2 e 3 local
global e1 e2 es ei nELinha nNosLinha

i0=1+(linhaE(elemento)-1)*nNosLinha; %index do no da esquerda inf da linha
                                     %do elemento
t=tipoE(elemento);
if ((t==e1) || (t==ei))
    n1=i0+fix(mod(elemento,nELinha)/2);
    n2=n1+1;
    n3=n1+nNosLinha;    
end
if ((t==e2) || (t==es))
    n1=i0+fix(mod(elemento,nELinha)/2);
    
    if mod(elemento,nELinha)==0
        n1=i0+nNosLinha-1;
    end
    
    n2=n1+nNosLinha;
    n3=n2-1;    
end

end

function [cxj,cyj,xj,yj]=vJ(elemento,U)
%retorna componentes em x e y (xj,yj) e coord
% (cxj,cyj) da corrente no elemento e de tipo e1
global h nNosLinha
[n1,n2,n3]=n123(elemento);
n4=n2+nNosLinha;
x2=xn(n2);
x3=xn(n3);
y2=yn(n2);
y3=yn(n3);

s=sigma(elemento);


x=(x3+x2)/2;
y=(y3+y2)/2;

dvdx=((U(n4)+U(n2))/2-(U(n3)+U(n1))/2)/h;
dvdy=((U(n4)+U(n3))/2-(U(n1)+U(n2))/2)/h;

cxj=x;
cyj=y;
xj=-s*dvdx;
yj=-s*dvdy;

end

function x=xn(i)
%retorna coord x do no i
global h nNosLinha
x=h*(mod(i,nNosLinha)-1);
if mod(i,nNosLinha)==0
    x=h*(nNosLinha-1);
end
end

function y=yn(i)
%retorna coord x do no i
global h nNosLinha
if mod(i,nNosLinha)==0
    i=i-1;
end
y=h*fix(i/nNosLinha);
end

function n=no(i,j)
%retorna index do no da linha i coluna j da malha
global nNosLinha
i=i-1;
j=j-1;
n=1+j+i*nNosLinha;
end





