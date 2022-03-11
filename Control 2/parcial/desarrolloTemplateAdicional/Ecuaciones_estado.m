clear;clc;close all;
A=[ -3 -2;
    1   0];
B=[ 1;
    0];
C=[ 1   3];
D=[0];
n=length(A);

%% Valor y vector propio
[vector_propio,val_propio]=eig(A);
val_propio=diag(val_propio);
display(vector_propio)
display(val_propio)

%% Funcion de transferencia
s=tf('s');
Gs=C*(s*eye(n)-A)^-1*B;%+D;
polos=roots(Gs.den{1});%Tambien: polos=val_propio;

%% Solucion Ec. Estados
x0=[0.18 0];
t=0:0.01:10;
u = zeros(size(t));
[~,xt] = lsim(A,B,C,D,u,t,x0);
temp='plot(';
for i=1:n
    eval(sprintf('x%dt = xt(:,%d);',i,i));
    temp=strcat(temp,'t,x',num2str(i),'t,');
end
temp=strcat(temp(1:end-1),');');
eval(temp);

%Matriz exponencial: expm(A*t)
