%% Set Matrices
A=[ 0 1 0;
    0 0 1;
    -4.4545 -31.5606    -0.1818];
B=[ 0;
    0;
    4.5455];
C=[ 1 0   0];
D=[0];
n=length(A);

%% Inicio de código
clear;
s=tf('s');
Gs=(s+2.5)/(s+2.5)/(s-1);
num=Gs.num{1};
den=Gs.den{1};
n=length(den)-1;



%% Canonica Controlable y Observable
Ac=[zeros(n-1,1) eye(n-1);
    flip(-den(2:end))];
Bc=[zeros(n-1,1);1];
Cc=[flip(num(2:end))];
Ao=Ac';
Bo=Cc';
Co=Bc';
Ac,Bc,Cc,Ao,Bo,Co
%%
Mc=[B];
Mo=[C];
for i=1:n-1
    Mc = [Mc,A^i*B];
    Mo = [Mo;C*A^i];
end
rank(Mc)
rank(Mo)
%% Canónica Diagonal o Modal
polos=roots(den)';
P=[ones(1,n)];
for i=1:n-1
    P(i+1,:)=polos.^i;
end
Ad=P^-1*Ac*P
Bd=P^-1*Bc
Cd=Cc*P
%% Canónica Jordan
[Tj,Aj]=jordan(Ac);
Aj
Bj=Tj^-1*Bc
Cj=Cc*Tj

