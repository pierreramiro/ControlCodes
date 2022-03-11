clear;clc;close all;
A=[ 0   1   0;
    -4 -0.4 0.8;
    0   0   -4];
B=[ 1;
    0;
    1.45];
C=[ 1   0   0];
n=length(A);

polos=roots(den);
P=[ones(1,n)];
for i=1:n-1
    P(i+1,:)=polos'.^i;
end
Ad=P^-1*A*P
Bd=P^-1*B
Cd=C*P