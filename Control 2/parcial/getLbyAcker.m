function [L,poly_caract,Mo]=getLbyAcker(A,C,pol_des,show)
if nargin<4
    show=true;
end
n=length(A);
[Mo,observable]=getMo(A,C);
if ~observable
    error('No es observable')
end
if show
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    fprintf('Procedemos a realizar el calculo de L por Ackerman\n')
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    fprintf('1. Debemos verificar si es observable hallando Mo y su rango\n')
    display(Mo)
    fprintf('Es observable =D\n')
end
syms s
poly_caract=poly(pol_des);
poly_caract=poly_caract/poly_caract(1);
beta=poly_caract(2:end);
if show
    fprintf('2. Hallamos el polinomio deseado para obtener phi:\n')
    disp(poly_caract)
end
phi=A^n;
if show
    fprintf('3. Mostramos a continuación los valores de [A^n A^(n-1) ... A^0]\n')
    disp(A^n)
    fprintf('\n')
end
for i=1:n
    phi=phi+A^(n-i)*beta(i);
    if show
        disp(A^(n-i))
        fprintf('\n')
    end
end
L=phi*Mo^-1*[zeros(1,n-1) 1]';
if show 
    fprintf('4. Mostramos el valor de phi:\n');
    display(phi)
    fprintf('5. Realizamos el calculo de L=phi*Mo^-1[0;0;...;1]\n')
    display(L)
end
end