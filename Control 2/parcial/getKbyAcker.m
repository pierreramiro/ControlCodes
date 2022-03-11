function [K,poly_caract,Mc]=getKbyAcker(A,B,pol_des,show)
if nargin<4
    show=true;
end
n=length(A);
[Mc,controlable]=getMc(A,B,false);
if ~controlable
    error('No es controlable :(')
end
if show
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    fprintf('Procedemos a realizar el calculo de K por Ackerman\n')
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    fprintf('1. Debemos verificar si es controlable hallando Mc y su rango\n')
    display(Mc)
    fprintf('Es controlable =D\n')
end
syms s
poly_caract=poly(pol_des);
poly_caract=poly_caract/poly_caract(1);
alfa=poly_caract(2:end);
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
    phi=phi+A^(n-i)*alfa(i);
    if show
        disp(A^(n-i))
        fprintf('\n')
    end
end
K=[zeros(1,n-1) 1]*Mc^-1*phi;
if show 
    fprintf('4. Mostramos el valor de phi:\n');
    display(phi)
    fprintf('5. Realizamos el calculo de K=[0 0 ... 1]*Mc^-1*phi\n')
    display(K)
end
end