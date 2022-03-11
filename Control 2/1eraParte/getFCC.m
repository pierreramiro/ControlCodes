function [Ac,Bc,Cc,T,Mc,W]=getFCC(A,B,C,show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Funcion: Permite obtener las matrices del modelo de estado en FCC, y
%otras matrices adicionales como la de transformación y controlabilidad
%Entradas
    %Matrices del modelo de estados: A,B,C
%Salidas
    %Ac,Bc,Cc   :Matrices del modelo de estado en FCC: 
    %T          :La matriz de transformación
    %Mc         :Matriz de controlabilidad
    %W          :Matriz que permite hallar T=Mc*W
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin <4
        show=true;
    end
    [Mc,controlable]=getMc(A,B,false);
    if ~controlable
        error('No es controlable no se puede hallar la FCC :(')
    end
    syms s
    n=length(A);
    poly_caract=eval(coeffs(det(s*eye(n)-A),'s','all'));
    poly_caract=poly_caract/poly_caract(1);
    W=zeros(n);
    W(1,:)=flip([1,poly_caract(2:end-1)]);
    for i=2:n
        W(i,:)=[W(i-1,2:end), 0];
    end
    T=Mc*W;
    Ac=T^-1*A*T;
    Bc=T^-1*B;
    Cc=C*T;
    if show
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        fprintf('Procedemos a obtener la FCC\n')
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        fprintf('1. Para poder obtener la FCC debemos ver si es controlable hallando Mc y su rango:\n')
        display(Mc)
        fprintf('es controlable\n')
        fprintf('2. Obtenemos el polinomio característico:\n')
        disp(poly_caract)
        fprintf('3. Con los coeficientes del polinomio caracteristico hallamos la matriz W:\n')
        disp(W)
        fprintf('4. Realizamos la operación T=Mc*W para obtener T\n')
        display(T)
        fprintf('5. Realizamos las operaciones Ac=T^-1*A*T\tBc=T^-1*B\tCc=C*T\nEl nuevo modelo de estados en FCC es el siguiente:\n')
        display(Ac);
        display(Bc);
        display(Cc);
    end
end