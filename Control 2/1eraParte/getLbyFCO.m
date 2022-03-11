function [L,Lo,T,Mo,W]=getLbyFCO(A,B,C,pol_des,show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Funcion: Halla los valores de L por el método de transformación a la FCO y
%multiplicar por la matriz T. 
%Entradas
    %A,B,C   :Matrices del modelo de estados
    %pol_des :polos deseados del estimador
%Salidas
    %L  :Valores del regulador del modelo inicial
    %Lo :Valores del regulador del modelo en FCC
    %T  :La matriz de transformación
    %Mo :Matriz de observabilidad
    %W  :Matriz que permite hallar T=Mc*W
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin <5
        show=true;
    end
    [~,~,~,T,Mo,W]=getFCO(A,B,C);
    syms s
    n=length(A);
    poly_caract=eval(coeffs(det(s*eye(n)-A),'s','all'));
    poly_caract=poly_caract/poly_caract(1);
    a=poly_caract(2:end);
    poly_des=poly(pol_des);
    beta=poly_des(2:end);
    Lo=zeros(n,1);
    for i=1:n
        Lo(i,1)=beta(i)-a(i);
    end
    Lo=flip(Lo);
    L=T*Lo;
    if show
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        fprintf('Procedemos a realizar el calculo de L por FCO\n')
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        fprintf('1. Para poder obtener la FCO debemos ver si es observable hallando Mo y su rango:\n')
        display(Mo)
        fprintf('Es observable =D\n')
        fprintf('2. Obtenemos el polinomio característico:\n')
        disp(poly_caract)
        fprintf('3. Con los coeficientes del polinomio caracteristico hallamos la matriz W:\n')
        disp(W)
        fprintf('4. Realizamos la operación T=(W*Mo)^-1 para obtener T\n')
        display(T)
        fprintf('5. Hallamos el polinomio deseado alfa:\n')
        disp(poly_des)
        fprintf('6. Hallamos los valores de Lo=[beta(n)-a(n),...,beta(1)-a(1)')
        display(Lo);
        fprintf('7. Obtenemos L=T*Lo')
        display(L);
    end
end