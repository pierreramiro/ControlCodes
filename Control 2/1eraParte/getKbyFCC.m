function [K,Kc,T,Mc,W]=getKbyFCC(A,B,C,pol_des,show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Funcion: Halla los valores de K por el método de transformación a la FCC y
%multiplicar por la matriz T. 
%Entradas
    %A,B,C   :Matrices del modelo de estados
    %pol_des :polos deseados del regulador
%Salidas
    %K  :Valores del regulador del modelo inicial
    %Kc :Valores del regulador del modelo en FCC
    %T  :La matriz de transformación
    %Mc :Matriz de controlabilidad
    %W  :Matriz que permite hallar T=Mc*W
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin <5
        show=true;
    end
    [~,~,~,T,Mc,W]=getFCC(A,B,C,false);
    syms s
    n=length(A);
    poly_caract=eval(coeffs(det(s*eye(n)-A),'s','all'));
    poly_caract=poly_caract/poly_caract(1);
    a=poly_caract(2:end);
    poly_des=poly(pol_des);
    alfa=poly_des(2:end);
    Kc=zeros(1,n);
    for i=1:n
        Kc(1,i)=alfa(i)-a(i);
    end
    Kc=flip(Kc);
    K=Kc*T^-1;
    if show
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        fprintf('Procedemos a realizar el calculo de K por FCC\n')
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        fprintf('1. Para poder obtener la FCC debemos ver si es controlable hallando Mc y su rango:\n')
        display(Mc)
        fprintf('Es controlable =D\n')
        fprintf('2. Obtenemos el polinomio característico:\n')
        disp(poly_caract)
        fprintf('3. Con los coeficientes del polinomio caracteristico hallamos la matriz W:\n')
        disp(W)
        fprintf('4. Realizamos la operación T=Mc*W para obtener T\n')
        display(T)
        fprintf('5. Hallamos el polinomio deseado alfa:\n')
        disp(poly_des)
        fprintf('6. Hallamos los valores de Kc=[alfa(n)-a(n),...,alfa(1)-a(1)')
        display(Kc);
        fprintf('7. Obtenemos K=Kc*T^-1')
        display(K);
    end
end