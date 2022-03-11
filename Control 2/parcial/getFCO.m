function [Ao,Bo,Co,T,Mo,W]=getFCO(A,B,C,show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Funcion: Permite obtener las matrices del modelo de estado en FCO, y
%otras matrices adicionales como la de transformación y obserrvabilidad
%Entradas
    %Matrices del modelo de estados: A,B,C
%Salidas
    %Ao,Bo,Co   :Matrices del modelo de estado en FCO 
    %T          :La matriz de transformación
    %Mo         :Matriz de observabilidad
    %W          :Matriz que permite hallar T=(W*Mo)^-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin<4
        show=true;
    end
    [Mo,observable]=getMo(A,C,false);
    if ~observable
        error('No es observable no se puede hallar la FCC :(')
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
    T=(W*Mo)^-1;
    Ao=T^-1*A*T;
    Bo=T^-1*B;
    Co=C*T;
    if show
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        fprintf('Procedemos a obtener la FCO\n')
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
        fprintf('5. Realizamos las operaciones Ac=T^-1*A*T\tBc=T^-1*B\tCc=C*T\nEl nuevo modelo de estados en FCO es el siguiente:\n')
        display(Ao);
        display(Bo);
        display(Co);
    end
end