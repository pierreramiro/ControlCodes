function [Mo,observable]=getMo(A,C,show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Funcion: Halla la Matriz de observabilidad e indica si es observable
%Entradas
    %A,B   :Matrices del modelo de estados
%Salidas
    %Mo             :Matriz de observabilidad.
    %observable     :variable de tipo booleana que indica si es observable
    %                o no.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin<3
        show=true;
    end
    n=length(A);
    Mo=zeros(n);
    Mo(1,:)=C;
    for i=1:n-1
        Mo(i+1,:)=C*A^i;
    end
    observable=rank(Mo)==n;
    if show
        if observable
            fprintf('La matriz es observable =D y Mo es:\n')
            disp(Mo)
        else
            fprintf('La matriz no es observable')
        end
    end
end