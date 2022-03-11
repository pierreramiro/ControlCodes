function [Mc,controlable]=getMc(A,B,show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Funcion: Halla la Matriz de controlabilidad e indica si es controlable
%Entradas
    %A,B   :Matrices del modelo de estados
%Salidas
    %Mc             :Matriz de controlabilidad.
    %controlable    :variable de tipo booleana que indica si es controlable 
    %                o no.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin<3
        show=true;
    end
    n=length(A);
    Mc=zeros(n);
    Mc(:,1)=B;
    for i=1:n-1
        Mc(:,i+1)=A^i*B;
    end
    controlable=rank(Mc)==n;
    if show
        if controlable
            fprintf('La matriz es controlable =D y Mc es:\n')
            disp(Mc)
        else
            fprintf('La matriz no es controlable')
        end
    end
end