function [Gsys]=getFT(A,B,C,K,L,Bw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Funcion: Halla la Funcion de transferencia segun el numero de entradas
%Entradas
    %A,B,C   :Matrices del modelo de estados.
    %K       :Matriz de los coeficientes del regulador.
    %L       :Matriz de los coeficientes del observador.
    %Bw      :Matriz que multiplica a la pertubación cuando esta se
    %         considera como entrada.
%Salidas
    %Gsys    :Funcion de transferencia según las entradas.
    %           i)  Si las entradas son A,B,C       -> Gsys es la FT de la planta
    %           ii) Si las entradas son A,B,C,K     -> Gsys es la FT que permite hallar el Ess 
    %                                                  de la planta-regulador, debido a la 
    %                                                  perturbacion
    %           iii)Si las entradas son A,B,C,K,L   -> Gsys es la FT del regulador-observador
    %           iv) Si las entradas son A,B,C,K,L,Bw-> Gsys es la FT que permite hallar el Ess 
    %                                                  de la planta-regulador-observador
    %                                                  debido a la perturbación
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    n=length(A);
    switch nargin
        case 3
            fprintf('FT de la planta\n')
            syms s
            den=eval(coeffs(det(s*eye(n)-A),'s','all'));
            num=eval(coeffs(C*adjoint(s*eye(n)-A)*B,'s','all'));
        case 4
            fprintf('FT para hallar el Ess de la planta con regulador, debido a la perturbación\n')
            syms s
            temp=s*eye(n)-A+B*K;
            den=eval(coeffs(det(temp),'s','all'));
            num=eval(coeffs(C*adjoint(temp)*B,'s','all'));
        case 5
            fprintf('FT del regulador y observador\n')
            syms s
            temp=s*eye(n)-A+B*K+L*C;
            den=eval(coeffs(det(temp),'s','all'));
            num=eval(coeffs(-K*adjoint(temp)*L,'s','all'));
            num=-num;%ver diap. 28 y 30 de sem5. El signo se debe a la 
                     %realimentación negativa de 'y'
        case 6
            fprintf('FT de la planta, regulador y observador para hallar el Ess debido a w\n')
            syms s
            Awlc=[A-B*K B*K;
                  zeros(n) A-L*C];
            Bwlc=[Bw;Bw];
            Cwlc=[C C*0];
            temp=s*eye(2*n)-Awlc;
            den=eval(coeffs(det(temp),'s','all'));
            num=eval(coeffs(Cwlc*adjoint(temp)*Bwlc,'s','all'));    
        otherwise
            error('Error en las entradas')
            %error('escoja uno de los siguientes casos:\n1:FT para la planta\n2:FT para hallar el Ess de la planta con regulador\n3:FT del regulador y observador\n4:FT de la planta, regulador y observador')
    end
    Gsys=tf(num,den);
end