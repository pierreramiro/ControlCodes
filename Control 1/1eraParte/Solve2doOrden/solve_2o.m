%solve_2o
%Esta función halla los parámetros de un sistema de 2do orden
%usar de ejemplo Gs=tf(120,[1 9 120]);
%                solve_2o(tf(120,[1 9 120]))

%Entradas: Gs--> variable class tf. contiene la Fdt del sistema de 2do
%orden

%Salidas: info--> variable class cell. Contiene la siguiente información
%         info={new_Gs, gain,   ezeta,  wn, wd;
%               X,      tr,     ts,     n,  0};
%En la primera fila
%                 -new_Gs:  El nuevo arreglo de la funcion de transferencia
%                 -Gain:    La ganancia del sistema de 2do orden
%                 -ezeta:   el valor de ezeta
%                 -wn:      La frecuencia natural no amortiguada
%                 -wd:      La freq natural amortiguada
%En la segunda fila
%                 -X: el máximo sobre impulso con respecto al valor de
%                     establecimiento
%                 -tr: el tiempo de subida
%                 -ts: el tiempo de establecimiento
%                 -n: el numero de oscilaciones hasta llegar al valor
%                 estable

%Realizado por: Pierre Pérez
%email: pierreramiro@gmail.com
function [info]=solve_2o (Gs)
    step(Gs);
    grid;
    num=Gs.num{1};
    den=Gs.den{1};
    num=num(end)/den(1,1);
    den=den/den(1,1);
    gain=num/den(1,end);    
    wn=sqrt(den(1,end));
    data=find_parameters(den(1,end),den(1,end-1));%(wn_2,coef)
    ezeta=data(1);
    new_Gs=tf(wn*wn,[1 2*wn*ezeta wn*wn]);
    X=gain*calc_sobreimpulso(data);
    fprintf("\nEl máximo sobre impulso es:\n");
    disp(X);
    tr=calc_tr(data);
    fprintf("\nEl tiempo de subida es:\n");
    disp(tr);
    ts=calc_ts(data);
    fprintf("\nEl tiempo de establecimiento es:\n");
    disp(ts);
    n=calc_n(data);
    fprintf("\nEl nro de oscilaciones es:\n");
    disp(n);
    
    info={new_Gs,gain,ezeta,wn,data(3);X,tr,ts,n,0};
    
end

%la variable data tiene las sigte información=[ezeta;teta;wd;A1;A2;P1;P2];
function data=find_parameters (wn_2,coef)
    wn=sqrt(wn_2);
    ezeta=coef/(2*wn);
    A1=-(0.5) +ezeta/(2*sqrt(ezeta*ezeta-1));
    A2=-(0.5) -ezeta/(2*sqrt(ezeta*ezeta-1));
    P1=-ezeta*wn +wn*sqrt(ezeta*ezeta-1);
    P2=-ezeta*wn -wn*sqrt(ezeta*ezeta-1); 
    theta=0;
    wd=0;
    fprintf('\nEl sistema es:\t')
    if ezeta==1
        fprintf('criticamente amortiguado\n');
        display(ezeta);display(A1);display(A2);display(P1);display(P2);
    elseif ezeta<1
        fprintf('sub-amortiguado\n');
        theta=atan(sqrt(1-ezeta*ezeta)/ezeta);
        wd=wn*sqrt(1-ezeta*ezeta);    
    else %ezeta >1
        fprintf('sobre-amortiguado\n');
        display(ezeta);display(A1);display(A2);display(P1);display(P2);
    end
    display(ezeta);display(theta);display(wd);display(A1);display(A2);display(P1);display(P2);
    data=[ezeta;theta;wd;A1;A2;P1;P2];
end

function X=calc_sobreimpulso(data,n)
    arguments
        data=0.7
        n =1
    end
    ezeta=data(1);
%     teta=data(2);
%     wd=data(3);
%     A1=data(4);
%     A2=data(5);
%     P1=data(6);
%     P2=data(7);
    
    X=exp(-ezeta*pi*n/sqrt(1-ezeta*ezeta));
end

function tr=calc_tr(data)
%     ezeta=data(1);
     teta=data(2);
     wd=data(3);
%     A1=data(4);
%     A2=data(5);
%     P1=data(6);
%     P2=data(7);
    tr=(pi-teta)/wd;

end

function ts=calc_ts(data)
     ezeta=data(1);
%     teta=data(2);
     wd=data(3);
%     A1=data(4);
%     A2=data(5);
%     P1=data(6);
%     P2=data(7);
    wn=wd/sqrt(1-ezeta*ezeta);
    ts=4/(ezeta*wn);
end

function n=calc_n(data)
%      ezeta=data(1);
%     teta=data(2);
     wd=data(3);
%     A1=data(4);
%     A2=data(5);
%     P1=data(6);
%     P2=data(7);
    ts=calc_ts(data);
    n=floor(wd*ts/pi);
end
