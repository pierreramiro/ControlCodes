function [K,L,Ki]=getServoBessel(A,B,C,Tes_c,Tes_o,show)
if nargin <6
    show=true;
end
%verificamos si el sistema es de tipo I
pol_des{1}=-4.62;
pol_des{2}=[-4.0530+2.34*1j;-4.0530-2.34*1j];
pol_des{3}=[-5.0093;-3.9668+1j*3.7845;-3.9668-1j*3.7845];
pol_des{4}=[-4.0156+1j*5.0723;-4.0156-1j*5.0723;-5.5281+1j*1.6552;-5.5281-1j*1.6552];
pol_des{5}=[-6.4480;-4.1104+1j*6.3142;-4.1104-1j*6.3142;-5.9268+1j*3.0813;-5.9268-1j*3.0813];
pol_des{6}=[-4.2169+1j*7.53;-4.2169-1j*7.53;-6.2613+1j*4.4018;-6.2613-1j*4.4018;-7.1205+1j*1.4540;7.1205-1j*1.4540];
        
n=length(A);
syms s
poly_caract=eval(coeffs(det(s*eye(n)-A),'s','all'));
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
fprintf('Procedemos a diseñar el Servosistema\n')
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
if poly_caract(end)==0
    %planta tipo I
    fprintf('La planta es de tipo I:Salidas[K,L]\n')
    K=getKbyAcker(A,B,pol_des{n}/Tes_c,show);
    Ki=[];
    L=getLbyAcker(A,C,pol_des{n}/Tes_o,show);
    if show
        fprintf('En resumen:\n')
        display(K)
        display(L)
    end
else
    %añadimos integrador
    fprintf('La planta no tiene integrador. Salidas:[K,L,Ki]\n')
    As=[A B*0;-C 0];
    Bs=[B;0];
    Cs=[C 0];
    if show
        fprintf('Se muestra el nuevo moledo de estados')
        display(As)
        display(Bs)
        display(Cs)
    end
    K=getKbyAcker(As,Bs,pol_des{n+1}/Tes_c,show);
    Ki=-K(end);
    K=K(1:end-1);
    L=getLbyAcker(A,C,pol_des{n}/Tes_o,show);
    if show
        fprintf('En resumen:\n')
        display(K)
        display(Ki)
        display(L)
    end
end
end