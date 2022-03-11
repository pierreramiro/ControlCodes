%% Diseño de un controlador y observador
%1. verificamos la controlabilidad
%2. Establecemos los polos deseados y el polinomio caaracteristico del
%mismo: Diseño de polos por 2do orden o por Prototipo
%3. Hallamos los valores de K
    %Existen distintas formas para hallar K
        %Si está en la FCC: ki=alfai-ai
        %De lo contrario, transformarlo con T a la FCC y los valores de K
        %serán: K=Kc*T^-1. Siendo Kc los valores de ganancia en la FCC
        %Otra forma es usando Ackerman
%4. Diseñamos el observador. Analogamente, podemos calcular L si el sistema
%está en FCO.
    %De lo contrario, podemos calcular la matriz de transformación T y
    %hallar los valores de L=T*Lc
    %Otra forma es usando la formula de ackerman
%5. Adicional, se puede analizar la rpta en freq del sistema

%polos bessel
% 1:[-4.62];
% 2:[-4.0530+2.34*1j;-4.0530-2.34*1j];
% 3:[-5.0093;-3.9668+1j*3.7845;-3.9668-1j*3.7845];
% 4:[-4.0156+1j*5.0723;-4.0156-1j*5.0723;-5.5281+1j*1.6552;-5.5281-1j*1.6552];
% 5:[-6.4480;-4.1104+1j*6.3142;-4.1104-1j*6.3142;-5.9268+1j*3.0813;-5.9268-1j*3.0813];
% 6:[-4.2169+1j*7.53;-4.2169-1j*7.53;-6.2613+1j*4.4018;-6.2613-1j*4.4018;-7.1205+1j*1.4540;7.1205-1j*1.4540];
% 7:[-8.0271
% [-4.4554
% [-9.6585
% [-4.6835

%Funciones
% [Mc,controlable]=getMc(A,B,show)
% [Mo,observable]=getMo(A,C,show)
% [Ac,Bc,Cc,T,Mc,W]=getFCC(A,B,C,show)
% [Ao,Bo,Co,T,Mo,W]=getFCO(A,B,C,show)
% [K,Kc,T,Mc,W]=getKbyFCC(A,B,C,pol_des,show)
% [L,Lo,T,Mo,W]=getLbyFCO(A,B,C,pol_des,show)
% [K,poly_caract,Mc]=getKbyAcker(A,B,pol_des,show)
% [L,poly_caract,Mo]=getLbyAcker(A,C,pol_des,show)
% [K,L,Ki]=getServoBessel(A,B,C,Tes_c,Tes_o,show)
% [Gsys]=getFT(A,B,C,K,L,Bw)
%%


