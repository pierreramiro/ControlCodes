%Esta función permite hallar el LGR de un sistema siguiendo los
%9 pasos que se vieron en la clase

%Sintaxis: LGR(GH)  
%Entradas:  GH  --> es un objeto que contiene la Funciond de transferencia
%                   de lazo abierto del sistema. 
%                   Ejemplo: >>GH=tf([2 1],[1 3 1 1]);
%                            >>LGR(GH)
%Salidas: ninguna. El programa imprime los pasos a seguir

%Realizado por: Pierre Pérez. 
%Ante consultas --> pierreramiro@gmail.com

function LGR(GH)
    %Paso 1. Hallamos los polos. Punto incio k=0
    polos=roots(GH.den{1});
    fprintf('\nPaso 1. Hallamos los puntos de inicio\n')
    display (polos);
    %Paso 2. Obtenemos el nro de curvas
    m=length(polos);
    fprintf("\nPaso 2. Hay un total de %d curvas\n",m);
    %Paso 3. Hallamos los zeros. Puntos finales
    ceros=roots(GH.num{1});
    fprintf("\nPaso 3. Hallamos los puntos finales\n");
    display(ceros);
    n=length(ceros);
    fprintf("Habran %d asíntotas\n",m-n);
    %Paso 4. Ubicamos sobre eje real. En este caso usaremos rlocus
    fprintf("\nPaso 4. Hallamos los pts en el eje real\nVer grafica rlocus\n");
    rlocus(GH);
    %Paso 5. Hallamos los ángulos de las asíntotas
    angulos_asintotas=zeros(m-n,1);
    for i=1:(m-n)
        angulos_asintotas(i,1)=(2*i-1)*180/(m-n);   
    end
    fprintf("\nPaso 5. Hallamos los ángulos de las asíntotas\n");
    disp(angulos_asintotas);
    %Paso 6. Hallamos el CG
    if m-n>1
        fprintf("\nPaso 6. Hallamos el pto de intersección de las asíntotas\n");
        CG= (sum(polos)-sum(ceros))/(m-n);
        display(CG);
    else
        fprintf("\nPaso 6. No hay intersección, solo hay 1 asíntota\n");
    end
    %Paso 7. Desarrollamos arreglo Routh-Hurwitz
    fprintf("\nPaso 7. Hallamos los valores de K para estabilidad\n")
    syms k EPS;
    poly=GH.den{1}+ k*GH.num{1};
    sol=routh(poly,EPS);
    simplify(sol)
    %Paso 8. Hallamos los puntos de ruptura
    fprintf("\nPaso 8. Hallamos los puntos de ruptura\n");
    [q,~]=polyder(GH.den{1},GH.num{1});
    ptos_ruptura=roots(q);
    display(ptos_ruptura);
    %Paso 9. Hallamos los angulos de inicio
    if sum(abs(imag(polos)))~=0
    fprintf("\nPaso 9. Hallamos los angulos de inicio\n");
    conjugado_flag=0;
       for i=1:m
           if(imag(polos(i,1))~=0)&&conjugado_flag==0
               ang_polos=0;
               ang_ceros=0;
               temp_polo=polos(i,1);
               for ii=1:m       
                   if i~=ii
                       %analizamos los angulos +90 o -90
                       if (real(polos(ii,1))==real(temp_polo))
                           if imag(polos(ii,1)) > imag(temp_polo)
                               ang_polos=ang_polos-90;
                           else
                               ang_polos=ang_polos+90;
                           end
                       %analizamos los ángulos +180 o -180
                       elseif(imag(polos(ii,1))==imag(temp_polo))
                           if real(polos(ii,1)) > real(temp_polo)
                               ang_polos=ang_polos+180;
                           else
                               ang_polos=ang_polos+0;
                           end               
                       %analizamos los angulos que generan los polos del eje real
                       elseif (imag(polos(ii,1))==0)
                           CA=real(temp_polo)-polos(ii,1);%
                           CO=abs(imag(temp_polo));
                           if CA>0
                               ang_polos=ang_polos+atan(CO/CA)*180/pi;
                           else
                               ang_polos=ang_polos+180-atan(-CO/CA)*180/pi;
                           end               
                       %analizamos otros ángulos
                       else
                           CA=real(temp_polo)-real(polos(ii,1));
                           CO=imag(temp_polo)-imag(polos(ii,1));

                           if CO>0
                               if CA>0
                                   ang_polos=ang_polos+atan(CO/CA)*180/pi;
                               else
                                   ang_polos=ang_polos+180-atan(-CO/CA)*180/pi;
                               end
                           else
                               if CA>0
                                   ang_polos=ang_polos-atan(-CO/CA)*180/pi;
                               else
                                   ang_polos=ang_polos+180+atan(CO/CA)*180/pi;
                               end                       
                           end
                       end
                   end
               end
               for ii=1:n       
                   %analizamos los angulos +90 o -90
                   if (real(ceros(ii,1))==real(temp_polo))
                       if imag(ceros(ii,1)) > imag(temp_polo)
                           ang_ceros=ang_ceros-90;
                       else
                           ang_ceros=ang_ceros+90;
                       end
                   %analizamos los ángulos +180 o -180
                   elseif(imag(ceros(ii,1))==imag(temp_polo))
                       if imag(ceros(ii,1)) > imag(temp_polo)
                           ang_ceros=ang_ceros-90;
                       else
                           ang_ceros=ang_ceros+90;
                       end

                   %analizamos los angulos que generan los polos del eje real
                   elseif (imag(ceros(ii,1))==0)
                       CA=real(temp_polo)-ceros(ii,1);%
                       CO=abs(imag(temp_polo));
                       if CA>0
                           ang_ceros=ang_ceros+atan(CO/CA)*180/pi;
                       else
                           ang_ceros=ang_ceros+180-atan(-CO/CA)*180/pi;
                       end 
                   %analizamos otros ángulos
                   else
                       CA=real(temp_polo)-real(ceros(ii,1));
                       CO=imag(temp_polo)-imag(ceros(ii,1));

                       if CO>0
                           if CA>0
                               ang_ceros=ang_ceros+atan(CO/CA)*180/pi;
                           else
                               ang_ceros=ang_ceros+180-atan(-CO/CA)*180/pi;
                           end
                       else
                           if CA>0
                               ang_ceros=ang_ceros-atan(-CO/CA)*180/pi;
                           else
                               ang_ceros=ang_ceros+180+atan(CO/CA)*180/pi;
                           end                       
                       end
                   end
               end
               ang_polo=180+ang_ceros-ang_polos;
               fprintf("\nEl angulo del polo %f%+fi es: %f\n",real(temp_polo),imag(temp_polo),ang_polo);
               conjugado_flag=1; 
               continue
           end

           if conjugado_flag==1
               conjugado_flag=0;
           end
       end
    else
        fprintf("\nPaso 9.Ya son conocidos los angulos de inicio\n")%es redundante
    end
end