%Código para hacer un scaling sobre traza y autovalores en los subespacios
%de misma energía
headerlines=0;
parts=[20 28 40 56 80 112 160 224 320];
%parts=[241];
for num = 1:length(parts) %Bucle para controlar las N, las N irán desde 3*1+1 hasta 3*(último número más uno)
    tic
	    N = parts(num); %La fórmula para los N, a lo mejor es muy exhaustiva yendo a casi todos los casos
    fprintf("N = %d\n", N);


    Base % Creo la Base "canónica" de nuestro sistema
    Basesim % Creo la Base de Rotación
    fprintf("Energia\n")
    Energiasim_mejorada % Creo la matriz H
    Corriente_mejorada % Genero la corriente en la base "canónica"

    Datos_a_guardar = zeros(1,10);
    filename = sprintf('ETH_Corriente_%d.txt', N);
    if exist(filename, 'file') == 2
        delete(filename);  % Elimina el archivo si existe para empezar desde cero
    end

fid = fopen(filename, 'w');
    for i = 1:D/3
        Datos_a_guardar(1,1:2) = [eid(i)/N, eR(i)/N]; % Guardo las energías
        matrizrot2 = zeros(2,2);
        matrizrot2(1,1) = operador2(D/3+i,D/3+i);
        matrizrot2(2,2) = operador2(2*D/3+i,2*D/3+i);
        matrizrot2(1,2) = operador2(i+D/3,i+2*D/3);
        matrizrot1(2,1) = operador2(2*D/3+i,i+D/3);
        Datos_a_guardar(1,3:5) = [trace(matrizrot2), eig(matrizrot2)'];
        if abs(eid(i)-eR(i))<10^(-6)
            matriztriple2 = zeros(3,3);
            matriztriple2(1,1) = operador2(i,i);
            matriztriple2(2,2) = operador2(D/3+i,D/3+i);
            matriztriple2(3,3) = operador2(2*D/3+i,2*D/3+i);
            matriztriple2(1,2) = operador2(i,i+D/3);
            matriztriple2(2,1) = operador2(i+D/3,i);
            matriztriple2(1,3) = operador2(i,2*D/3+i);
            matriztriple2(3,1) = operador2(2*D/3+i,i);
            matriztriple2(2,3) = operador2(i+D/3,i+2*D/3);
            matriztriple2(3,2) = operador2(2*D/3+i,i+D/3);
            Datos_a_guardar(1,6:9) = [trace(matriztriple2), eig(matriztriple2)'];
            Datos_a_guardar(10)=3;
        else
            Datos_a_guardar(1,6:9) = [trace(matrizrot2), eig(matrizrot2)', operador2(i,i)];
            Datos_a_guardar(10)=2;
        end
        save(filename,'Datos_a_guardar','-ascii','-append');
    end

    clear operador2
      fclose(fid);
    
    CreacDest12_Mejorado % Genero el a1(dagger)a2+a2(dagger)a1 en la base "canónica

      Datos_a_guardar = zeros(1,10);
filename = sprintf('ETH_Hopping_%d.txt', N);
    if exist(filename, 'file') == 2
        delete(filename);  % Elimina el archivo si existe para empezar desde cero
    end

fid = fopen(filename, 'w');

    for i = 1:D/3
        Datos_a_guardar(1,1:2) = [eid(i)/N, eR(i)/N]; % Guardo las energías
        matrizrot3 = zeros(2,2);
        matrizrot3(1,1) = operador3(D/3+i,D/3+i);
        matrizrot3(2,2) = operador3(2*D/3+i,2*D/3+i);
        matrizrot3(1,2) = operador3(i+D/3,i+2*D/3);
        matrizrot3(2,1) = operador3(2*D/3+i,i+D/3);
        Datos_a_guardar(1,3:5) = [trace(matrizrot3), eig(matrizrot3)'];
        if abs(eid(i)-eR(i))<10^(-6)
            matriztriple3 = zeros(3,3);
            matriztriple3(1,1) = operador3(i,i);
            matriztriple3(2,2) = operador3(D/3+i,D/3+i);
            matriztriple3(3,3) = operador3(2*D/3+i,2*D/3+i);
            matriztriple3(1,2) = operador3(i,i+D/3);
            matriztriple3(2,1) = operador3(i+D/3,i);
            matriztriple3(1,3) = operador3(i,2*D/3+i);
            matriztriple3(3,1) = operador3(2*D/3+i,i);
            matriztriple3(2,3) = operador3(i+D/3,i+2*D/3);
            matriztriple3(3,2) = operador3(2*D/3+i,i+D/3);
            Datos_a_guardar(1,6:9) = [trace(matriztriple3), eig(matriztriple3)'];
            Datos_a_guardar(10)=3;
        else
            Datos_a_guardar(1,6:9) = [trace(matrizrot3), eig(matrizrot3)', operador3(i,i)];
            Datos_a_guardar(10)=2;
        end
        save(filename,'Datos_a_guardar','-ascii','-append');
    end

    clear operador3
      fclose(fid);
      
    Memoria
    clear B Bsim
    %Matriz_Ni_autotrios_mejorada % Genero las matrices N1,N2,N3 en la base de rotación y paso todo a la base conjunta de rotación y hamiltoniano

      Datos_a_guardar = zeros(1,15);
filename = sprintf('ETH_Rotacion_%d.txt', N);
    if exist(filename, 'file') == 2
        delete(filename);  % Elimina el archivo si existe para empezar desde cero
    end

fid = fopen(filename, 'w');

    for i = 1:D/3
        Datos_a_guardar(1,1:2) = [eid(i)/N, eR(i)/N]; % Guardo las energías
        matrizrot1 = zeros(2,2);
        matrizrot1(1,1) = operador1(D/3+i,D/3+i);
        matrizrot1(2,2) = operador1(2*D/3+i,2*D/3+i);
        matrizrot1(1,2) = operador1(i+D/3,i+2*D/3);
        matrizrot1(2,1) = operador1(2*D/3+i,i+D/3);
        Datos_a_guardar(1,3:7) = [trace(matrizrot1), real(eig(matrizrot1))', imag(eig(matrizrot1))'];   
        if abs(eid(i)-eR(i))<10^(-6)
            matriztriple1 = zeros(3,3);
            matriztriple1(1,1) = operador1(i,i);
            matriztriple1(2,2) = operador1(D/3+i,D/3+i);
            matriztriple1(3,3) = operador1(2*D/3+i,2*D/3+i);
            matriztriple1(1,2) = operador1(i,i+D/3);
            matriztriple1(2,1) = operador1(i+D/3,i);
            matriztriple1(1,3) = operador1(i,2*D/3+i);
            matriztriple1(3,1) = operador1(2*D/3+i,i);
            matriztriple1(2,3) = operador1(i+D/3,i+2*D/3);
            matriztriple1(3,2) = operador1(2*D/3+i,i+D/3);
            Datos_a_guardar(1,8:14) = [trace(matriztriple1), real(eig(matriztriple1))', imag(eig(matriztriple1))'];   
            Datos_a_guardar(15)=3;
        else
	    Datos_a_guardar(1,8:14) = [trace(matrizrot1), real(eig(matrizrot1))',imag(eig(matrizrot1))', real(operador1(i,i)), imag(operador1(i,i))];   
            Datos_a_guardar(15)=2;
        end
        save(filename,'Datos_a_guardar','-ascii','-append');
    end
    clear operador1 eid eR
    fclose(fid);
    toc
end
