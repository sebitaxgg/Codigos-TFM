ValoresN = [20 28 40 56 80 112 160 224 320];
%Correr esto
colores = [
    0.000, 0.447, 0.741;  % azul (por defecto de MATLAB)
    0.850, 0.325, 0.098;  % naranja
    0.929, 0.694, 0.125;  % amarillo
    0.494, 0.184, 0.556;  % morado
    0.466, 0.674, 0.188;  % verde
    0.301, 0.745, 0.933;  % celeste
    0.635, 0.078, 0.184;  % rojo oscuro
    0.25,  0.25,  0.25;   % gris oscuro
    0.75,  0.75,  0;      % amarillo oscuro
    0.75,  0,     0.75;   % magenta fuerte
    0,     0.5,   0.5;    % verde-azulado
    1.0,   0.6,   0.0;    % naranja fuerte
    0.6,   0.4,   0.8;    % lavanda
    0.8,   0.2,   0.2;    % rojo claro
    0.2,   0.8,   0.2;    % verde claro
    0.2,   0.2,   0.8;    % azul claro
    0.5,   0.5,   0.5;    % gris medio
    0.8,   0.6,   0.7;    % rosa viejo
    0.4,   0.6,   0.2;    % verde oliva
    0.1,   0.1,   0.1;    % casi negro
    0.9,   0.4,   0.7     % fucsia
];
%ValoresN = 20;
fichero_hopping = 'ETH_Scaling/ETH_Hopping_';
fichero_autoval_hopping = 'ETH_Scaling/ETH_Hopping_Autoval_Full_';
fichero_autoval_corriente = 'ETH_Scaling/Corriente_Positivos_';
fichero_autoval_rotacion = 'ETH_Scaling/Orden_';
fichero_micro = 'ETH_Scaling/micro_hopping.dat';
Micro = load(fichero_micro);
energia_min = -3.1;
energia_max = -2.3;
%Cotas_autos = [0.2,0.1,0.05,0.025];
%Cotas_Hopp = [0.08, 0.04, 0.02,0.01];
Cotas_autos=zeros(1,6);
Cotas_Hopp=zeros(1,6);
leyenda_Auto = cell(1,length(Cotas_autos));
leyenda_Hopp = cell(1,length(Cotas_autos));
for i = 1:length(Cotas_autos)
    Cotas_autos(i) = 0.08 + 0.02*i; %Para auto Hopping funciona a partir de 0.05
    leyenda_Auto{i} = ['Cota = ' num2str(Cotas_autos(i))];
    Cotas_Hopp(i)= 0.015 + 0.005*i;
    leyenda_Hopp{i} = ['Cota = ' num2str(Cotas_Hopp(i))];
end


Subespacios = zeros(length(ValoresN),length(Cotas_autos));
Subespacios_Term_Corr = zeros(length(ValoresN),length(Cotas_autos));
Probabilidad_Term_Corr = zeros(length(ValoresN),length(Cotas_autos));
Subespacios_Term_Rot = zeros(length(ValoresN),length(Cotas_autos));
Probabilidad_Term_Rot = zeros(length(ValoresN),length(Cotas_autos));
Subespacios_Term_Traza_Hopp = zeros(length(ValoresN),length(Cotas_autos));
Probabilidad_Term_Traza_Hopp = zeros(length(ValoresN),length(Cotas_autos));
Subespacios_Term_Auto_Hopp = zeros(length(ValoresN),length(Cotas_autos));
Probabilidad_Term_Auto_Hopp = zeros(length(ValoresN),length(Cotas_autos));
%Maximo_Autoval_Hopping = zeros(length(ValoresN),length(Cotas_autos));
for i = 1:length(ValoresN)  
    N = ValoresN(i);
    Fichero_hopping = fichero_hopping + string(ValoresN(i))+'.txt'; 
    Fichero_autoval_hopping = fichero_autoval_hopping + string(ValoresN(i))+'.dat'; 
    Fichero_autoval_corriente = fichero_autoval_corriente + string(ValoresN(i))+'.dat'; 
    Fichero_autoval_rotacion = fichero_autoval_rotacion + string(ValoresN(i))+'.dat'; 
    Hopping = load(Fichero_hopping);
    Autoval_Hopping = load(Fichero_autoval_hopping);
    Corriente = load(Fichero_autoval_corriente);
    Rotacion = load(Fichero_autoval_rotacion);
    %Corriente
    Eid_Corriente = Corriente(:,1);
    Auto_Corriente = Corriente(:,2);
    Auto_Rot = abs(Rotacion(:,2));
    Traza_Hopp = Hopping(:,6);
    Micro_interpolado = interp1(Micro(:,1),Micro(:,3),Eid_Corriente,'linear');
    Residuo_Traza = abs(Micro_interpolado-Traza_Hopp/3);
    Autovalor_Hopp1 = Hopping(:,7);
    Autovalor_Hopp2 = Hopping(:,9);
    Autovalor_Hopp3 = Hopping(:,9);
    Residuo_Auto_Hopp1 = abs(Micro_interpolado-Autovalor_Hopp1);
    Residuo_Auto_Hopp2 = abs(Micro_interpolado-Autovalor_Hopp2);
    Residuo_Auto_Hopp3 = abs(Micro_interpolado-Autovalor_Hopp3);
    for j = 1:length(Cotas_autos)
        %Estudiar si autoval Corriente cercano a 0
        for k = 1:length(Eid_Corriente)
            if(Eid_Corriente(k)<energia_max && Eid_Corriente(k)>energia_min)
                Subespacios(i,j) = Subespacios(i,j)+1;
                if (Auto_Corriente(k) > Cotas_autos(j))
                    Subespacios_Term_Corr(i,j) = Subespacios_Term_Corr(i,j) + 1;
                end
                if (Auto_Rot(k) > Cotas_autos(j))
                    Subespacios_Term_Rot(i,j) = Subespacios_Term_Rot(i,j) + 1;
                end
                if (Residuo_Traza(k) > Cotas_Hopp(j))
                    Subespacios_Term_Traza_Hopp(i,j) = Subespacios_Term_Traza_Hopp(i,j) + 1;
                end
                Maximo_Autoval_Hopping = max([Residuo_Auto_Hopp1(k),Residuo_Auto_Hopp2(k),Residuo_Auto_Hopp3(k)]);
                if (Maximo_Autoval_Hopping > Cotas_autos(j))
                    Subespacios_Term_Auto_Hopp(i,j) = Subespacios_Term_Auto_Hopp(i,j) + 1;
                end

                %{
                if (Residuo_Auto_Hopp1(k) < Cotas_Hopp(j))
                    Subespacios_Term_Auto_Hopp(i,j) = Subespacios_Term_Auto_Hopp(i,j) + 1;
                end
                if (Residuo_Auto_Hopp2(k) < Cotas_Hopp(j))
                    Subespacios_Term_Auto_Hopp(i,j) = Subespacios_Term_Auto_Hopp(i,j) + 1;
                end
                if (Residuo_Auto_Hopp3(k) < Cotas_Hopp(j))
                    Subespacios_Term_Auto_Hopp(i,j) = Subespacios_Term_Auto_Hopp(i,j) + 1;
                end
                %}
            end
        end
        Probabilidad_Term_Corr(i,j) = Subespacios_Term_Corr(i,j)/Subespacios(i,j);
        Probabilidad_Term_Rot(i,j) = Subespacios_Term_Rot(i,j)/Subespacios(i,j);
        Probabilidad_Term_Traza_Hopp(i,j) = Subespacios_Term_Traza_Hopp(i,j)/Subespacios(i,j);
        Probabilidad_Term_Auto_Hopp(i,j) = Subespacios_Term_Auto_Hopp(i,j)/(Subespacios(i,j));
    end
    %{
    disp(N)
    disp(Subespacios_Term_Auto_Hopp(i,:))
    disp(3*Subespacios(i,:))
    figure
    grid on
    hold on
    %scatter(Eid_Corriente,Residuo_Auto_Hopp1,10,'r','filled')
    %scatter(Eid_Corriente,Residuo_Auto_Hopp2,10,'g','filled')
    scatter(Eid_Corriente,Residuo_Auto_Hopp3,10,'b','filled')
    xlim([energia_min energia_max]);
    ylim([0 0.2]);
    yline(0.08,'k')
    yline(0.04,'r')
    yline(0.02,'b')
    yline(0.01,'g')
    title(N);  % Agrega el título
    hold off
    %}
end

%colores = ['r','g','b','m'];
%Gráfica Corriente
figure
hold on
grid on
for i = 1:length(Cotas_autos)
    plot(ValoresN, Subespacios_Term_Corr(:,i)', '-', 'Color', colores(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    scatter(ValoresN,Subespacios_Term_Corr(:,i)',10,colores(i,:),'filled')
end
%ylim([0 1.05]);
legend(leyenda_Auto)
title('Subespacios por encima cota Auto Corriente');
xlabel('N')
ylabel('# Subespacios')
hold off

%Gráfica Rotación
figure
hold on
grid on
for i = 1:length(Cotas_autos)
    plot(ValoresN, Subespacios_Term_Rot(:,i)', '-', 'Color', colores(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');

    scatter(ValoresN,Subespacios_Term_Rot(:,i)',10,colores(i,:),'filled')
end
%ylim([0 1.05]);
legend(leyenda_Auto)
title('Subespacios por encima cota Auto Rotacion');
xlabel('N')
ylabel('# Subespacios')
hold off

%Gráfica Traza 
figure
hold on
grid on
for i = 1:length(Cotas_autos)
    plot(ValoresN, Subespacios_Term_Traza_Hopp(:,i)', '-', 'Color', colores(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');

    scatter(ValoresN, Subespacios_Term_Traza_Hopp(:,i)',10,colores(i,:),'filled')
end
%ylim([0 1.05]);
legend(leyenda_Hopp)
title('Subespacios por encima cota Traza Hopping');
xlabel('N')
ylabel('# Subespacios')
hold off

%Gráfica Auto Hopp
figure
hold on
grid on
for i = 1:length(Cotas_autos)
    plot(ValoresN, Subespacios_Term_Auto_Hopp(:,i)', '-', 'Color', colores(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');

    scatter(ValoresN,Subespacios_Term_Auto_Hopp(:,i)',10,colores(i,:),'filled')
end
%ylim([0 1.05]);
legend(leyenda_Auto)
title('Subespacios por encima cota Auto Hopping');
xlabel('N')
ylabel('# Subespacios')
hold off

