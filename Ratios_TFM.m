fichero_N1000 = 'ETH_Scaling/espectro.t1_mas.u-5.n1000.dat';
Energy = load(fichero_N1000);
Todos = Energy(:,3);
Energia_min = -3.9;
Energia_max = -3.85;
Autovalores = Todos(Todos > Energia_min & Todos < Energia_max);
num_Autovalores = length(Autovalores);
sn = zeros(num_Autovalores-1,1);
rn = zeros(num_Autovalores-1,1);
rn2 = zeros(num_Autovalores-1,1);
for i = 1:num_Autovalores-1
    sn(i) = Autovalores(i+1)-Autovalores(i);
    if i >1
        rn(i) = sn(i)/sn(i-1);
        rn2(i) = min(rn(i),1/rn(i));
    end
end
rn = rn(2:length(rn),1);
rn2 = rn2(2:length(rn2),1);
figure
h = histogram(rn2);
h.BinWidth = 1/35;
hold on
area = sum(h.Values)*h.BinWidth;
h.BinCounts = h.BinCounts/area;
f = @(x) 2./(x + 1).^2;
g = @(x) (27/4) * (x + x.^2) ./ (1 + x + x.^2).^(5/2);
modelo = @(theta,x) cos(theta)*f(x)+sin(theta)*g(x);
fplot(f,[0,1],'r')
fplot(g,[0,1],'g')
centros = h.BinEdges(1:end-1) + h.BinWidth/2;
valores = h.BinCounts;
legend('Histograma', 'Poisson', 'GOE')
% Mostrar ángulo y coeficientes
a = cos(theta_opt);
b = sin(theta_opt);
legend(' r̃','Poisson', 'GOE')
xlabel(' r̃')
title(['Ratios de los autoestados con r = 1,  s = 1 y energías ', num2str(Energia_min),'< e <', num2str(Energia_max)])
data = [centros(:), valores(:)];
save('ETH_Scaling/histograma_zona6.dat', 'data', '-ascii');