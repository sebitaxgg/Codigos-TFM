%tic
%N = 5; Número de partículas
U = -5;
L = 3; % Número de pozos
J0 = [1 1 1]; %J12 J13 J23
J = zeros(L);
J(1,2) = J0(1); J(2,1) = J0(1);
J(1,3) = J0(2); J(3,1) = J0(2);
J(2,3) = J0(3); J(3,2) = J0(3);
D = (N+1)*(N+2)/2;
%D = factorial(N+L-1)/(factorial(L-1)*factorial(N));
%Crear la base
B = zeros(L,(N+2)*(N+1)/2);
j=1;
for n1 = 0:N % Estos bucles funcionan para cualquier N y me crean la base pero solo para L = 3, pero se puede generalizar con una línea más
    for n2 = 0:(N-n1)
        B(:,j) = [n1 n2 N-n1-n2];
        j=j+1;
    end
end
%toc
%Para pasar de cano a auto V^(-1)
%Para pasar de auto a cano V
%Para pasar de cano a triple Cambiotrio^(-1)
%Para pasar de triple a cano Cambiotrio
%Para pasar de triple a eigtriple Vtrios^(-1)
%Para pasar de eigtriple a triple Vtrios