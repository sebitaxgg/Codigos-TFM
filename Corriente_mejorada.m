%tic
M =[0 1 -1; -1 0 +1 ; 1 -1 0];
Corrienteid = zeros(D/3);
CorrienteR = zeros(D/3);
CorrienteR2 = zeros(D/3);
Numerosopuestos = [2 3 1 ;3 1 2]; %Lo creo para que cuando j valga 1 tener el 2 y 3 y pongo el primero el que va con +i
f = @(x) exp((x-1)*(-3*x+10)*1i*pi/3);
for i = 1:D/3
    for l = 1:3
        for j = 1:3
            if B(j,Bsim(l,i))~=0
                c = 1;
                b = 1;
                Numeros = Numerosopuestos(:,j);
                k = i;
                while (k <= D/3 && (c+b)>0)  %c+b hacen de bandera para parar el bucle antes.
                    for n = 1:3
                        if B(Numeros(1),Bsim(n,k)) == B(Numeros(1),Bsim(l,i)) + 1 && B(j,Bsim(n,k)) == B(j,Bsim(l,i))-1
                            Corrienteid(i,k) = Corrienteid(i,k) + Corrsim(B,Bsim,l,n,i,k,M)/3;
                            CorrienteR(i,k) = CorrienteR(i,k) + conj(f(l))*f(n)*Corrsim(B,Bsim,l,n,i,k,M)/3;
                            CorrienteR2(i,k) = CorrienteR2(i,k) + f(l)*conj(f(n))*Corrsim(B,Bsim,l,n,i,k,M)/3;
                            if k~=i
                                Corrienteid(k,i) = conj(Corrienteid(i,k));%Porque la parte real es 0 asi que conj = -
                                CorrienteR(k,i) = conj(CorrienteR(i,k));
                                CorrienteR2(k,i) = conj(CorrienteR2(i,k));
                            end
                            c=0;
                        end
                        if B(Numeros(2),Bsim(n,k)) == B(Numeros(2),Bsim(l,i)) + 1 && B(j,Bsim(n,k)) == B(j,Bsim(l,i))-1
                            Corrienteid(i,k) = Corrienteid(i,k) + Corrsim(B,Bsim,l,n,i,k,M)/3;
                            CorrienteR(i,k) = CorrienteR(i,k) + conj(f(l))*f(n)*Corrsim(B,Bsim,l,n,i,k,M)/3;
                            CorrienteR2(i,k) = CorrienteR2(i,k) + f(l)*conj(f(n))*Corrsim(B,Bsim,l,n,i,k,M)/3;
                            if k~=i
                                Corrienteid(k,i) = conj(Corrienteid(i,k));%Porque la parte real es 0 asi que conj = -
                                CorrienteR(k,i) = conj(CorrienteR(i,k));
                                CorrienteR2(k,i) = conj(CorrienteR2(i,k));
                            end
                            b=0;
                        end 
                    end
                    k = k+1;
                end
            end
        end
    end
end
%Pasamos a la base de H,R la Corriente
operador2 = zeros(D);
operador2(1:D/3,1:D/3) = Vid'*Corrienteid*Vid/N;
operador2(1+D/3:2*D/3,1+D/3:2*D/3) = VR'*CorrienteR*VR/N;
operador2(1+2*D/3:D,1+2*D/3:D) = VR2'*CorrienteR2*VR2/N;
clear Corrienteid CorrienteR CorrienteR2

%toc
