%tic
Operador12trio = zeros(D);
Numerosopuestos = [2 3 1 ;3 1 2]; %Lo creo para que cuando j valga 1 tener el 2 y 3 y pongo el primero el que va con +i
f = @(x) exp((x-1)*(-3*x+10)*1i*pi/3);
for i = 1:D/3
    for l =1:3
        for j = 1:2
            if B(j,Bsim(l,i)) ~= 0
                c=1;
                b=1;
                Numeros = Numerosopuestos(:,j);
                k = i;
                while (k <= D/3 && (c+b)>0) %c+b hacen de bandera para parar el bucle antes.
                    for n = 1:3
                        if j ==1
                            if B(Numeros(1),Bsim(n,k)) == B(Numeros(1),Bsim(l,i)) + 1 && B(j,Bsim(n,k)) == B(j,Bsim(l,i))-1
                                Operador12trio(i,k) = Operador12trio(i,k)+ sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(1),Bsim(l,i))+1)/3;
                                Operador12trio(i+D/3,k+D/3) = Operador12trio(i+D/3,k+D/3)+sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(1),Bsim(l,i))+1)*conj(f(l))*f(n)/3;
                                Operador12trio(i+2*D/3,k+2*D/3) = Operador12trio(i+2*D/3,k+2*D/3) +sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(1),Bsim(l,i))+1)*conj(f(n))*f(l)/3;
                                Operador12trio(i,k+D/3) = Operador12trio(i,k+D/3)+ sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(1),Bsim(l,i))+1)*f(n)/3;
                                Operador12trio(i+D/3,k) = Operador12trio(i+D/3,k)+ sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(1),Bsim(l,i))+1)*conj(f(l))/3;
                                Operador12trio(i,k+2*D/3) = Operador12trio(i,k+2*D/3) +sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(1),Bsim(l,i))+1)*conj(f(n))/3;
                                Operador12trio(i+2*D/3,k) = Operador12trio(i+2*D/3,k) +sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(1),Bsim(l,i))+1)*f(l)/3;
                                Operador12trio(i+D/3,k+2*D/3) = Operador12trio(i+D/3,k+2*D/3) + sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(1),Bsim(l,i))+1)*conj(f(n))*conj(f(l))/3;
                                Operador12trio(i+2*D/3,k+D/3) = Operador12trio(i+2*D/3,k+D/3) + sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(1),Bsim(l,i))+1)*f(n)*f(l)/3;
                                if k ~= i
                                    Operador12trio(k,i) = conj(Operador12trio(i,k));
                                    Operador12trio(k+D/3,i+D/3) = conj(Operador12trio(i+D/3,k+D/3));
                                    Operador12trio(k+2*D/3,i+2*D/3) = conj(Operador12trio(i+2*D/3,k+2*D/3));
                                    Operador12trio(k+D/3,i) = conj(Operador12trio(i,k+D/3));
                                    Operador12trio(k+D/3,i+2*D/3) = conj(Operador12trio(i+2*D/3,k+D/3));
                                    Operador12trio(k+2*D/3,i+D/3) = conj(Operador12trio(i+D/3,k+2*D/3));
                                    Operador12trio(k,i+2*D/3) = conj(Operador12trio(i+2*D/3,k));
                                    Operador12trio(k+2*D/3,i) = conj(Operador12trio(i,k+2*D/3));
                                    Operador12trio(k,i+D/3) = conj(Operador12trio(i+D/3,k));
                                end
                                c=0;
                            end
                        end
                        if j == 2
                            if B(Numeros(2),Bsim(n,k)) == B(Numeros(2),Bsim(l,i)) + 1 && B(j,Bsim(n,k)) == B(j,Bsim(l,i))-1
                                Operador12trio(i,k) = Operador12trio(i,k)+ sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(2),Bsim(l,i))+1)/3;
                                Operador12trio(i+D/3,k+D/3) = Operador12trio(i+D/3,k+D/3)+  sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(2),Bsim(l,i))+1)*conj(f(l))*f(n)/3;
                                Operador12trio(i+2*D/3,k+2*D/3) = Operador12trio(i+2*D/3,k+2*D/3) + sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(2),Bsim(l,i))+1)*conj(f(n))*f(l)/3;
                                Operador12trio(i,k+D/3) = Operador12trio(i,k+D/3)+ sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(2),Bsim(l,i))+1)*f(n)/3;
                                Operador12trio(i+D/3,k) = Operador12trio(i+D/3,k)+ sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(2),Bsim(l,i))+1)*conj(f(l))/3;
                                Operador12trio(i,k+2*D/3) = Operador12trio(i,k+2*D/3) + sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(2),Bsim(l,i))+1)*conj(f(n))/3;
                                Operador12trio(i+2*D/3,k) = Operador12trio(i+2*D/3,k) +sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(2),Bsim(l,i))+1)*f(l)/3;
                                Operador12trio(i+D/3,k+2*D/3) = Operador12trio(i+D/3,k+2*D/3) + sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(2),Bsim(l,i))+1)*conj(f(n))*conj(f(l))/3;
                                Operador12trio(i+2*D/3,k+D/3) = Operador12trio(i+2*D/3,k+D/3) + sqrt(B(j,Bsim(l,i)))*sqrt(B(Numeros(2),Bsim(l,i))+1)*f(n)*f(l)/3;
                                if k ~= i
                                     Operador12trio(k,i) = conj(Operador12trio(i,k));
                                    Operador12trio(k+D/3,i+D/3) = conj(Operador12trio(i+D/3,k+D/3));
                                    Operador12trio(k+2*D/3,i+2*D/3) = conj(Operador12trio(i+2*D/3,k+2*D/3));
                                    Operador12trio(k+D/3,i) = conj(Operador12trio(i,k+D/3));
                                    Operador12trio(k+D/3,i+2*D/3) = conj(Operador12trio(i+2*D/3,k+D/3));
                                    Operador12trio(k+2*D/3,i+D/3) = conj(Operador12trio(i+D/3,k+2*D/3));
                                    Operador12trio(k,i+2*D/3) = conj(Operador12trio(i+2*D/3,k));
                                    Operador12trio(k+2*D/3,i) = conj(Operador12trio(i,k+2*D/3));
                                    Operador12trio(k,i+D/3) = conj(Operador12trio(i+D/3,k));
                                end
                                b=0;
                                
                            end 
                        end
                    end
                    k = k+1;
                end
            end
        end
    end
end
operador3= zeros(D);
operador3(1:D/3,1:D/3) = Vid'*Operador12trio(1:D/3,1:D/3)*Vid/N;
operador3(1+D/3:2*D/3,1+D/3:2*D/3) = VR'*Operador12trio(1+D/3:2*D/3,1+D/3:2*D/3)*VR/N;
operador3(1+2*D/3:D,1+2*D/3:D) = VR2'*Operador12trio(1+2*D/3:D,1+2*D/3:D)*VR2/N;
operador3(1:D/3,1+D/3:2*D/3) = Vid'*Operador12trio(1:D/3,1+D/3:2*D/3)*VR/N;
operador3(1+D/3:2*D/3,1:D/3) = operador3(1:D/3,1+D/3:2*D/3)';
operador3(1:D/3,1+2*D/3:D) = Vid'*Operador12trio(1:D/3,1+2*D/3:D)*VR2/N;
operador3(1+2*D/3:D,1:D/3) = operador3(1:D/3,1+2*D/3:D)';
operador3(1+D/3:2*D/3,1+2*D/3:D) = VR'*Operador12trio(1+D/3:2*D/3,1+2*D/3:D)*VR2/N;
operador3(1+2*D/3:D,1+D/3:2*D/3) = operador3(1+D/3:2*D/3,1+2*D/3:D)';
clear  Operador12trio 
%toc
