%tic
%Elementos de Matriz
Hid = zeros(D/3,D/3); %Hamiltoniano identidad
HR = zeros(D/3,D/3);%Hamiltoniano rotación
HR2 = zeros(D/3,D/3);%Hamiltoniano 2 rotación
%Puede estar mal HR y HR2 y sean las traspuestas?
Numerosopuestos = [2 3 1 ;3 1 2];
f = @(x) exp((x-1)*(-3*x+10)*1i*pi/3);
for i = 1:D/3
    for l =1:3
        Hid(i,i) = Hid(i,i) + Hsim(B,Bsim,U,N,l,l,i,i)/3;
        HR(i,i) = HR(i,i) + Hsim(B,Bsim,U,N,l,l,i,i)/3;
        HR2(i,i) = HR2(i,i) + Hsim(B,Bsim,U,N,l,l,i,i)/3;
        for j =1:3
            if (B(j,Bsim(l,i)) ~= 0)
               b = 1;
                c = 1;
                Numeros = Numerosopuestos(:,j);
                k = i;
                while (k <= D/3 && (c+b)>0)
                    for n = 1:3
                        if B(Numeros(1),Bsim(n,k)) == B(Numeros(1),Bsim(l,i)) + 1 && B(j,Bsim(n,k)) == B(j,Bsim(l,i))-1
                                Hid(i,k) = Hid(i,k) + Hsim(B,Bsim,U,N,l,n,i,k)/3;
                                HR(i,k) = HR(i,k) + conj(f(l))*f(n)*Hsim(B,Bsim,U,N,l,n,i,k)/3;
                                HR2(i,k) = HR2(i,k) + f(l)*conj(f(n))*Hsim(B,Bsim,U,N,l,n,i,k)/3;
                                if k ~= i
                                    Hid(k,i) = Hid(i,k);
                                    HR(k,i) = conj(HR(i,k));
                                    HR2(k,i) = conj(HR2(i,k));
                                end %Porque la parte real es 0 asi que conj = -
                                c=0;
                        end
                        if B(Numeros(2),Bsim(n,k)) == B(Numeros(2),Bsim(l,i)) + 1 && B(j,Bsim(n,k)) == B(j,Bsim(l,i))-1
                            Hid(i,k) = Hid(i,k) + Hsim(B,Bsim,U,N,l,n,i,k)/3;
                            HR(i,k) = HR(i,k) + conj(f(l))*f(n)*Hsim(B,Bsim,U,N,l,n,i,k)/3;
                            HR2(i,k) = HR2(i,k) + f(l)*conj(f(n))*Hsim(B,Bsim,U,N,l,n,i,k)/3;
                            if k ~= i
                                Hid(k,i) = Hid(i,k);
                                HR(k,i) = conj(HR(i,k));
                                HR2(k,i) = conj(HR2(i,k));
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
HR(D/3,D/3) = real(HR(D/3,D/3));
HR2(D/3,D/3) = real(HR2(D/3,D/3));
clear H
[Vid,Did] = eig(Hid);
eid = diag(Did);
%invVid = Wid';
clear Did Hid
[VR,DR] = eig(HR);
eR = diag(DR);
%invVR = WR';
clear DR HR
[VR2,DR2] = eig(HR2);
eR2 = diag(DR2);
%invVR2 = WR2';
clear DR2 HR2
%toc
