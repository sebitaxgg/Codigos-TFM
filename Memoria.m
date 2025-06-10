Nexptrio = zeros(D);
for i = 1:D/3 %Para n1
    Nexptrio(D/3+i,i) = (B(1,Bsim(1,i)) + exp(2i*pi/3)*B(3,Bsim(1,i)) + exp(4i*pi/3)*B(2,Bsim(1,i)))/3;
    Nexptrio(2*D/3+i,i) = (B(1,Bsim(1,i)) + exp(4i*pi/3)*B(3,Bsim(1,i)) + exp(2i*pi/3)*B(2,Bsim(1,i)))/3;
    Nexptrio(2*D/3+i,D/3+i) = (B(1,Bsim(1,i)) + exp(2i*pi/3)*B(3,Bsim(1,i)) + exp(4i*pi/3)*B(2,Bsim(1,i)))/3;
    Nexptrio(i,D/3+i) = conj(Nexptrio(D/3+i,i));
    Nexptrio(i,2*D/3+i) = conj(Nexptrio(2*D/3+i,i)) ;
    Nexptrio(D/3+i,2*D/3+i) = conj(Nexptrio(2*D/3+i,D/3+i));
    Nexptrio(D/3+i,i) = Nexptrio(D/3+i,i)+ exp(2i*pi/3)*(B(2,Bsim(1,i)) + exp(2i*pi/3)*B(1,Bsim(1,i)) + exp(4i*pi/3)*B(3,Bsim(1,i)))/3;
    Nexptrio(2*D/3+i,i) = Nexptrio(2*D/3+i,i)+ exp(2i*pi/3)*(B(2,Bsim(1,i)) + exp(4i*pi/3)*B(1,Bsim(1,i)) + exp(2i*pi/3)*B(3,Bsim(1,i)))/3;
    Nexptrio(2*D/3+i,D/3+i) = Nexptrio(2*D/3+i,D/3+i)+ exp(2i*pi/3)*(B(2,Bsim(1,i)) + exp(2i*pi/3)*B(1,Bsim(1,i)) + exp(4i*pi/3)*B(3,Bsim(1,i)))/3;
    Nexptrio(i,D/3+i) = Nexptrio(i,D/3+i) + exp(2i*pi/3)*conj(B(2,Bsim(1,i)) + exp(2i*pi/3)*B(1,Bsim(1,i)) + exp(4i*pi/3)*B(3,Bsim(1,i)))/3;
    Nexptrio(i,2*D/3+i) = Nexptrio(i,2*D/3+i) + exp(2i*pi/3)*conj(B(2,Bsim(1,i)) + exp(4i*pi/3)*B(1,Bsim(1,i)) + exp(2i*pi/3)*B(3,Bsim(1,i)))/3;
    Nexptrio(D/3+i,2*D/3+i) = Nexptrio(D/3+i,2*D/3+i)+ exp(2i*pi/3)*conj(B(2,Bsim(1,i)) + exp(2i*pi/3)*B(1,Bsim(1,i)) + exp(4i*pi/3)*B(3,Bsim(1,i)))/3;
    Nexptrio(D/3+i,i) = Nexptrio(D/3+i,i) + exp(4i*pi/3)*(B(3,Bsim(1,i)) + exp(2i*pi/3)*B(2,Bsim(1,i)) + exp(4i*pi/3)*B(1,Bsim(1,i)))/3;
    Nexptrio(2*D/3+i,i) = Nexptrio(2*D/3+i,i) + exp(4i*pi/3)*(B(3,Bsim(1,i)) + exp(4i*pi/3)*B(2,Bsim(1,i)) + exp(2i*pi/3)*B(1,Bsim(1,i)))/3;
    Nexptrio(2*D/3+i,D/3+i) = Nexptrio(2*D/3+i,D/3+i)+exp(4i*pi/3)*(B(3,Bsim(1,i)) + exp(2i*pi/3)*B(2,Bsim(1,i)) + exp(4i*pi/3)*B(1,Bsim(1,i)))/3;
    Nexptrio(i,D/3+i) = Nexptrio(i,D/3+i)+ exp(4i*pi/3)*conj(B(3,Bsim(1,i)) + exp(2i*pi/3)*B(2,Bsim(1,i)) + exp(4i*pi/3)*B(1,Bsim(1,i)))/3;
    Nexptrio(i,2*D/3+i) = Nexptrio(i,2*D/3+i)+ exp(4i*pi/3)*conj(B(3,Bsim(1,i)) + exp(4i*pi/3)*B(2,Bsim(1,i)) + exp(2i*pi/3)*B(1,Bsim(1,i)))/3 ;
    Nexptrio(D/3+i,2*D/3+i) = Nexptrio(D/3+i,2*D/3+i)+ exp(4i*pi/3)*conj(B(3,Bsim(1,i)) + exp(2i*pi/3)*B(2,Bsim(1,i)) + exp(4i*pi/3)*B(1,Bsim(1,i)))/3;
end
%Pasamos a la base de H,R N1
operador1 = zeros(D);
operador1(1:D/3,1:D/3) = Vid'*Nexptrio(1:D/3,1:D/3)*Vid/N;
operador1(1+D/3:2*D/3,1+D/3:2*D/3) = VR'*Nexptrio(1+D/3:2*D/3,1+D/3:2*D/3)*VR/N;
operador1(1+2*D/3:D,1+2*D/3:D) = VR2'*Nexptrio(1+2*D/3:D,1+2*D/3:D)*VR2/N;
operador1(1:D/3,1+D/3:2*D/3) = Vid'*Nexptrio(1:D/3,1+D/3:2*D/3)*VR/N;
operador1(1+D/3:2*D/3,1:D/3) = VR'*Nexptrio(1+D/3:2*D/3,1:D/3)*Vid/N;
operador1(1:D/3,1+2*D/3:D) = Vid'*Nexptrio(1:D/3,1+2*D/3:D)*VR2/N;
operador1(1+2*D/3:D,1:D/3) = VR2'*Nexptrio(1+2*D/3:D,1:D/3)*Vid/N;
operador1(1+D/3:2*D/3,1+2*D/3:D) = VR'*Nexptrio(1+D/3:2*D/3,1+2*D/3:D)*VR2/N;
operador1(1+2*D/3:D,1+D/3:2*D/3) =  VR2'*Nexptrio(1+2*D/3:D,1+D/3:2*D/3)*VR/N;
clear Nexptrio Vid VR VR2 
%Pasamos a la base de H,R la Creac

%toc
