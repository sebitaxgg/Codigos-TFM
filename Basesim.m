%tic
Index = zeros(1,D);
Bsim = zeros(3, D/3); %Me dice los índices que vienen de rotaciones
j=1; %Me coloca los pares de vectores en Basesime
for i = 1:D
    if Index(i)==1
    else
        Bsim(1,j)=i;
        Vector = B(:,i);
        a=1;
        b=1;
        for k = i+1:D
            if Index(k)==1
            else
                if B(1,k)==Vector(3) && B(2,k) == Vector(1) && a==1
                    Bsim(2,j) = k;
                    Index(k)=1;
                    a=0;
                elseif B(1,k)==Vector(2) && B(3,k) == Vector(1) && b==1
                    Bsim(3,j) = k;
                    Index(k)=1;
                    b=0;
                end
            end
        end
        j= j+1;
    end 
end
clear Index
%toc
