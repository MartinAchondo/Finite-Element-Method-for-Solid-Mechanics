%Parametros:%
Area=250;
E=200e3;
matrizL=ones(4);
MatrizK=zeros(8);

%Elemento 1%
V1=[1,2];  % nodo 1 a 2%
Ang1=126.87; %Angulo respecto a x%
Le1=750;
Klocal1=((E*Area)/Le1)/1e5;
[ML1,l1,m1]=Matriztrans(Ang1,matrizL);
K1=Klocal1*ML1;
KT1=transfglobal(MatrizK,V1,K1); %Integración de matrix 4x4 a 8x8%
disp(['Para elemento 1: l=', num2str(l1), ' y m=', num2str(m1)]);


%Elemento 2%
V2=[1,3]; %nodo 1 a 3%
Ang2=36.87; %Angulo respecto a x%
Le2=1000;
Klocal2=((E*Area)/Le2)/1e5;
[ML2,l2,m2]=Matriztrans(Ang2,matrizL);
K2=Klocal2*ML2;
KT2=transfglobal(MatrizK,V2,K2);%Integración de matrix 4x4 a 8x8
disp(['Para elemento 2: l=', num2str(l2), ' y m=', num2str(m2)]);


%Elemento 
V3=[1,4]; %nodo 1 a 4;
Ang3=53; %Angulo respecto a x
13;
Le3=750;
Klocal3=((E*Area)/Le3)/1e5;
[ML3,l3,m3]=Matriztrans(Ang3,matrizL);
K3=Klocal3*ML3;
KT3=transfglobal(MatrizK,V3,K3); %Integración de matrix 4x4 a 8x8%
disp(['Para elemento 3: l=', num2str(l3), ' y m=', num2str(m3)]);

%Matriz de rigidez global%
KTransf=(KT1+KT2+KT3);


%Resolviendo sistema de ecuaciones reducido:%
A=(1e5)*[0.8,0.24;0.24,1.04];
B=[0;-18000];
Q=inv(A)*B;
disp(['Desplazamiento global: Q1=', num2str(Q(1)), ' y Q2=', num2str(Q(2))]);

%Calculo de esfuerzos:%
Sigma1=(E/Le1)*[-l1,-m1,l1,m1]*[0.0558;-0.1860;0;0];
Sigma2=(E/Le2)*[-l2,-m2,l2,m2]*[0.0558;-0.1860;0;0];
Sigma3=(E/Le3)*[-l3,-m3,l3,m3]*[0.0558;-0.1860;0;0];
disp(['Esfuerzo en los elemento: Sigmal=', num2str(Sigma1), ' , Sigma2=', num2str(Sigma2), ' y Sigma3=', num2str(Sigma3)]);

%Deformación unitaria%
e1=Sigma1*E;
e2=Sigma2*E;
e3=Sigma3*E;
disp(['Deformación unitaria: el=', num2str(e1), ' , e2=', num2str(e2), ' y e3=', num2str(e3)]);

%Fuerza en cada elemento:%
F1=Sigma1*Area;
F2=Sigma2*Area;
F3=Sigma3*Area;
disp(['Fuerza por elemento: Fl=', num2str(F1), ' , F2=', num2str(F2), ' y F3=', num2str(F3)]);

%reacciones%
C=[0.0558;-0.1860;0;0;0;0;0;0];
FF=[0;-18000;0;0;0;0;0;0];
R_prev=((1e5)*KTransf)*C-FF;
Reacciones=R_prev(3:8)

function [ML,ll,mm] =Matriztrans(Ang,matrizL)
    L=matrizL;
    l=cosd(Ang);
    m=sind(Ang);
    for i=1:4
        if i==1
            L(i,1)=l^2;
            L(i,2)=l*m;
            L(i,3)=-l^2;
            L(i,4)=-l*m;
        elseif i==2
            L(i,1)=l*m;
            L(i,2)=m^2;
            L(i,3)=-l*m;
            L(i,4)=-m^2 ;  
        elseif i==3
            L(i,1:4)=-L(1,1:4);
        elseif i==4
            L(i,1:4)=-L(2,1:4);
        end
    end
    ML=L;
    ll=l;
    mm=m;
end

function Trans=transfglobal(MatrizK,V,K)
    MM=MatrizK;
    i=V(1);
    j=V(2);
    if i==1 & j==2
        MM(1:2,1:2)=K(1:2,1:2);
        MM(1:2,3:4)=K(1:2,3:4);
        MM(3:4,1:2)=K(3:4,1:2);
        MM(3:4,3:4)=K(3:4,3:4);
    elseif i==1 & j==3
        MM(1:2,1:2)=K(1:2,1:2);
        MM(1:2,5:6)=K(1:2,3:4);
        MM(5:6,1:2)=K(3:4,1:2);
        MM(5:6,5:6)=K(3:4,3:4);
    elseif i==1 & j==4
        MM(1:2,1:2)=K(1:2,1:2);
        MM(1:2,7:8)=K(1:2,3:4);
        MM(7:8,1:2)=K(3:4,1:2);
        MM(7:8,7:8)=K(3:4,3:4);
    end
    Trans=MM;
end