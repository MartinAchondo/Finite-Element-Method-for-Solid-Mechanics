format compact
clc

if ~exist('results', 'dir')
    mkdir('results');
end

% Inputs del problema %
P = 50e3;
L = 500e-3;
r_0 = 50e-3;

% Propiedades del material %
rho = 7850;
g = 9.81;
gamma = rho*g;
E = 200e9;

A_0 = pi*r_0^2;
sigma_0 = P/A_0;

prop = [rho,gamma,E];
prob = [P,L,r_0,A_0,sigma_0];

% Discretizacion %
N = 500; 
y = linspace(0,L,N);

% Solucion teorica %
delta_teo = (P/(E*A_0))*L - (P/(E*A_0))*y;
def_teo = -sigma_0/E*ones(1,N-1);
sigma_teo = def_teo*E;

% Calculos %
A = Area(y,prop,prob);
r = Radius(y,prop,prob);

H = [1 -1;-1 1];
le = y(2)-y(1);

% Matriz de conectividad %
desordenar = false;
ordenar = false;

Nodos = 1:N;
if desordenar
    Nodos = Nodos(randperm(N));
end

M = zeros(N-1,2);
for i=1:N-1
    M(i,1) = Nodos(i);
    M(i,2) = Nodos(i+1);
end


% Matriz de rigidez K %
K = zeros(N,N);
for k=1:N-1
    ke = (1/le^2)*E*H*sigma_0/gamma*(Area(y(k+1),prop,prob)-Area(y(k),prop,prob));
    for i=1:2
        for j=1:2
            I = M(k,i);
            J = M(k,j);
            K(I,J) = K(I,J) + ke(i,j);
        end
    end
end

% Vector de fuerzas F %
F = zeros(N,1);
for k=1:N-1
    dA_e = sigma_0/gamma*(Area(y(k+1),prop,prob)-Area(y(k),prop,prob));
    dyA_e = (sigma_0/gamma)^2*((gamma/sigma_0*y(k+1)-1)*Area(y(k+1),prop,prob)-(gamma/sigma_0*y(k)-1)*Area(y(k),prop,prob));
    fe_1 = -gamma/le*(y(k)*dA_e - dyA_e);
    fe_2 = gamma/le*(y(k+1)*dA_e - dyA_e);
    I = M(k,1);
    J = M(k,2);
    F(I) = F(I) + fe_1;
    F(J) = F(J) + fe_2;
end
F(M(1,1)) = F(M(1,1)) + P;

% ordenar según algoritmo symrcm %
if ordenar
    r = symrcm(K);
    K = K(r,r);
    F = F(r);
end

% Eliminar ultimo nodo %
K(M(end,2),:) = [];
K(:,M(end,2)) = [];
F(M(end,2)) = [];


% Resolucion sistema lineal %
Q = zeros(1,N-1);
Q = linsolve(K,F);

% Agregar nodo eliminado %
Q = [Q(1:M(end,2)-1); 0; Q(M(end,2):end)];

% Calculo de deformaciones y esfuerzos %
epsilon = zeros(1,N-1);
B = 1/le*[-1 1];
for k=1:N-1
    q = [Q(M(k,1)) ; Q(M(k,2))];
    epsilon(k) = B*q;
end
sigma = E*epsilon;

% Calculo de errores %
eQ = calc_L2_error(Q,delta_teo);
eS = calc_L2_error(sigma,sigma_teo);

save_files(sigma,Q,N);

disp("Resultados guardados en carpeta results")

function A = Area(y,prop,prob)
    r = num2cell(prop);
    [~,gamma,~] = deal(r{:});
    r = num2cell(prob);
    [~,~,~,A_0,sigma_0] = deal(r{:});
    A = A_0*exp(gamma/sigma_0*y);
end

function R = Radius(y,prop,prob)
    r = num2cell(prop);
    [~,gamma,~] = deal(r{:});
    r = num2cell(prob);
    [~,~,r_0,~,sigma_0] = deal(r{:});
    R = r_0*exp(gamma/(2*sigma_0)*y);
end

function [error] = calc_L2_error(u_aprox,u_teo)
    dif_u = (u_aprox - u_teo);
    error = sqrt(sum(dif_u.^2)/sum(u_teo.^2));
end

function save_files(sigma,Q,N)
    directory = 'results';
    filename = fullfile(directory, sprintf('Sigma_%d.dat', N));
    save(filename,"sigma", "-ascii")
    
    filename = fullfile(directory, sprintf('Q_%d.dat', N));
    save(filename,"Q", "-ascii")
end