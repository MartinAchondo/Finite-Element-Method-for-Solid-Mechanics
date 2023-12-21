format compact
clc

% Dimension geometria %
a = 1;
b = 2;

% Propiedad material %
G_mu = 5;

% Termino fuente %
f0 = 2*G_mu;

% Creacion malla %
Hmax = 0.01;
[coordinates,elements,dirichlet] = create_mesh(a,b,Hmax);
[centroids] = get_centroids(coordinates,elements);

% Matriz de rigidez %
K = sparse(size(coordinates, 1), size(coordinates, 1));
for j = 1:size(elements, 1)
    I = elements(j, [1 2 3]);
    K(I, I) = K(I, I) + local_K(coordinates(elements(j,:),:));
end

% Fuerzas de cuerpo %
F = zeros(size(coordinates, 1), 1);
for j = 1:size (elements, 1)
    I = elements(j, [1 2 3]);
    F(I) = F(I) + local_F(coordinates(elements(j,:),:),centroids(j),f0);
end

% Condiciones de dirichlet %
DirichletNodes = unique(dirichlet);
d_values = u_d(coordinates(DirichletNodes,:));

% Sistema reducido %
excludeIdx = true(size(K,1),1);
excludeIdx(DirichletNodes) = false;
K_modified = K(excludeIdx,excludeIdx);
F_modified = F(excludeIdx);

% Resolucion %
x = K_modified \ F_modified;

% Construccion solucion %
u = zeros(size(F));
u(excludeIdx) = x;
u(DirichletNodes) = d_values(:);

% Calculo del esfuerzo cortante %
tau = zeros(size(coordinates,1),2);
for j = 1:size(elements, 1)
    I = elements(j, [1 2 3]);
    tau(j,:) = get_tau(coordinates(elements(j,:),:),u(I));
end
tau_yz = -tau(:,1);
tau_xz = tau(:,2);
tau = sqrt(tau_yz.^2+tau_xz.^2);

% Calculo del torque %
Torque = get_Torque(coordinates,elements,u);


% Postprocessing %

disp(" ")

disp("Valor maximo u")
disp(max(u))

disp("Valor máximo tau")
disp(max(tau))

disp("Torque")
disp(Torque)

disp("Error L2 u")
u_an = analytic_u(coordinates,a,b,G_mu);
disp(strcat("  ",num2str(calc_L2_error(u,u_an),'%.3e')))

disp("Error L2 tau")
tau_an = analytic_tau(centroids,a,b,G_mu);
disp(strcat("  ",num2str(calc_L2_error(tau,tau_an),'%.3e')))

disp("Error L2 torque")
Torque_an = analytic_torque(a,b,G_mu);
disp(strcat("  ",num2str(calc_L2_error(Torque,Torque_an),'%.3e')))


plot_solutions(coordinates,elements,u,tau,b)


% Funciones %

% Condiciones de Dirichlet %
function [values] = u_d(x)
    values = zeros(size(x,1),1);
end

% Fuerzas de cuerpo %
function volforce = f(x,f0)
    volforce = ones(size(x,1),3)*f0;
end

% Matriz de rigidez local %
function local_K = local_K(vertices)
    x = vertices(:,1);
    y = vertices(:,2);
    A = det([1,1,1;vertices'])/2;
    det_J = A*2;
    B = 1/det_J*[y(2)-y(3), y(3)-y(1), y(1)-y(2);
                x(3)-x(2), x(1)-x(3), x(2)-x(1)];
    k = A*B'*B;
    local_K = k;
end

% Fuerza en elementos %
function local_F = local_F(vertices,centroid,f0)
    A = det([1,1,1;vertices'])/2;
    fs = f(centroid,f0)';
    local_F = fs*A/3;
end

% Calcular esfuerzo en elemento %
function [tau] = get_tau(vertices,u)
    A = det([1,1,1;vertices'])/2;
    x = vertices(:,1);
    y = vertices(:,2);
    det_J = A*2;
    B = 1/det_J*[y(2)-y(3), y(3)-y(1), y(1)-y(2);
                x(3)-x(2), x(1)-x(3), x(2)-x(1)];
    tau = B*u;
end

% Calcular el torque %
function [Torque] = get_Torque(coordinates,elements,u)
    Torque = 0;
    for j = 1:size(elements, 1)
        I = elements(j, [1 2 3]);
        vertices = coordinates(elements(j,:),:);
        A = det([1,1,1;vertices'])/2;
        Torque = Torque + A*mean(u(I));
    end
    Torque = Torque*4*2;
end

% Soluciones analiticas % 
function [u] = analytic_u(X,a,b,G_mu)
   x = X(:,1);
   y = X(:,2);
   u = G_mu*a^2*b^2/(a^2+b^2)*(1-x.^2/a^2-y.^2/b^2);
end

function [tau] = analytic_tau(X,a,b,G_mu)
    x = X(:,1);
    y = X(:,2);
    T = analytic_torque(a,b,G_mu);
    tau = 2*T/(pi*a*b)*sqrt(x.^2/a^4+y.^2/b^4);
end

function [torque] = analytic_torque(a,b,G_mu)
    torque = pi*a^3*b^3*G_mu/(a^2+b^2);
end

% Calculo de error L2 %
function [error] = calc_L2_error(u_aprox,u_teo)
    dif_u = (u_aprox - u_teo);
    error = sqrt(sum(dif_u.^2)/sum(u_teo.^2));
end

% Funciones para crear la malla %
function [centroids] = get_centroids(coordinates,elements)
    centroids = zeros(size(elements,1),2);
    for j = 1:size(elements, 1)
        centroids(j,:) = mean(coordinates(elements(j,:),:),1);
    end
end

function [coordinates,elements3,dirichlet] = create_mesh(a,b,Hmax) 
    R1 = [3,4,-a,0,0,-a,-b,-b,b,b]';
    R2 = [3,4,-a,a,a,-a,-b,-b,0,0]';
    C1 = [4,0,0,a,b]'; 
    C1 = [C1;zeros(length(R1)-length(C1),1)];
    gm = [R1,R2,C1]; 
    ns = char('R1','R2','C1'); ns=ns';
    sf = 'C1-R1-R2';
    d1 = decsg(gm,sf,ns); 
    model = createpde;
    geometryFromEdges(model,d1); 
    Hgrad = 1.5;
    generateMesh(model,'Hmax',Hmax,'Hgrad',Hgrad,'Jiggle','on', 'GeometricOrder','linear')
    [p,e,t] = meshToPet(model.Mesh);
    coordinates = p';
    elements3 = t(1:3,:)';
    e2 = e';
    dirichlet = e2(e2(:,5)==3, 1:2);
    %pdegplot(d1,"EdgeLabels","on","FaceLabels","on")%
end

% Funciones para graficar la solucion %
function plot_solutions(coordinates,elements,u,tau,b)

    if exist('plots', 'dir') ~= 7
        mkdir('plots');
    end

    title = "Función esfuerzo de Prandtl u";
    show(elements,coordinates,u,b,title)
    filename = fullfile('plots', sprintf('u_n%d_b%.1f.png', size(coordinates, 1), b));
    saveas(gcf,filename)  

    title = "Esfuerzo cortante \tau";
    show2(elements,coordinates,tau,b,title)
    filename = fullfile('plots', sprintf('tau_n%d_b%.1f.png', size(coordinates, 1), b));
    saveas(gcf,filename) 
end

function show(elements3, coordinates, X,b,plot_title)
    figure    
    colormap('jet')
    trisurf(elements3, coordinates(:, 1),coordinates(: ,2),zeros(size(coordinates, 1), 1), X, 'facecolor', 'interp');
    view (0,90)
    colorbar('vert')
    xlabel('x')
    ylabel('y')
    xlim([0,b])
    title(plot_title)
end

function show2(elements, coordinates, X,b,plot_title)
    figure
    colormap('jet')
    patch('Faces',elements,'Vertices',coordinates,'FaceVertexCData',X,'FaceColor','flat');
    view(0,90)
    colorbar('vert')
    xlabel('x')
    ylabel('y')
    xlim([0,b])
    title(plot_title)
end