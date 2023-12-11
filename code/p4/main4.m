format compact
clc

% Propiedades %
E = 1e4; nu = 0.25;
mu = E/ (2* (1+nu)) ; lambda = E*nu/((1+nu) * (1-2*nu)) ;

% Cargar malla, se debe correr codigo save_variables.m %
load vars/coordinates.dat;
eval("load vars/elements3.dat;", "elements3 = [];");
eval("load vars/elements4.dat;", "elements4 = [];");
eval("load vars/neumann.dat;", "neumann = [];");
load vars/dirichlet.dat;

if ~exist('plots', 'dir')
    mkdir('plots');
end

% Matriz de rigidez %
K = sparse(2*size(coordinates, 1),2*size(coordinates, 1));
for j = 1:size(elements3, 1)
    I = 2*elements3(j, [1 1 2 2 3 3]) -[1 0 1 0 1 0];
    K(I,I) = K(I,I) + stima3(coordinates(elements3(j,:),:), lambda, mu);
end

F = zeros(2*size(coordinates, 1), 1);

% Fuerzas de cuerpo %
for j = I:size (elements3, 1)
    I = 2*elements3(j,[1 1 2 2 3 3]) - [1 0 1 0 1 0];
    fs = f(sum(coordinates(elements3(j,:),:))/3)';
    F(I) = F(I) + det([1, 1, 1 ; coordinates(elements3(j,:),:)'])*[fs;fs;fs] /6;
end

% Condiciones de neumann %
if ~isempty (neumann)
    n = (coordinates (neumann (: ,2) , :) -coordinates (neumann (:, 1),:))* [0, -1;1,0];
    for j = 1: size (neumann, 1)
        I = 2*neumann (j, [1,1,2,2]) - [1,0,1,0];
        gm = g (sum(coordinates (neumann (j,:) , :)) /2, n(j, :) /norm (n(j,:)))';
        F(I) = F(I) + norm(n(j, :)) * [gm; gm] /2;
    end
end

% Condiciones de dirichlet %
DirichletNodes = unique(dirichlet);
[W,M] = u_d(coordinates(DirichletNodes, :));
B = sparse(size(W, 1) ,2*size(coordinates, 1));
for k = 0:1
    for l = 0:1
        B(1+l:2:size(M, 1),2*DirichletNodes-1+k) = diag(M(1+l:2:size(M, 1), 1+k));
    end
end

mask = find(sum(abs(B)'));
K = [K, B(mask,:)'; B(mask,:), sparse(length(mask), length(mask))];
F = [F; W(mask, :)];

% Ordenar segun symrcm %
ordenar = true;
if ordenar
    r = symrcm(K);
    K = K(r,r);
    F = F(r);
end

% Resolver sistema lineal % 
x = K \ F;

% Devolver numeracion si se reordeno con symrcm %
if ordenar
    [~,r_2] = sort(r); 
    x = x(r_2);
end

% Obtener desplazamientos %
u = x(1:2*size(coordinates,1));

% Calcular esfuerzos y errores %
[AvE, Eps3, Eps4, AvS, Sigma3, Sigma4] = avmatrix(coordinates,elements3,elements4,u,lambda,mu);
[Sigma_principal,von_mises,AvC] = extra_stresses(AvS,lambda,mu);
[L2_error_u, L2_error_von_mises, mesh_error] = calculate_errors(coordinates, elements3, elements4, AvE, Eps3, Eps4, AvS, Sigma3, Sigma4, u, lambda, mu);

% Graficar solucion %
%plot_solutions(elements3, elements4, coordinates, AvS, Sigma_principal, von_mises, AvC, u, lambda, mu);

% Imprimir resultados %
disp("Numero de nodos")
disp(size(coordinates,1))
disp("Von mises max")
disp( max(von_mises))
disp("un max")
disp(max(sqrt(u(1:2:end,1).^2 + u(2:2:end,1).^2))*10^3)
L2_error_u
L2_error_von_mises
mesh_error

disp("Gráficos guardados en plots")

% Condicion de borde de dirichlet %
function [W,M] = u_d(x)
    M = zeros(2*size(x,1),2);
    W = zeros(2*size(x,1),1);
    M(1:2:end,1) = 1;
    temp = find(x(:,1)<1e-9);
    M(2*temp,2) = 1;
    value = u_value(x);
    W(1:2:end,1)= value(:,1);
    W(2*temp,1) = value(temp,2);
end

% Fuerzas de cuerpo %
function volforce = f(x)
    volforce = zeros(size(x,1),2);
end

% Fuerzas de superficie %
function sforce = g(x,n)
    sforce = zeros(size(x,1),2);
end

% Matriz rigidez local %
function stima3=stima3(vertices, lambda, mu)
    PhiGrad = [1, 1, 1; vertices']\ [zeros(1,2) ; eye(2)];
    R = zeros (3,6);
    R([1,3], [1,3,5]) = PhiGrad';
    R([3,2], [2,4,6]) = PhiGrad';
    C = mu*[2,0,0;0,2, 0; 0, 0, 1] + lambda* [1, 1,0;1, 1, 0;0,0,0];
    stima3 = det([1,1,1;vertices'])/2*R'*C*R;
end


% Solucion teorica para desplazamientos %
function value = u_value(x)
    [phi,r] = cart2pol(x(:,1),x(:,2));
    E = 1e4; nu = 0.25;
    mu = E/ (2* (1+nu)) ; lambda = E*nu/((1+nu) * (1-2*nu)) ;
    a = 5; b = 10; u0 = -0.01;
    ab = a^2 + b^2;
    PN = -u0*E/(pi*ab);

    K = PN/E*( 0.5*(1-3*nu)*a^2 - b^2*(1+nu)/2 - ab*(1-nu)*log(a) ); 
    ur = PN/E*( (0.5*(1-3*nu)*r.^2 - a^2*b^2*(1+nu)/(2*r.^2) - ab*(1-nu)*log(r)).*sin(phi) + ab*(2*phi-pi).*cos(phi)) - K*sin(phi) ;
    ut = -PN/E*( (0.5*(5+nu)*r.^2 - a^2*b^2*(1+nu)/(2*r.^2) + ab*(1-nu)*log(r) + (1+nu) ).*cos(phi) + ab*(2*phi-pi).*sin(phi)) - K*cos(phi) ;
    value = [ur.*cos(phi)-ut.*sin(phi), ur.*sin(phi)+ut.*cos(phi)];
end


% Solucion teorica para esfuerzos %
function sig_cartesian = analytic_Sig(x)
    [phi,r] = cart2pol(x(:,1),x(:,2));
    E = 1e4; nu = 0.25;
    mu = E/ (2* (1+nu)) ; lambda = E*nu/((1+nu) * (1-2*nu)) ;
    a = 5; b = 10; u0 = -0.01;
    ab = a^2 + b^2;
    PN = -u0*E/(pi*ab);
    
    srr = PN*(r + a^2*b^2./r.^3 - ab./r).*sin(phi);
    stt = PN*(3*r - a^2*b^2./r.^3 - ab./r).*sin(phi);
    str = -PN*(r + a^2*b^2./r.^3 - ab./r).*cos(phi);
    Sig_polar = [srr str str stt];

    cos_phi = cos(phi);
    sin_phi = sin(phi);
    Q = [cos_phi -sin_phi sin_phi cos_phi];
       
    SQ = zeros(size(r,1),4);
    SQ(:,1) = Sig_polar(:,1).*Q(:,1)+Sig_polar(:,2).*Q(:,3);
    SQ(:,2) = Sig_polar(:,1).*Q(:,2)+Sig_polar(:,2).*Q(:,4);
    SQ(:,3) = Sig_polar(:,3).*Q(:,1)+Sig_polar(:,4).*Q(:,3);
    SQ(:,4) = Sig_polar(:,3).*Q(:,2)+Sig_polar(:,4).*Q(:,4);
    SC = zeros(size(r,1),4);
    SC(:,1) = Q(:,1).*SQ(:,1)+Q(:,3).*SQ(:,3);
    SC(:,2) = Q(:,1).*SQ(:,2)+Q(:,3).*SQ(:,4);
    SC(:,3) = Q(:,2).*SQ(:,1)+Q(:,4).*SQ(:,3);
    SC(:,4) = Q(:,2).*SQ(:,2)+Q(:,4).*SQ(:,4); 

    sig_cartesian = SC;
end


% Calculo de esfuerzos y deformaciones % 
function [AvE, Eps3, Eps4, AvS, Sigma3, Sigma4] = avmatrix(coordinates, elements3, elements4, u, lambda, mu)
    Eps3 = zeros(size(elements3, 1) ,4);
    Sigma3 = zeros(size(elements3, 1),4);
    Eps4 = zeros(size (elements4, 1) ,4);
    Sigma4 = zeros(size(elements4, 1), 4);
    AreaOmega = zeros(size(coordinates, 1), 1);
    AvS = zeros(size(coordinates, 1),4);
    AvE = zeros(size(coordinates, 1),4);
    for j = 1:size(elements3, 1)
        area3 = det ([1,1,1;coordinates(elements3(j,:),:)'])/2;
        AreaOmega(elements3(j,:)) = AreaOmega (elements3 (j,:)) +area3;
        PhiGrad = [1, 1,1; coordinates(elements3(j,:) , :)'] \ [zeros(1,2) ; eye(2)];
        U_Grad = u([1;1]*2*elements3(j,:)-[1;0]*[1, 1, 1])*PhiGrad;
        Eps3(j,:) = reshape((U_Grad+U_Grad') /2, 1,4);
        Sigma3 (j,:) = reshape(lambda*trace(U_Grad)*eye(2) +2*mu* (U_Grad+U_Grad') /2,1,4) ;
        AvE(elements3(j,:),:) = AvE(elements3(j,:),:) +area3* [1;1;1]*Eps3(j,:);
        AvS(elements3(j,:),:) = AvS(elements3(j,:) , :) +area3* [1;1;1]*Sigma3(j,:);
    end
    AvE = AvE./ (AreaOmega*[1, 1, 1, 1]);
    AvS = AvS./ (AreaOmega*[1,1,1,1]);
end

% Calculo de esfuerzos principales y von mises %
function [Sigma_principal,von_mises,AvC] = extra_stresses(AvS,lambda,mu)
    n = size(AvS);
    sum_Sig = AvS(:,1) + AvS(:,4);
    dif_Sig = AvS(:,1) - AvS(:,4);
    Sigma_principal = zeros(size(AvS, 1),3);
    Sigma_principal(:,1) = 1/2*sum_Sig + sqrt(1/4*dif_Sig.^2 + AvS(:,2).^2);
    Sigma_principal(:,2) = 1/2*sum_Sig - sqrt(1/4*dif_Sig.^2 + AvS(:,2).^2);
    Sigma_principal(:,3) = lambda/(2*(mu+lambda))*sum_Sig;
    von_mises = sqrt( (Sigma_principal(:,1)-Sigma_principal(:,2)).^2 + (Sigma_principal(:,2)-Sigma_principal(:,3)).^2 + (Sigma_principal(:,1)-Sigma_principal(:,3)).^2)/sqrt(2); 
    AvC=(mu/ (24* (mu+lambda)^2)+1/(8*mu))*(AvS(:, 1) + AvS(:,4)).^2+1/(2*mu)*(AvS(:,2).^2-AvS(:,1).*AvS(:,4)) ;
end


% Calculo de errores % 
function [L2_error_u, L2_error_von_mises, mesh_error] = calculate_errors(coordinates, elements3, elements4, AvE, Eps3, Eps4, AvS, Sigma3, Sigma4, u, lambda, mu)
    ux = u(1:2:end,1);
    uy = u(2:2:end,1);
    u_n = sqrt(ux.^2 + uy.^2)*10^3;
    u_teo = u_value(coordinates);
    ux_teo = u_teo(:,1);
    uy_teo = u_teo(:,2);
    u_n_teo = sqrt(ux_teo.^2 + uy_teo.^2)*10^3;
    L2_error_u = calc_L2_error(u_n,u_n_teo);

    sig_teo = analytic_Sig(coordinates);
    [~,von_mises_teo,AvC_teo] = extra_stresses(sig_teo,lambda,mu);
    [~,von_mises,AvC] = extra_stresses(AvS,lambda,mu);
    L2_error_von_mises = calc_L2_error(von_mises,von_mises_teo);

    mesh_error = aposteriori(coordinates, elements3, elements4, AvE, Eps3, Eps4, AvS, Sigma3, Sigma4, u, lambda, mu);
end

function [error] = calc_L2_error(u_aprox,u_teo)
    dif_u = (u_aprox - u_teo);
    error = sqrt(sum(dif_u.^2)/sum(u_teo.^2));
end

% Calculo del error de la malla %
function estimate=aposteriori(coordinates, elements3, elements4,AvE, Eps3, Eps4, AvS, Sigma3, Sigma4, u, lambda, mu);
    eta3 = zeros(size(elements3,1),1) ;
    eta4 = zeros(size(elements4,1),1);
    e3 = zeros(size(elements3,1),1);
    e4 = zeros(size (elements4,1), 1);
    for j=1:size(elements3,1)
        Area3=det([1,1,1;coordinates(elements3(j,:),:)'])/2;
        for i=1:4
            eta3(j) = eta3(j) + Area3 * (Sigma3(j,i)*Eps3(j,i) + AvS(elements3(j,:),i)'* [2,1,1;1,2,1;1,1,2]*AvE(elements3(j,:),i)/12 - AvS(elements3(j,:),i)'* [1;1;1]*Eps3(j,i)/3 - AvE(elements3 (j,:) ,i)'* [1;1;1]*Sigma3 (j,i) /3) ;
            e3(j) = e3(j) + Area3 * AvS(elements3(j,:),i)'*[2,1,1;1,2,1;1,1,2]*AvE(elements3(j,:),i)/12;
        end
    end
    estimate = sqrt (sum (eta3))/ sqrt (sum(e3));
end


% Funcion para graficar solido deformado %
function show(elements3, elements4, coordinates, X, u, lambda, mu, titulo)
    factor=60;
    colormap("jet")
    trisurf(elements3, factor*u(1:2:size(u, 1))+coordinates(:, 1),factor*u(2:2:size(u, 1))+coordinates(: ,2),zeros(size(coordinates, 1), 1), X, 'facecolor', 'interp');
    hold on
    view (0,90)
    hold off
    colorbar("vert")
    xlabel('x')
    ylabel('y')
    title(titulo)
end

% Creacion de graficos %
function plot_solutions(elements3, elements4, coordinates, AvS, sigma_principal, von_mises, AvC, u, lambda, mu)
    n = size(coordinates,1);    
    ux = u(1:2:end,1);
    uy = u(2:2:end,1);
    sx = AvS(:,1);
    sy = AvS(:,4);
    txy = AvS(:,2);

    var_plot = von_mises;
    titulo = 'Esfuerzo de Von Mises';
    show(elements3, elements4, coordinates, var_plot, u, lambda, mu,titulo);
    filename = fullfile('plots', sprintf('von_%d.png', n));
    saveas(gcf,filename)

    var_plot = ux;
    titulo = 'Desplazamiento ux';
    show(elements3, elements4, coordinates, var_plot, u, lambda, mu,titulo);
    filename = fullfile('plots', sprintf('ux_%d.png', n));
    saveas(gcf,filename)
    var_plot = uy;
    titulo = 'Desplazamiento uy';
    show(elements3, elements4, coordinates, var_plot, u, lambda, mu,titulo);
    filename = fullfile('plots', sprintf('uy_%d.png', n));
    saveas(gcf,filename)

    var_plot = sx;
    titulo = 'Esfuerzo \sigma_{xx}';
    show(elements3, elements4, coordinates, var_plot, u, lambda, mu,titulo);
    filename = fullfile('plots', sprintf('sx_%d.png', n));
    saveas(gcf,filename)
    var_plot = sy;
    titulo = 'Esfuerzo \sigma_{yy}';
    show(elements3, elements4, coordinates, var_plot, u, lambda, mu,titulo);
    filename = fullfile('plots', sprintf('sy_%d.png', n));
    saveas(gcf,filename)

    var_plot = txy;
    titulo = 'Esfuerzo \tau_{xy}';
    show(elements3, elements4, coordinates, var_plot, u, lambda, mu,titulo);
    filename = fullfile('plots', sprintf('txy_%d.png', n));
    saveas(gcf,filename)  
end
