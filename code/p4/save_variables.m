format compact

% Cargar y guardar malla de PDETool %
save('vars/vars.mat','p','t','e')
load('vars/vars.mat','p','t','e')

% Vertices %
cords = p';
save vars/coordinates.dat cords -ascii

% Elementos %
el3 = t(1:3,:)';
save vars/elements3.dat el3 -ascii

% Condiciones borde dirichlet %
e2 = e';
de = e2(e2(:, 5) == 1 | e2(:, 5) == 2, 1:2);
dirich = de;
save vars/dirichlet.dat dirich -ascii

% Condiciones borde neumann %
ne = e2(e2(:, 5) == 3 | e2(:, 5) == 4, 1:2);
neum = ne;
save vars/neumann.dat neum -ascii