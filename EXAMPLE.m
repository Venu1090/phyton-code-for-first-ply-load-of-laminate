clc;
clear all;

% Properties of lamina
n = input('Number of plies: ');
t = input('Thickness of each ply (in m): ');
D = input('Enter the properties of lamina in GPa ( E1,E2,V12,G12 ):   ', 's');
D = strrep(D, ',', ' ');
P = sscanf(D, '%f');

% strength parameters
S_1 = input('Enter the strength properties of lamina in MPa in MATERIAL DIRECTION - 1 ( (σ_T),(σ_C) ): ', 's');
S_1 = strrep(S_1, ',', ' ');
sigma_1 = sscanf(S_1, '%f');
S_2 = input('Enter the strength properties of lamina in MPa in MATERIAL DIRECTION - 2 ( (σ_T),(σ_C) ): ', 's');
S_2 = strrep(S_2, ',', ' ');
sigma_2 = sscanf(S_2, '%f');
tau = input('Shear strength: ');
sigma_1 = sigma_1 * 1.0e6; % For direction - 1
sigma_2 =  sigma_2 * 1.0e6; % For direction - 2
tau = tau*1.0e6;

% Angle of lamina
ang = input('Enter the angles of lamina in degrees (theta_1, theta_2, theta_3... ): ','s');
ang = strrep(ang, ',', ' ');
A_dash = sscanf(ang, '%f');

disp('===========================================================================================================')

x = input('Whether temerature is also changing (yes = 1 or no = 2): ');
if x == 1
    T = input('Temperature change (in C): ');
    alpha = input('Enter alpha values for lamina in 10 ^ -6 * m/m/C ( a_1,a_2,a_12):   ', 's');
    alpha = strrep(alpha, ',', ' ');
    Alpha = sscanf(alpha, '%f');
end

% converting angles to Radian
Ang = (pi/180) * A_dash;

% Finding Q matrix
Q = Q(P);

% Finding the Q_bar matix
Q_bar = zeros(3,3*n);
for i = 1:n
    theta = Ang(i);
    Q_bar(:,3*i-2:3*i) = Qbar(P,theta);
end

% Modifying 'P' to be in Pa
for i = [1,2,4]
    P(i) = P(i) * 1.0e9;
end

% finding the z_k , z_k-1
Z = zeros(n+1,1);
Z(1) = - (n * t) / 2;
for i = 2:n+1
    Z(i) = Z(i-1) + t;
end

% Finding the 'A','B' and 'D' matrix
A = 0;
B = 0;
D = 0;
for i = 1:n
    A = A + ( Z(i+1) - Z(i) ) *  Q_bar(:,3*i-2:3*i);
    B = B + 0.5 * ( (Z(i+1))^2 - (Z(i))^2 ) * Q_bar(:,3*i-2:3*i);
    D = D + 0.333 * ( (Z(i+1))^3 - (Z(i))^3 ) * Q_bar(:,3*i-2:3*i);
end

% Finding the ABBD matrix
ABBD = [A,B;B,D];

% Finding the strain in global co-ordinates
N = [100000;0;0;0;0;0];
e__xy = inv(ABBD) *  N;
e_xy = zeros(3,n);
for i = 1:n
    e_xy(:,i) = e__xy(1:3) + Z(i) * e__xy(4:6) ;  % strains for each ply in Global axis
end

% Finding the strains in material's axis
e_12 = zeros(3,n);
R = [1,0,0;0,1,0;0,0,2];
for i = 1:n
    e_12(:,i) = (R) * transformation(Ang(i)) * inv(R) * e_xy(:,i);
end

% Finding the stresses in each ply in material axis direction
str_dash = zeros(3,n);
for i = 1:n
    str_dash(:,i) = Q * e_12(:,i);
end

str = str_dash';

% Applying the maximum stress theory
SR = zeros(n,3);
for i = 1:n

    for j = 1:3
        if j == 1
            if str(i,j) > 0
                SR(i,j) = str(i,j) / sigma_1(1);

            else
                SR(i,j) = str(i,j) / sigma_1(2);
            end
        end
        if j == 2
            if str(i,j) > 0
                SR(i,j) = str(i,j) / sigma_2(1);
            else
                SR(i,j) = str(i,j) / sigma_2(2);
            end
        end
        if j == 3
            SR(i,j) = str(i,j) / tau;
        end
    end
end

% Finding the First ply failure mode
[~, idx] = max(abs(SR(:)));
[layer, dir] = ind2sub(size(SR), idx);
N_x = 100000 / SR(layer , dir);


disp('=======================================================================================================================')
% Now if laminate is subjected to temp change

if x == 1

    % Modifying Alpha
    Alpha = Alpha * 1.0e-6;

    % Finding the Alpha_xy for each layer
    Alpha__xy = zeros(3,n);
    for i = 1:n
        Alpha__xy(:,i) = inv(transformation(Ang(i))) * Alpha;
    end
    Alpha_xy = [Alpha__xy(1,:);Alpha__xy(2,:);2 * Alpha__xy(3,:)];

    % Finding the equivalent thermal load and moment due to T
    n_t = zeros(3,1);
    m_t = zeros(3,1);
    for i = 1:n
        n_t = n_t + T * ( Z(i+1) - Z(i) ) *  Q_bar(:,3*i-2:3*i) * Alpha_xy(:,i);
    end

    % Finding the strains in all the layers
    N_t = [n_t;m_t];
    e__xyt = inv(ABBD) * N_t; % strains due to temperature at the mid surface
    e_xyt = zeros(3,n);
    for i = 1:n
        e_xyt(:,i) = e__xyt(1:3) + Z(i) * e__xyt(4:6) ;
    end

    % Finding free thermal strains in all the layers
    e_xyf = zeros(3,n);
    for i = 1:n
        e_xyf(:,i) = T * Alpha_xy(:,i);
    end

    % Finding the residual strains in all the layers
    e_xyr = e_xyt - e_xyf;

    % Residual stresses in plies in xy co odrdinate
    str_r = zeros(3,n);
    for i = 1:n
        str_r(:,i) = Q_bar(:,3*i-2:3*i) * e_xyr(:,i);
    end

    % Residual stresses in plies in material directions
    str_r12 = zeros(3,n);
    for i = 1:n
        str_r12(:,i) = transformation(Ang(i)) * str_r(:,i);
    end
    %         disp('=============================================================================================')
    %         fprintf('For Temp. change of %d degree celcius, Residual thermal stress in %d degree layer: %f Pa\n',T,A_dash(layer),str_r12(dir,layer));

    % Finding the total stress in concerned ply
    Sigma = str(layer,dir) + str_r12(dir,layer);

    % Calculation for FPF load for thermal load
    if (1 <= idx) && (idx <= n) && str(layer,dir) > 0
        sigma_u = sigma_1(1);
    elseif (1 <= idx) && (idx <= n) && str(layer,dir) < 0
        sigma_u = sigma_1(2);
    elseif (n+1 <= idx) && (idx <= 2*n) && str(layer,dir) > 0
        sigma_u = sigma_2(1);
    elseif (n+1 <= idx) && (idx <= 2*n) && str(layer,dir) < 0
        sigma_u = sigma_2(2);
    elseif (2*n + 1 <= idx) && (idx <= 3 * n)
        sigma_u = tau;
    end
    sigma2 = sigma_u - str_r12(dir,layer);
    % FPF load at which sigma2 will attain its value
    N_xt = 100000 / ( ( str(layer,dir) ) / ( sigma2 ) );
    %         disp('=============================================================================================')
    fprintf('%d degree layer(s) will fail with temperature change of %d at load of %f N/m\n',A_dash(layer),T,N_xt);
elseif x == 2
    disp('No temperature change occured');
    fprintf('%d degree layer(s) will fail with load of %f N/m\n',A_dash(layer),N_x);
else
    disp('Enter correct data as directed');
end
