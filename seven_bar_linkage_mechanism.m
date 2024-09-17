clear all
clc

% Define the matrix representing the links connected to other links
A = [0 1 0 1 0 0 1;
     1 0 1 0 0 0 0;
     0 1 0 1 1 0 0;
     1 0 1 0 0 0 0;
     0 0 1 0 0 1 0;
     0 0 0 0 1 0 1;
     1 0 0 0 0 1 0];

% Calculate different parameters
totalSum = sum(A(:));
j1 = totalSum / 2;
order = size(A, 1);
dof = (3 * (order - 1)) - 2 * j1;

% Link lengths
AO2 =2;
AB = 3;
AC = 2;
BC = 3;
PD = 2;
DO7 = 1;
BO4 = 3;

% Fixed angles of triangular link 3
alpha1 = acosd((AB^2 + AC^2 - BC^2) / (2 * AB * AC));
alpha2 = acosd((AB^2 + BC^2 - AC^2) / (2 * AB * BC));
alpha3 = acosd((AC^2 + BC^2 - AB^2) / (2 * BC * AC));

% Create video writer object
writerobj = VideoWriter('simulation_of_mechanism', 'MPEG-4');
open(writerobj);

% Initial positions of fixed points
O2x = 0; O2y = 0;
O4x = 3; O4y = 0;
O7x = 7; O7y = 0;

% Lengths between fixed points
O2_O4 = sqrt((O4x - O2x)^2 + (O4y - O2y)^2);
O4_O7 = sqrt((O7x - O4x)^2 + (O7y - O4y)^2);

% Loop for theta2 and theta7
theta2 = 0;
theta7 = 180;

fig1 = figure(1);

while (theta2 <= 360 && theta7 >=-180)
    % Calculations for the positions and angles
    AO4 = sqrt(AO2^2 + O2_O4^2 - 2 * AO2 * O2_O4 * cosd(theta2));
    beta = acosd((O2_O4^2 + AO4^2 - AO2^2) / (2 * O2_O4 * AO4));
    psi = acosd((AB^2 + AO4^2 - BO4^2) / (2 * AB * AO4));
    lamda = acosd((BO4^2 + AO4^2 - AB^2) / (2 * BO4 * AO4));
    theta4 = 180 - lamda - beta;
    BO7 = sqrt(BO4^2 + O4_O7^2 - 2 * BO4 * O4_O7 * cosd(theta4));
    gama = acosd((BO7^2 + O4_O7^2 - BO4^2) / (2 * BO7 * O4_O7));
    delta = 180 - theta4 - gama;
    kai = 180 - gama - theta7;
    DB = sqrt(BO7^2 + DO7^2 - 2 * BO7 * DO7 * cosd(kai));
    theta3 = psi - beta;

    % Defining joints positions
    Ax = O2x + AO2 * cosd(theta2);
    Ay = O2y + AO2 * sind(theta2);
    Bx = O4x + BO4 * cosd(theta4);
    By = O4y + BO4 * sind(theta4);
    Cx = Ax + AC * cosd(alpha1 + theta3);
    Cy = Ay + AC * sind(alpha1 + theta3);
    Dx = O7x + DO7 * cosd(theta7);
    Dy = O7y + DO7 * sind(theta7);

    theta6 = atan2d(Cy - Dy, Cx - Dx);

    Px = Dx + PD*cosd(theta6);
    Py = Dy + PD*sind(theta6);
    
    PC = sqrt((Px-Cx)^2+(Py-Cy)^2);


     % Solve for P coordinates3
    % syms x y;
    % eq1 = (x - Dx)^2 + (y - Dy)^2 == DP^2;
    % sol = solve(eq1, [x, y]);


    % Px = double(sol.x)
    % Py = double(sol.y)
    % disp(Px)

    % Ensure to take the correct solution for D
    % dx=Dx(1);
    % dy=0;
    % for i=1:length(Dy)
    %     if real(Dy(i))<0
    %         dy=Dy(i);
    %     end
    % end
    
    if length(Px) > 1
       Px = Px(1);
       Py = Py(1);
    end

    % Dx=dx;
    % Dy=dy;

    % Plot the mechanism
   


    plot([O2x Ax], [O2y Ay], [Ax Bx], [Ay By], [Ax Cx], [Ay Cy], ...
         [Cx Bx], [Cy By], [Bx O4x], [By O4y], [Cx Px], [Cy Py], ...
          [Px Dx], [Py Dy],[Dx O7x], [Dy O7y], 'LineWidth', 3);

    grid on;
    axis equal;
    axis([-5 15 -5 10]);
    drawnow;
    hold off;

    % Write frame to video
    F = getframe(fig1);
    writeVideo(writerobj, F);

    % Increment theta values
    theta2 = theta2 + 5;
    theta7 = theta7 - 5;
end

close(writerobj);
disp("Video File Written");





%% velocity analysis 


% give your  inputs 
w2=5;
w7=6;
theta2=160;
theta7=50;
k1 = [0 0 1];
k2 = [0 0 -1];
 

    AO4 = sqrt(AO2^2 + O2_O4^2 - 2 * AO2 * O2_O4 * cosd(theta2));
    beta = acosd((O2_O4^2 + AO4^2 - AO2^2) / (2 * O2_O4 * AO4));
    psi = acosd((AB^2 + AO4^2 - BO4^2) / (2 * AB * AO4));
    lamda = acosd((BO4^2 + AO4^2 - AB^2) / (2 * BO4 * AO4));
    theta4 = 180 - lamda - beta;
    BO7 = sqrt(BO4^2 + O4_O7^2 - 2 * BO4 * O4_O7 * cosd(theta4));
    gama = acosd((BO7^2 + O4_O7^2 - BO4^2) / (2 * BO7 * O4_O7));
    delta = 180 - theta4 - gama;
    kai = 180 - gama - theta7;
    DB = sqrt(BO7^2 + DO7^2 - 2 * BO7 * DO7 * cosd(kai));
    theta3 = psi - beta;

    %% Defining joints positions

    Ax = O2x + AO2 * cosd(theta2);
    Ay = O2y + AO2 * sind(theta2);
    Bx = O4x + BO4 * cosd(theta4);
    By = O4y + BO4 * sind(theta4);
    Cx = Ax + AC * cosd(alpha1 + theta3);
    Cy = Ay + AC * sind(alpha1 + theta3);
    Dx = O7x + DO7 * cosd(theta7);
    Dy = O7y + DO7 * sind(theta7);

    theta6 = atan2d(Dy-Cy, Dx - Cx);

    Px = Dx - PD*cosd(theta6);
    Py = Dy - PD*sind(theta6);
    
    PC = sqrt((Px-Cx)^2+(Py-Cy)^2);

    %% position vectors
  

rAO2=[AO2*cosd(theta2) AO2*sind(theta2) 0];

rBA=[AB*cosd(theta3) AB*sind(theta3) 0];

rBO4=[BO4*cosd(theta4) BO4*sind(theta4) 0];

rCA=[AC*cosd(theta3+alpha1) AC*sind(theta3+alpha1) 0];

rBC=[BC*cosd(180+theta3-alpha2) BC*sind(180+theta3-alpha2) 0];

rDO7=[DO7*cosd(theta7) DO7*sind(theta7) 0];

rPD=[PD*cosd(theta6) PD*sind(theta6)  0];

rCP=[Px-Cx Py-Cy  0];

rCD=[Dx-Cx Dy-Cy 0];

rDC=[Cx-Dx Cy-Dy 0];



%% velocity calculation for point A

V_A = w2*cross(k1,rAO2);
V_A_mag=norm(V_A);


%% velocity calculation for point B

syms w3 w4

V_B1 = V_A + w3*cross(k1,rBA);
V_B2 = w4*cross(k1,rBO4);

eqn_vel = V_B1 - V_B2 == 0;
S = solve(eqn_vel,[w3 w4]);

w3 = double(S.w3);
w4 = double(S.w4);



%% velocity calculation for point C

V_C=V_A + w3*cross(k1,rCA);

%% velocity calculation for point D
V_D=w7*cross(k2,rDO7);

%%  radial velocity calculation for point P
V1=[V_C(1) V_C(2) 0];
V2=rCD;
thetaVC_CD=acosd(dot(V1,V2)/(norm(V1)*norm(V2)));


V3=[V_D(1) V_D(2) 0];
V4=rCD;
thetaVD_CD=acosd(dot(V3,V4)/(norm(V3)*norm(V4)));

V_C_CD=(norm(V_C)*cosd(thetaVC_CD))*(rCD/norm(rCD));
V_D_CD=(norm(V_D)*cosd(thetaVD_CD))*(rCD/norm(rCD));

V_PfromC_rad=V_D_CD-V_C_CD;








%% normal  velocity calculation for point p  

sym w6;

V_D_norm=V_D-V_D_CD;
V_C_norm=V_C-V_C_CD;



% Calculate the difference in normal velocities
V_diff = V_D_norm - V_C_norm;

% Construct the matrix for the cross product operation
M = [
    0, -rCD(3), rCD(2);
    rCD(3), 0, -rCD(1);
    -rCD(2), rCD(1), 0
];

% Solve for w6 using the pseudo-inverse
w6 = pinv(M) * V_diff';


V_P5fromC=V_C_norm+w6*cross(k1,rCP);







%% plot all the velocity vectors 

% Define the vector names
vector_names = {'V-A', 'V-B', 'V-C', 'V-D','V-P6-P5','V-P5-C'};

% Convert symbolic expressions to double
V_A = double(V_A);
V_B= double(subs(V_B2,w4))
V_C = double(V_C);
V_D = double(V_D)
V_P6_P5=double(V_PfromC_rad);
V_P5_C=double(V_P5fromC);

% Extract x and y components for quiver plot (ignoring z-component)
U = [V_A(1), V_B(1), V_C(1), V_D(1),V_P6_P5(1),V_P5_C(1)];
V = [V_A(2), V_B(2), V_C(2), V_D(2),V_P6_P5(2),V_P5_C(2)];

% Define the origin points (all vectors start from the same origin)
origin_x = [0,0,0,0,0,0];
origin_y = [0,0,0,0,0,0];

% Plot the vectors
figure;
quiver(origin_x, origin_y, U, V, 0, 'AutoScale', 'off'); % The last 0 disables automatic scaling
hold on;
scatter(0, 0, 'filled'); % Mark the origin
xlim([-10 10]);
ylim([-10 10]);
xlabel('X');
ylabel('Y');
title('Velocity Vectors from Origin');
grid on;
axis equal; % Ensure the aspect ratio is equal
% Add labels to the vectors
for i = 1:length(U)
    text(U(i), V(i), vector_names{i}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end
hold off;









  














