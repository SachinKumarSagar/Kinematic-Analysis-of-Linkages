% this all stuff can be done using a function also


% define a matrix representing the links is connected to other links


         A=[0 1 0 1 0 0 1;
            1 0 1 0 0 0 0;
            0 1 0 1 1 0 0;
            1 0 1 0 0 0 0;
            0 0 1 0 0 1 0;
            0 0 0 0 1 0 1;
            1 0 0 0 0 1 0];

 % calculating different parameters which are imp for rest of calculation

totalSum = sum(A(:));
j1=totalSum/2;
order = size(A, 1);
dof=(3*(order-1))-2*j1;


%disp(totalSum);
%disp(j1);
%disp(order);
%disp(dof);

% link lengths can also be taken as input in above matrix itself like 3D
% matrix in which 3rd dimention is representing the length of respective
% link


AO2=1;
AB=3;
AC=1;
BC=1;
CD=5;
DP=2;
PO7=1;
BO4=4;


%w2=omega2;
%w7=omega7;

% considering link 2 as crank it can rotate 360 


% theta2 = 0:2:360;
% theta7=0:2:360;
writerobj=VideoWriter('simulation of mechanism','MPEG-4');
open(writerobj);




% fixed angles of triangular link 3

alpha1=acosd((AB^2+AC^2-BC^2)/(2*AB*AC));

alpha2=acosd((AB^2+BC^2-AC^2)/(2*AB*BC));

alpha3=acosd((AC^2+BC^2-AB^2)/(2*BC*AC));


%  dynamic links and angles have to  be inside loop

% for i = 1:length(theta2) 

fig1 = figure(1);


theta2=0;
theta7=0;

x1=0;
y1=0;
x2=3;
y2=0;
x3=7;
y3=0;

 O2_O4=sqrt((x2-x1)^2+(y2-y1)^2);

 O4_O7=sqrt((x3-x2)^2+(y3-y2)^2);



 % define a function to find the unknown variables x nad y 



while (theta2<=180 && theta7<=180)
    
 % if theta2(i) > 360
 % theta2(i) = theta2(i)-360;
 % end





AO4=sqrt(AO2^2+O2_O4^2-2*AO2*O2_O4*cosd(theta2));

beta=acosd((O2_O4^2+AO4^2-AO2^2)/(2*O2_O4*AO4));

psi = acosd((AB^2 + AO4^2 - BO4^2) / (2*AB*AO4));

lamda = acosd((BO4^2 + AO4^2 - AB^2) / (2*BO4*AO4));

theta4=180-lamda-beta;

BO7=sqrt(BO4^2+O4_O7^2-2*BO4*O4_O7*cosd(theta4));

gama= acosd((BO7^2 + O4_O7^2 - BO4^2) / (2*BO7*O4_O7));

delta=180-theta4-gama;

kai=180-gama-theta7;

PB=sqrt(BO7^2+PO7^2-2*BO7*PO7*cosd(kai));





theta3=psi-beta;







 %if theta2 > 180
 %theta3 = psi + beta;
 %theta4 = 180 - lamda + beta;
 %end



 % defining joints positions
 
 O2x = 0;
 O2y = 0;

 O4x= x1;
 O4y = y1;

 O7x = x2;
 O7y = y2;



 Ax = O2x + AO2*cosd(theta2);
 Ay = O2y + AO2*sind(theta2);

 Bx = O4x + BO4*cosd(theta4);
 By = O4y+ BO4*sind(theta4);

 Cx=Ax+AC*cosd(alpha1+theta3);
 Cy=Ay+AC*sind(alpha1+theta3);

 Px=O7x+PO7*cosd(theta7);
 Py=O7y+PO7*sind(theta7);

 % x and y are  coordinates of joint between link 5 ans 6

syms x y;

eq1=(x-Px)^2+(y-Py)^2==DP^2;
eq2=(x-Cx)^2+(y-Cy)^2==CD^2;

sol=solve([eq1,eq2],[x,y]);

%valid_solutions = [sol.x, sol.y];

Dx = double(sol.x);
Dy = double(sol.y);

%disp(valid_solutions);

plot( [O2x Ax], [O2y Ay], [Ax Bx], [Ay By] ...
 , [Ax Cx], [Ay Cy], [Cx Bx], [Cy By],[Bx O4x], [By O4y]...
 ,[Cx Dx], [Cy Dy],[Dx Px], [Dy Py],[Px O7x], [Py O7y], 'LineWidth',3)
 % hold on;
 
 % plot( [Ax(i) Px(i)], [Ay(i) Py(i)],'LineWidth',3)
 % plot( [Bx(i) Px(i)], [By(i) Py(i)],'LineWidth',3)

 grid;
 axis equal;
 axis([-5 15 -5 10]);
 drawnow;
 hold off;

 F=getframe(fig1);
 writeVideo(writerobj, F);


theta2=theta2+2;
theta7=theta7+2;
end

close(writerobj);
disp("Video File Written");











