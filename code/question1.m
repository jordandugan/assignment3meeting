%simulation constants
L = 200e-9;
W = 100e-9;
K = 1.3806e-23;
m = 0.26*9.1093e-31;
q = 1.60217662e-19;

%Calculate force on electrons
V = 0.1;
Fx = -q*V/200e-9

%Place electrons in Boundary
x = L*rand(1000,1);
y = W*rand(1000,1);

% Calculate how many electrons each "particle" in the simulation represents
electronConcentration = 10^15*1e4;

electronsPerParticle = electronConcentration*W*L/1000;

% Assign electron velocity based on Maxwell Boltzman Distribution
T = 300;
std = sqrt(K*T/m);
vth = sqrt(2*K*T/m);
Vx = normrnd(0,std,[1000,1]);
Vy = normrnd(0,std,[1000,1]);
V = sqrt(Vx.^2 + Vy.^2); %Now it is important to check that the velocicties were assigned correctly by confirming that the average velocity is close to the thermal velocity and by plotting a histogram of the velocity.
mean(V)
dt = 100e-9/vth/100;
figure(1)
histogram(V);
xlabel('Velocity')
ylabel('Number of Electrons')

%Initialize variables to plot
Tplot = zeros(300);
Ix = zeros(300);


for i =1:300
    xold = x;
    yold = y;

    % Define region boundaries and rules for interacting with boundaries
    xboundRight = x > L;
    xboundLeft = x < 0;
    ybound = (y > W) | (y <0);
    x(xboundRight) = x(xboundRight) - L;
    x(xboundLeft) = x(xboundLeft) + L;
    xold(xboundRight | xboundLeft) = x(xboundRight | xboundLeft);
    Vy(ybound) = -Vy(ybound);

    %Update Position
    x = x + Vx*dt + 0.5*Fx*dt^2/m;
    y = y + Vy*dt;
    Vx = Vx + Fx*dt/m;

    % Determine Witch electrons scatter and update velocity
    scatter = rand(1000,1) < (1 - exp(-dt/0.2e-12));
    Vx(scatter) = normrnd(0,std,size(Vx(scatter)));
    Vy(scatter) = normrnd(0,std,size(Vy(scatter)));

    xplot = transpose([xold(1:10) x(1:10)]);
    yplot = transpose([yold(1:10) y(1:10)]);
    Tplot(i) = (1/(2*K))*mean(Vx.^2 + Vy.^2)*m;
    
    %Current density and Current Calculations
    Jx = mean(Vx)*1000*electronsPerParticle*(-q);
    Ix(i) = Jx*W*L;
    figure(2)
    subplot(3,1,1);
    plot(xplot,yplot)
    xlim([0 L])
    ylim([0 W])
    title('Electron Trajectory')
    xlabel('x')
    ylabel('y')
    hold on
    subplot(3,1,2)
    plot(Tplot(1:i))
    title('Temperature vs Time Step')
    xlabel('Number of Time Steps')
    ylabel('Temperature (K)')
    hold on
    subplot(3,1,3)
    plot(Ix(1:i))
    title('Current vs Time Step')
    xlabel('Number of Time Steps')
    ylabel('Current (I)')
    hold on
    drawnow


end


% Show electron distribution
figure(3);
hist3([x y],'CdataMode','auto'); 
view(2);
title('Electron Density');
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('Electron Density Heat Map');

temp_sum_x = zeros(20,10);
temp_sum_y = zeros(20,10);
temp_num = zeros(20,10);

for i=1:1000
 % Find which "bin" it belongs in:
 x1 = floor(x(i)/1e-8);
 y1 = floor(y(i)/1e-8);
 if(x1<=0)
 x1 = 1;
 end
 if(y1<=0)
 y1= 1;
 end
 if(y1>10)
     y1 = 10;
 end
 if(x1>20)
     x1=20;
 end
 % Add its velocity components to the cumulative count:
 temp_sum_y(x1,y1) = temp_sum_y(x1,y1) + Vy(i).^2;
 temp_sum_x(x1,y1) = temp_sum_x(x1,y1) + Vx(i).^2;
 temp_num(x1,y1) = temp_num(x1,y1) + 1;

end

temp = (temp_sum_x + temp_sum_y).*m./K./2./temp_num;
temp(isnan(temp)) = 0;
temp = transpose(temp);


figure(4)
%Plot temperature distribution
[X Y] = meshgrid(linspace(0,200,20),linspace(0,100,10));
surf(X,Y,temp)
title('Temperature Map');
xlabel('x (nm)');
ylabel('y (nm)');
zlim([0 1500]);

 

