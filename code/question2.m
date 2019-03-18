clear all

L = 200e-9;
W = 100e-9;
Lb = 40e-9;
Wb = 40e-9;
meshspace = 1e-9;
nx = round(L/meshspace + 1);
ny = round(W/meshspace + 1);

cond1 = 1;
cond2 = 1e-2;


%This section creates a conductivity map witch assigns a conductivity to
%each node in the mesh. This will be used to create the G matrix
condMap = zeros(nx,ny);

for i = 1:nx
   for j = 1:ny
       if (i-1)>0.5*(L-Lb)/meshspace&&(i-1)<0.5*(L+Lb)/meshspace&&((j-1)<Wb/meshspace||(j-1)>(W-Wb)/meshspace)
           condMap(i,j) = cond2;
       else
           condMap(i,j) = cond1;
       end
   end
end

figure(1)
imagesc([0 W],[0 L],condMap);
colorbar
xlabel('y')
ylabel('x')
title('conductivity vs position')

G = sparse(nx*ny);
B = zeros(1,nx*ny);
for i = 1:nx
    for j = 1:ny
        n = j +(i-1)*ny;
        
        %V=1 @ x=0 BC
        if i == 1
            G(n,n) = 1;
            B(n) = 0.1;
        
        %V=0 @ x=L BC
        elseif i == nx 
            G(n,n) = 1;
        
        %Absorbing BC @ y=0
        elseif j == 1 
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nyp = j+1 + (i-1)*ny;
            
            %Resistor Values from conduction map
            rxm = (condMap(i,j) + condMap(i-1,j))/2;
            rxp = (condMap(i,j) + condMap(i+1,j))/2;
            ryp = (condMap(i,j) + condMap(i,j+1))/2;
            
            %node equations from resistor values
            G(n,n) = -(rxm + rxp + ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp; 
        
         %Absorbing BC @ y=W
        elseif j == ny
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nym = j-1 + (i-1)*ny;
            
            %Resistor Values from conduction map
            rxm = (condMap(i,j) + condMap(i-1,j))/2;
            rxp = (condMap(i,j) + condMap(i+1,j))/2;
            rym = (condMap(i,j) + condMap(i,j-1))/2;
            
            %node equations from resistor values
            G(n,n) = -(rxm + rxp + rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
        
        %internal nodes
        else
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            %Resistor Values from conduction map
            rxm = (condMap(i,j) + condMap(i-1,j))/2;
            rxp = (condMap(i,j) + condMap(i+1,j))/2;
            ryp = (condMap(i,j) + condMap(i,j+1))/2;
            rym = (condMap(i,j) + condMap(i,j-1))/2;
            
            %node equations from resistor values
            G(n,n) = -(rxm + rxp + rym + ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp; 
        
        end
    end
end

%solving and plotting
V = G\B';
Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j +(i-1)*ny;
        Vmap(i,j) = V(n);
    end
end

[X, Y] = meshgrid(0:meshspace:L,0:meshspace:W);
figure(2)
surf(X',Y',Vmap)
colorbar
hold on
imagesc([0 L],[0 W],Vmap')
xlabel('x')
ylabel('y')
zlabel('potential')
title('potential vs position')
hold off

[Ey, Ex] = gradient(Vmap,meshspace);
Ex = -Ex;
Ey = -Ey;

figure(3)
quiver(X',Y',Ex,Ey)
xlim([0 L])
ylim([0 W])
xlabel('x')
ylabel('y')
title('Electric Feild')