clc
clear
% Putting down constants
n=1000;
ax=0; bx=2*pi; ay=-0; by=2*pi;
% Generating x and y values
x=linspace(ax,bx,n);
ix=1:n;
y=linspace(ay,by,n);
iy=1:n;
% The boundary conditions
% 

fa=(y-ay).^2.*cos(y); ga= y.*(y-ay).^2;
ubx=ga;
uax= fa; 
uby = (by-ay).^2.*cos(by) + (x-ax)/(bx-ax)*(by*(by-ay)^2-(by-ay)^2*cos(by));
hx= 2*pi/(n-1); hy=hx; h=hy;
u(:,1)=uax; u(:,n)=ubx;
u(n,:)=uby; 
F = 0;
% Loop time

for j=2:n-1
    for i=2:n-1
        %F(i,j) = sin(pi.*(x(i)-ax)./(bx-ax)).*cos(pi/2*(y(j)-ay)./(by-ay)+1);
        u(i,j)=1/4*(u(i-1,j)+u(i,j-1)+u(i+1,j)+u(i,j+1));%+h^2*F(i,j);

    end
end 
% Ghost node Neumann Conditions du/dy(y=ay)= 0
u(1,:)=u(2,2);
U=u
surf(x,y,U)

