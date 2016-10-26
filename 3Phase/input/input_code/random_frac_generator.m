

DX=100;
DY=100;
NX=50;
NY=50;
Lx=5000;
Ly=4500;
mu=1800;
phi=100;
n=100;
a=rand(n,2);
a(:,1)=a(:,1)*Lx;
a(:,2)=a(:,2)*Ly;
theta=pi/9*randn(n,1);
length=mu+5*randn(n,1);
for i=1:n
   a(i,3)=a(i,1)+length(i)*cos(theta(i));
   a(i,4)=a(i,2)+length(i)*sin(theta(i));
   if(a(i,3)<0 || a(i,3)>Lx)
      a(i,3)=a(i,3)-2*length(i)*cos(theta(i));
   end
   if(a(i,4)<0 || a(i,4)>Ly)
      a(i,4)=a(i,4)-2*length(i)*sin(theta(i));
   end
end
a(:,5)=2000;
a(:,6)=0.05;
a(:,7)=0.5;
for i=1:n
    line([a(i,1),a(i,3)],[a(i,2),a(i,4)],'LineWidth',2);
    hold on;
end

Dx=10;
Dy=10;
for i = 0:n  %Horizontal Lines
    plot([0 NX*DX], [DY*i DY*i], 'k');  
end 

for i = 0:n   %Vertical Lines
    plot([DX*i DX*i], [0 NY*DY], 'k');    
end
axis([0,5000,0,5000]);