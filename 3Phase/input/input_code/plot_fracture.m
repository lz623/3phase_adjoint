% 
%  a(10,3)=a(10,3)-600;
%  a(10,1)=a(10,1)-600;
figure
for i=1:size(a,1)
    line([a(i,1),a(i,3)],[a(i,2),a(i,4)],'LineWidth',2);
    hold on;
end