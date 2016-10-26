close all;
figure;
hold on;
title('comparation of injecting BHP');
xlabel('production date (day)');
% ylabel('production rate mcf/day');
ylabel('injecting BHP');
plot(x(:,1),x(:,2));

plot(y(:,1),y(:,2),'d');
legend('eclipse result', 'my_result');

