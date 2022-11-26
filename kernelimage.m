x=linspace(0,1,100);
y_1=exp(-abs(x)/0.001);
y_2=exp(-abs(x)/0.01);
y_3=exp(-abs(x)/0.1);
figure(1)
plot(x,y_1,'LineStyle','-','LineWidth',1.5);
hold on
plot(x,y_2,'LineStyle','-','LineWidth',1.5);
hold on
plot(x,y_3,'LineStyle','-','LineWidth',1.5);
xlabel('$|X_i-X_j|$','Interpreter','latex','FontSize',16);
ylabel('$w_{j,h}(X_i)$','Interpreter','latex','FontSize',16)
legend('$h=0.001$','$h=0.01$','$h=0.1$','Interpreter','latex')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
saveas(gcf,'kernel.pdf');