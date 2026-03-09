function f=ellipse(center_x, center_y, R1, R2, clr)

x = center_x - R1;
y = center_y - R2;
h1 = R1*2;
h2 = R2*2;

rectangle('Position',[x,y,h1,h2],...
          'Curvature',[1,1],...
         'LineWidth',2,'LineStyle','-','Edgecolor',clr)
daspect([1,1,1])