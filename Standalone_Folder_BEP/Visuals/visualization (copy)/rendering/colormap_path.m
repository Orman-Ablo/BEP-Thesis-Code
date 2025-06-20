 % See if colormap is cyclic

 n_color = 500;
 cmap = hsv(n_color);
 my_colormap = colormap(cmap(1:floor(n_color*5/6),:));
 %my_colormap = cmap;

 figure;
 h = scatter3(my_colormap(:,1),my_colormap(:,2),my_colormap(:,3),150,'filled','MarkerEdgeColor','k');
 xlabel('Red')
 ylabel('Green')
 zlabel('Blue')
 colormap(my_colormap)
 colorbar
 ax=gca;
 ax.FontSize = 24;
 h.CData = my_colormap;
 set(gca,'colormap',my_colormap)
 ax.GridAlpha = 0.3;