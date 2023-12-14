
function GroupsAll=PlotClustering(X, sol)

    % Cluster Centers
    m = sol.Position;
    k = size(m,1);
    
    % Cluster Indices
    ind = sol.Out.ind;
    
    Colors = hsv(k);
    
    for j=1:k
        Xj = X(ind==j,:);
        GroupsAll(j).Coords=Xj;
        plot3(Xj(:,1),Xj(:,2),Xj(:,3),'o','LineWidth',2,'Color',Colors(j,:));
        hold on;
         plot3(m(j,1),m(j,2),m(j,3),'sk','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','k');
  text(m(j,1)+5, m(j,2),m(j,3), ['CH' num2str(j)]);

    end
    
     
%     hold off;
    grid on;
    
end