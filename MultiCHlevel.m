
function GroupsAll=MultiCHlevel(X, sol,LevelN_CH)

    % Cluster Centers
    k = size(sol,2);
    
    Colors=[0 0 255;128 255 0;0 255 255;255 0 128]./255;
%     Colors =['bgyc'];
    Colors2=hsv(k);
      
    for j=1:k
        m = sol(j).position;
    % Cluster Indices
    ind = sol(j).sensingNodes;
    
        Xj = X(ind,:);
        GroupsAll(j).Coords=Xj;
         plot3(Xj(:,1),Xj(:,2),Xj(:,3),'ko','LineWidth',2,'MarkerFaceColor',Colors2(j,:));
        for xx=1:size(Xj,1)
            plot3([m(1) Xj(xx,1)],[m(2) Xj(xx,2)],[m(3) Xj(xx,3)],'--','LineWidth',1.2,'Color',Colors2(j,:));
        end
        hold on;
        if nnz(ismember(LevelN_CH(1).pos,j))==1
         plot3(m(1),m(2),m(3),'s','LineWidth',2,'Color',Colors(1,:),'MarkerSize',10,'MarkerFaceColor',Colors(1,:));
        elseif nnz(ismember(LevelN_CH(2).pos,j))==1
         plot3(m(1),m(2),m(3),'s','LineWidth',2,'Color',Colors(2,:),'MarkerSize',10,'MarkerFaceColor',Colors(2,:));
        elseif nnz(ismember(LevelN_CH(3).pos,j))==1
         plot3(m(1),m(2),m(3),'s','LineWidth',2,'Color',Colors(3,:),'MarkerSize',10,'MarkerFaceColor',Colors(3,:));
        elseif nnz(ismember(LevelN_CH(4).pos,j))==1
         plot3(m(1),m(2),m(3),'s','LineWidth',2,'Color',Colors(4,:),'MarkerSize',10,'MarkerFaceColor',Colors(4,:));
        end
  text(m(1)+5, m(2),m(3), ['CH' num2str(j)]);

    end
    
     
%     hold off;
    grid on;
    
end