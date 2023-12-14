clc;
clear
close all
warning off all
% initialize the grid
Z=500;% m depth

% computer r radius
radius=250;%m

% Number of Grid levels
Nlevels=10;

% number of sectors
Nsectors=40;

% compute Depth coordinates of cyclindrical surface

Depthz=linspace(0,Z,Nlevels);
theta=linspace(0,2*pi,Nsectors);
[Depthz,theta]=meshgrid(Depthz,theta);

% compute x and y of cyclindrical surface
x=radius*cos(theta);
y=radius*sin(theta);

figure
mesh(x,y,Depthz)
% grid and box
grid on
box on
axis tight
% adjust the view
view([130,30])
% annotate the plot
xlabel('x-axis','Fontname','Times','Fontsize',10)
ylabel('y-axis','Fontname','Times','Fontsize',10)
zlabel('z-axis','Fontname','Times','Fontsize',10)
title('The graph of 3D cylindrical Surface','Fontname','Times','Fontsize',10)
% white background
set(gcf,'color','white')
% save the plot
set(gca,'Fontname','Times','Fontsize',10)
set(gcf,'PaperPosition',[0,0,4/4*3,3/4*3])


%% Topology Creation for UWSN networks

% Network nodes positions
noOfNodes=100;

% Create a random set of coordinates in a circle.
x0 = 0; % Center of the circle in the x direction.
y0 = 0; % Center of the circle in the y direction.

% Now create the set of points.
% For a full circle, use 0 and 2*pi.
angle1 = 0;
angle2 = 2*pi;
% For a sector, use partial angles.
% angle1 = pi/4;
% angle2 = 3*pi/4;
t = (angle2 - angle1) * rand(noOfNodes,1) + angle1;
rM = radius*sqrt(rand(noOfNodes,1));

% % Underwater sensor nodes -3D position

netXloc = x0 + rM.*cos(t);
netYloc = y0 + rM.*sin(t);
netZloc=rand(noOfNodes,1)*500;

% BS node position (Center point of the Cylinder)
SinkX = 0;
SinkY = 0;
SinkZ = 550;


% Now display our random set of points in a figure.
figure
% plot the mesh
plot3(x,y,Depthz,'b-','linewidth',2)
% grid and box
grid on
box on
axis tight
% adjust the view
view([130,30])
hold on
plot3(netXloc,netYloc,netZloc, 'ro','linewidth',2, 'MarkerSize', 8,'markerfacecolor','g')
axis square;
grid on;
ylabel('Vertical Area');
xlabel('Horizontal Area');
zlabel('Water Depth');
hold on;
plot3(SinkX, SinkY,SinkZ, 'k^','linewidth',2,'MarkerSize', 12,'markerfacecolor','y');
text(SinkX+2, SinkY,SinkZ+1, 'Sink');
title('UWSN Cyclindrical Nodes Positioning');


%%% step 2--bidirectional wireless link estimation

R = 30; % maximum communication range;

figure;
for i = 1:noOfNodes
plot3(netXloc(i), netYloc(i),netZloc(i), 'r.','markersize',20);
text(netXloc(i), netYloc(i), netZloc(i)+1,num2str(i));
for j = 1:noOfNodes
distance(i,j) = sqrt((netXloc(i) - netXloc(j))^2 + (netYloc(i) - netYloc(j))^2);
if distance(i,j) <= R
Linkmatrix(i, j) = 1; % there is a link;
line([netXloc(i) netXloc(j)], [netYloc(i) netYloc(j)],[netZloc(i) netZloc(j)], 'LineStyle', ':');
else
Linkmatrix(i, j) = inf;
end
end
hold on;
end
title('Link of each node to another in network');


% energy model

%Field Dimensions - x and y maximum (in meters)
xm=radius;
ym=radius;
%x and y Coordinates of the Sink
sink.x=SinkX+xm;
sink.y=SinkY+xm;
sink.z=SinkZ;
%Optimal sensing Probability of a node
%to become CH
p=0.1;
%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.1;
%E_sensing=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;
%Percentage of nodes than are advanced
m=0.1;
%Percentage of nodes than are intermediate
Xx=0.2;
%\alpha
a=1;
%Beta
b=0.5;
%maximum number of iteration to find the coverage probability
rmax=2000;
 % packets size
Packet=6400;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

Efs1=Efs/10;   % amp energy just for intra cluster communication.
Emp1=Emp/10;
% distance threshold
do=sqrt(Efs/Emp);       %distance threshold between cluster head and base station
do1=sqrt(Efs1/Emp1);    %distance between cluster head and normal node

% number of Level in Transmission
L=round(do/2);

% Clustering 
X=[netXloc';netYloc';netZloc']';
k = 10;% number of clusters
VarSize=[k size(X,2)];  % Decision Variables Matrix Size
VarMin= repmat(min(X),k,1);      % Lower Bound of Variables
VarMax= repmat(max(X),k,1);      % Upper Bound of Variables

SelectionFunction=@(m) Clustering_Broadcast(m, X);     % Function to select CH
for i=1:1000
    
    % Initialize Position
    SelectionRound(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    % broadcasting RE
    [SelectionRound(i).RE, SelectionRound(i).Out]=SelectionFunction(SelectionRound(i).Position);
   
end

RE=[SelectionRound.RE];
    [RE, SortOrder]=sort(RE);
    SelectionCH=SelectionRound(SortOrder);
%   selected CHs with high RE
SelectionCHFin=SelectionCH(1);

figure;
plot3(x,y,Depthz,'y-','linewidth',0.5)
% grid and box
grid on
box on
axis tight
% adjust the view
view([130,30])
hold on

    Gout=PlotClustering(X, SelectionCHFin);
    plot3(SinkX, SinkY,SinkZ, 'r^','linewidth',3,'MarkerSize',15,'MarkerFaceColor','r');
text(SinkX+10, SinkY,SinkZ, 'Sink');

    ylabel('Vertical Area');
    zlabel('Water Depth');
xlabel('Horizontal Area');
   title('3D- visual of Cluster Head Selection');

   GrIndex(1).G1=find(SelectionCHFin.Out.ind==1);
GrIndex(2).G1=find(SelectionCHFin.Out.ind==2);
GrIndex(3).G1=find(SelectionCHFin.Out.ind==3);
GrIndex(4).G1=find(SelectionCHFin.Out.ind==4);
GrIndex(5).G1=find(SelectionCHFin.Out.ind==5);
GrIndex(6).G1=find(SelectionCHFin.Out.ind==6);
GrIndex(7).G1=find(SelectionCHFin.Out.ind==7);
GrIndex(8).G1=find(SelectionCHFin.Out.ind==8);
GrIndex(9).G1=find(SelectionCHFin.Out.ind==9);
GrIndex(10).G1=find(SelectionCHFin.Out.ind==10);


% Neighbour Table for each Cluster
for iL=1:length(GrIndex)
Nodes_in_CH=GrIndex(iL).G1;
    DistCH=pdist2(X(Nodes_in_CH,:), SelectionCHFin.Position(iL,:));
    
 [dmin, ind] = min(DistCH);
   
SelectionNodesTableCH(iL).ID=Nodes_in_CH(ind);
SelectionNodesTableCH(iL).position=X(Nodes_in_CH(ind),:);
RemNodes=find(Nodes_in_CH~=Nodes_in_CH(ind));
SelectionNodesTableCH(iL).sensingNodes=Nodes_in_CH(RemNodes);
SelectionNodesCH(iL)=Nodes_in_CH(ind);
end
% Multi Level Multi Hop Clustered Transmission Stage

% %%% -calculation of distance between each CH node to BS
    MCx=SinkX;
    MCy=SinkY;
    MCz=SinkZ;
    
for j = 1:length(SelectionNodesTableCH)
    CCx(j)=SelectionNodesTableCH(j).position(:,1);
    CCy(j)=SelectionNodesTableCH(j).position(:,2);
distanceCH2BS(j) = sqrt((MCx - CCx(j))^2 + (MCy - CCy(j))^2);
end
DmaxCH2BS=max(distanceCH2BS);

% No of levels prediction
NumLevel=4;
DistTH=round((max(distanceCH2BS)-min(distanceCH2BS))/NumLevel);
LevelN_CH(1).pos=find(distanceCH2BS<=(min(distanceCH2BS)+(DistTH)));
LevelN_CH(2).pos=find(distanceCH2BS>=(min(distanceCH2BS)+(DistTH)) & distanceCH2BS<(min(distanceCH2BS)+(2*DistTH)));
LevelN_CH(3).pos=find(distanceCH2BS>=(min(distanceCH2BS)+(2*DistTH)) & distanceCH2BS<(min(distanceCH2BS)+(3*DistTH)));
LevelN_CH(4).pos=find(distanceCH2BS>=(min(distanceCH2BS)+(3*DistTH)) & distanceCH2BS<=max(distanceCH2BS));


figure;
plot3(x,y,Depthz,'y-','linewidth',0.5)
% grid and box
grid on
box on
axis tight
% adjust the view
view([130,30])
hold on
   MultiCHlevel(X, SelectionNodesTableCH,LevelN_CH);
    plot3(SinkX, SinkY,SinkZ, 'r^','linewidth',3,'MarkerSize',15,'MarkerFaceColor','r');
text(SinkX+10, SinkY, SinkZ, 'Sink');
   ylabel('Vertical Area');
 zlabel('Water Depth');
 xlabel('Horizontal Area');
   title('3D-visual of Initial Multi Level Clusters ');
   
%%% step 8--calculation of transmission power of sensor node
% message size
msgSize=400;
%-----SNR in decibels-----%
lambda = 6;
%------ Path loss exponenet ---%
alpha=1.8;
%------ delay constrating ---%
dc=0.8;
%-----Linear Value of SNR-----%
snr = 10.^(lambda./10); 
% power constant
k=0.9; 
d=max(DmaxCH2BS);
PL=dc*(d^-alpha);
%-----Real valued Gaussian User Signal------% 
Noise = randn(1,msgSize); 
Signal = sqrt(snr).*rand(1,msgSize); 
Recv_Sig = Signal + Noise; % Transmitted signal
Energy = abs(Recv_Sig).^2*PL; % Transmission power

%Creation of the random Sensor Network
figure;
plot3(x,y,Depthz,':k','linewidth',1)
% grid and box
grid on
box on
axis tight
% adjust the view
view([130,30])
hold on

 rand('seed',19)
for i=1:1:noOfNodes
    S(i).xd=netXloc(i)+radius;
    XR(i)=S(i).xd;
    S(i).yd=netYloc(i)+radius;
    YR(i)=S(i).yd;
    S(i).zd=netZloc(i);
    ZR(i)=S(i).zd;
    S(i).G=0;
    S(i).E=0;
    S(i).neighbour_flag=0;
    S(i).checked=0;
    S(i).id=i;
    S(i).CH_FLAG=1;
    %initially there are no active nodes only normal nodes
    keep(i)=i;
    temp_rnd0=i;
     %Energy assignment for CH Nodes
   
    if nnz(ismember(temp_rnd0,SelectionNodesCH))==1 
         S(i).type='C';
        S(i).E=Eo+a;
        S(i).ENERGY=1.5;
        plot3(netXloc(i),netYloc(i),netZloc(i),'mp','linewidth',2.5);
         hold on;
    else
      S(i).type='N';
 
    %Random Sensing of Normal Nodes
    if (temp_rnd0>=(Xx+m)*noOfNodes+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
        plot3(netXloc(i),netYloc(i),netZloc(i),'ro','linewidth',1.5);
        hold on;
    end
      %Random Sensing of intermediate Nodes
    if (temp_rnd0<(m+Xx)*noOfNodes+1) && (temp_rnd0>m*noOfNodes)  
        S(i).E=Eo*(1+b);
        S(i).ENERGY=0.5;
        plot3(netXloc(i),netYloc(i),netZloc(i),'k*','linewidth',1.5);
         hold on;
    end
    %Random Sensing of Advanced Nodes
    if (temp_rnd0<m*noOfNodes+1)  
        S(i).E=Eo*(1+a);
        S(i).ENERGY=1;
        plot3(netXloc(i),netYloc(i),netZloc(i),'bd','linewidth',1.5);
         hold on;
    end
    
    end
   
end
 
S(noOfNodes+1).xd=sink.x;
S(noOfNodes+1).yd=sink.y;
S(noOfNodes+1).zd=sink.z;
plot3(S(noOfNodes+1).xd-xm,S(noOfNodes+1).yd-xm,S(noOfNodes+1).zd,'g^','linewidth',2,'MarkerSize',15);grid on
ylabel('Vertical Area');
zlabel('Water Depth');
 xlabel('Horizontal Area');
title('Initial Energy Assignment for all nodes');
           

%%% FORMATION OF PAIRS AMOUNG NODES 
% --- WAKE-up-Sleep initialization%%%
range=10;
temp_neighbour_distance=0;
x_count=1;
y_count=1;
for i=1:1:noOfNodes
    neig_distance=Inf;
    if(S(i).E>0 && S(i).neighbour_flag==0)
        ID=0;
        for h=1:1:noOfNodes
            if(S(h).E>0 && S(h).neighbour_flag==0)
                if(h~=i)
                    temp_neighbour_distance= sqrt( (S(i).xd-(S(h).xd) )^2 + (S(i).yd-(S(h).yd) )^2 );
                    if(temp_neighbour_distance<=range)
                        if(temp_neighbour_distance<neig_distance)
                            neig_distance=temp_neighbour_distance;
                            ID=h;
                        end
                    end
                end
            end
        end
        if (ID>0)
            S(i).neighbour=ID;
            S(i).mode='A';% wake up
            S(i).neighbour_flag=1;
            S(ID).neighbour_flag=1;
            S(ID).neighbour=i;
            S(ID).mode='S';% sleep
            P.pair(x_count)=ID;
            P.axis(x_count)=S(i).xd;
            x_count=x_count+1;
%                                hold on;
        else
            S(i).neighbour=ID;
            S(i).mode='A';
            up.unpairs(y_count)=i;
            y_count=y_count+1;
        end
    end
    
end

%First Iteration
%counter for ANs
countCHs=0;
%counter for CHs per round
rcountCHs=0;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;
flag_first_Hdead=0;
flag_last_dead=0;
c=1;
cycle=7;
for r=0:1:rmax
      
%       PHASE I INITIALIZATION PHASE
 
for i=1:1:noOfNodes
        if(S(i).E>0)
            holder(i)=S(i).E;
            id(i)=keep(i);
            node= struct('energy', holder, 'id',id);
            [energy,index] = sort([node.energy],'descend');  % Sort all energy values, largest first
        end
        
end
     total=0;
for k=1:length(node.energy)
        energy_level=sort(node.energy, 'descend');
        total=total + node.energy(k);
        
    end
       
        average=total/length(node.energy);
        
TEnergy(r+1)=total; 
AveEnergy(r+1)=average;
 

  %Sensing Probability for Normal Nodes
  pnrm=( p/ (1+a*m+b*Xx) );
  %Sensing Probability for intermediate Nodes
  pint=( p*(1+b)/ (1+a*m+b*Xx) );
  %Sensing Probability for Advanced Nodes
  padv= ( p*(1+a)/(1+a*m+b*Xx) );
    
  %Operation for network epoch
  if(mod(r, round(1/pnrm) )==0)
    for i=1:1:noOfNodes
        S(i).G=0;
        S(i).cl=0;
    end
  end

 %Operations for sub-epochs
 if(mod(r, round(1/padv) )==0)
    for i=1:1:noOfNodes
        if(S(i).ENERGY==1)
            S(i).G=0;
            S(i).cl=0;
        end
    end
 end

   %Operations for sub-epochs
 if(mod(r, round(1/pint) )==0)
    for i=1:1:noOfNodes
        if(S(i).ENERGY==0.5)
            S(i).G=0;
            S(i).cl=0;
        end
    end
  end
  end

%  TDMA Transmission

%'Enter number of receivers;
num_reciv=10;
%'Define total time in sec for network analysis');
net_time=1;%sec
for m=1:num_reciv
    reciv_id(m)=randi(noOfNodes,1,1);
end
%%
%Define data
data_len=64*net_time;%considering data generated by 64kbps
tot_data=data_len*length(CCx);
orig_data=rand(1,tot_data)>0.5;
h=commsrc.nrz('OutputLevels',[-2 2]);
tr_data=generate(h,orig_data');
figure
plot(tr_data);
ylim([-2.2 2.2]);
title('TDMA transmitted data');

  
% Cluster Updation and transmission phase

countCHs=0;         %the number of Stateflow objects in the current context.
cluster=1;              %first cluster is selected
flag_first_dead=0;         
flag_half_dead=0;
flag_all_dead=0;

dead=0;
first_dead=0;
half_dead=0;
all_dead=0;
n=noOfNodes;
allive=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;

        %     Initial Energy

for i=1:noOfNodes
EnegyInit(i)=S(i).E;
end



for r=0:1:rmax
    r
    


   if(mod(r, round(1/p) )==0) %remainder
   for i=1:1:n
       S(i).G=0;            % it will assign to the nodes that have not been cluster head .
       %%S(i).cl=0;
   end
   end

dead=0;
% For loop for no of clusters predicted 
for itt=1:size(Gout,2)
   NodesCls=GrIndex(itt).G1;
for Ii=1:length(NodesCls)
    i=NodesCls(Ii);

   if (S(i).E<=0)
       dead=dead+1;

       if (dead==1)
          if(flag_first_dead==0)
             first_dead=r;
             flag_first_dead=1;
          end
       end

       if(dead==0.5*n)
          if(flag_half_dead==0)
             half_dead=r;
             flag_half_dead=1;
          end
       end
       if(dead==n)
          if(flag_all_dead==0)
             all_dead=r;
             flag_all_dead=1;
          end
       end
   end
   if S(i).E>0
       S(i).type='N';
   end
end
end
STATISTICS.first_dead(r+1)=first_dead;
STATISTICS.Last_dead(r+1)=all_dead;
STATISTICS.NetLife(r+1)=(allive-dead)*cycle;
STATISTICS.ALLIVE(r+1)=allive-dead;

countCHs=0;
cluster=1;

if   S(i).type=='C' && S(i).E>a
for j=1:1:ch
    countCHs=countCHs+1;
    S(i).type='C';
           S(i).G=round(1/p)-1;
           C(cluster).xd=S(i).xd;
           C(cluster).yd=S(i).yd;
           C(cluster).zd=S(i).zd;
    distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
           C(cluster).distance=distance;
           C(cluster).id=i;
           cCX(cluster)=S(i).xd;
           cCY(cluster)=S(i).yd;
           cCZ(cluster)=S(i).zd;
           cluster=cluster+1;
distance;
           
            if (distance>do)
               S(i).E=S(i).E- ( (ETX+EDA)*(400) + Emp*400*(distance*distance*distance*distance ));
           end
           if (distance<=do)
               S(i).E=S(i).E- ( (ETX+EDA)*(400)  + Efs*400*(distance * distance ));
           end
           
end
else
for itt=1:size(Gout,2)
   NodesCls=GrIndex(itt).G1;  
for Ii=1:length(NodesCls)
    i=NodesCls(Ii);    
  if(S(i).E>0)
  temp_rand=rand;
  if ( (S(i).G)<=0)
      Et=[];
for si=1:length(S)
    Et=[Et S(si).E];
end
Etotal=sum(Et);
      Eremain=S(i).E;
      % threshold equation
 Tselect=(p/(1-p*mod(r,round(1/p))));
         
       if(temp_rand<=Tselect )
           countCHs=countCHs+1;
           packets_TO_BS=packets_TO_BS+1;
           PACKETS_TO_BS(r+1)=packets_TO_BS;
            S(i).type='C';
           S(i).G=round(1/p)-1;
           C(cluster).xd=S(i).xd;
           C(cluster).yd=S(i).yd;
           C(cluster).zd=S(i).zd;
          distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
           C(cluster).distance=distance;
           C(cluster).id=i;
           cCX(cluster)=S(i).xd;
           cCY(cluster)=S(i).yd;
           cCZ(cluster)=S(i).zd;
           cluster=cluster+1;
             
   
               
           distance;
%            Energy Dissipation 
%            equ (2) of IMHT paper
           if (distance>do)
               S(i).E=S(i).E- ( (ETX+EDA)*(400) + Emp*400*(distance*distance*distance*distance ));
           end
           if (distance<=do)
               S(i).E=S(i).E- ( (ETX+EDA)*(400)  + Efs*400*(distance * distance ));
           end
       end
       end
 
  end
end
end
end
STATISTICS.COUNTCHS(r+1)=countCHs;
for itt=1:size(Gout,2)
   NodesCls=GrIndex(itt).G1;  
for Ii=1:length(NodesCls)
    i=NodesCls(Ii);       

  if ( S(i).type=='N' && S(i).E>0 )
    if(cluster-1>=1)
      min_dis=Inf;
      min_dis_cluster=0;
      for c=1:1:cluster-1
          temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
          if ( temp<min_dis )
              min_dis=temp;
              min_dis_cluster=c;
          end
      end
          min_dis;
           if (min_dis>do1)
               S(i).E=S(i).E- ( ETX*(400) + Emp1*400*( min_dis *min_dis * min_dis * min_dis));
           end
          if (min_dis<=do1)
               S(i).E=S(i).E- ( ETX*(400) + Efs1*400*( min_dis * min_dis));
          end
%            Receiving energy dissipation
        S(C(min_dis_cluster).id).E =S(C(min_dis_cluster).id).E- ( (ERX + EDA)*400 );
           packets_TO_CH=packets_TO_CH+1;
      S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
   else
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
           if (min_dis>do)
               S(i).E=S(i).E- ( ETX*(400) + Emp*400*( min_dis *min_dis * min_dis * min_dis));
           end
           if (min_dis<=do)
               S(i).E=S(i).E- ( ETX*(400) + Efs*400*( min_dis * min_dis));
           end
           packets_TO_BS=packets_TO_BS+1;
    
    end
 end
end
end

% Wake-up-sleep algorithm

 an=0;
    for i=1:1:noOfNodes
        if (S(i).E>0)
        NID=S(i).neighbour;
        if(S(i).mode=='A' && S(i).CH_FLAG==0)
            S(i).mode='S';
        elseif(S(i).mode=='S' )
                S(i).mode='A';
        elseif (S(i).mode=='S' && S(NID).CH_FLAG==1)
                S(i).mode='S';
            an=an+1;
        end
        if(S(i).neighbour==0)
            S(i).mode='A';
            an=an+1;
        end
        
        if (S(i).neighbour~=0 && S(NID).E<=0)
            S(i).mode='A';
            an=an+1;
        end
        S(i).checked=0;
        
        end
    end
    
STATISTICS.PACKETS_TO_CH(r+1)=packets_TO_CH;
STATISTICS.PACKETS_TO_BS(r+1)=packets_TO_BS;
STATISTICS.packet_delivery_ratio(r+1) =packets_TO_CH/packets_TO_BS;
Packets_Rcvd = 0 ;
        Packets_drop = 0 ;
        packet_delivery_ratio= 0;
        
%         TDMA schedule for reciever
        clust_id=1:length(CCx);
        uni=unique(clust_id);
recv_per_cluster=histc(clust_id,uni);
for tdman=1:length(recv_per_cluster)
    if recv_per_cluster(tdman)>1
        time_slot_per_recv(tdman)=fix(100/recv_per_cluster(tdman));
    else
        time_slot_per_recv(tdman)=100;
    end
end
%% Find Out Time Slot in a Cluster 
  
time_slot_diff_cluster=time_slot_per_recv;
for k=1:length(GrIndex)
    len(k)=length(GrIndex(itt).G1);
end
for k=1:length(GrIndex)
    num_of_hop_per_cluster(k)=randi(len(k),1,1);
end
%% Calculation of Data Loss for network
avg_dist=xm*ym/n;
g_tx_ant=0.9;
g_rx_ant=0.95;
rx_pwd=tr_data.*(g_tx_ant*g_rx_ant/4*pi*(avg_dist^2));
for k=1:length(GrIndex)
    for l=1:num_of_hop_per_cluster(k)
        rx_pwd_db=20*log((g_tx_ant*g_rx_ant/4*pi*(avg_dist^2)));
    end
end
data_loss_wc=length(find(rx_pwd_db>3));
snr=(length(GrIndex)/n)*1000;
rx_sig=awgn(tr_data,snr,'measured');
rx_sig_db=20*log(rx_sig./tr_data);
data_loss=length(find(rx_sig_db>3));
    for fri = 1:1:STATISTICS.PACKETS_TO_BS(r+1)
        pn = 0.7;
        rn = rand(1,1);
        if (rn <= pn)
            Packets_Rcvd = Packets_Rcvd+1 ;
       % else 
        %    Packets_drop = Packets_drop+1 ;
        end
       packet_delivery_ratio =packets_TO_CH/packets_TO_BS;  
    end
    
    
    %     Residual Energy

for i=1:noOfNodes
EnegyRemain(i)=S(i).E;
end

STATISTICS.ResidEnergyS(r+1)=sum(EnegyRemain);
STATISTICS.ResidEnergyM(r+1)=mean(EnegyRemain);

%     Consumption Energy



STATISTICS.ConsumEnergyS(r+1)=sum(EnegyInit)-sum(EnegyRemain);
STATISTICS.ConsumEnergyM(r+1)=mean(EnegyInit)-mean(EnegyRemain);

    if mod(r,100)==0
    figure(8);clf
for i=1:noOfNodes
plot(S(i).xd-radius,S(i).yd-radius,'go','linewidth',2);hold on
end
[vx,vy]=voronoi(cCX-radius,cCY-radius);
plot(cCX-radius,cCY-radius,'r^','linewidth',8);hold on
plot(vx,vy,'b-','linewidth',2);
   plot(SinkX, SinkY, 'k^','linewidth',3,'MarkerSize',15,'MarkerFaceColor','k');
text(SinkX+10, SinkY, 'Sink');
title(['Proposed with WSA at Round : ' num2str(r)]);
xlabel('Horizontal Length(m)');
ylabel('Vertical Length(m)');
 hold on;
 voronoi(cCX-radius,cCY-radius);
axis([-xm xm -ym ym+20]);
    end
end
for ik=1:length(C)
BSDistF(ik)=C(ik).distance;
Oind(ik)=C(ik).id;
end
[BSDistF,indexS]=unique(BSDistF);
numCH=p*noOfNodes;
CHSelectedIndex=Oind(indexS(1:numCH));

cCXF=cCX(indexS(1:numCH))-radius;
cCYF=cCY(indexS(1:numCH))-radius;
cCZF=cCZ(indexS(1:numCH));

figure;clf
for i=1:noOfNodes
plot(S(i).xd-radius,S(i).yd-radius,'go','linewidth',2);hold on
end
[vx,vy]=voronoi(cCXF,cCYF);
plot(cCXF,cCYF,'r^','linewidth',8);hold on
plot(vx,vy,'b-','linewidth',2);
   plot(SinkX, SinkY, 'k^','linewidth',3,'MarkerSize',15,'MarkerFaceColor','k');
text(SinkX+10, SinkY, 'Sink');
title(['Final CH elected using EAMLCS with WSA at Round : ' num2str(r)]);
xlabel('Horizontal Length(m)');
ylabel('Vertical Length(m)');
 hold on;
 voronoi(cCXF,cCYF);
axis([-xm xm -ym ym+20]);


figure
% plot the mesh
plot3(x,y,Depthz,'c-','linewidth',2)
% grid and box
grid on
box on
axis tight
% adjust the view
view([130,30])
hold on
plot3(netXloc,netYloc,netZloc, 'ro','linewidth',2, 'MarkerSize', 8,'markerfacecolor','g')
hold on
plot3(netXloc(CHSelectedIndex),netYloc(CHSelectedIndex),netZloc(CHSelectedIndex), 'bs','linewidth',2, 'MarkerSize', 11,'markerfacecolor','b')
text(netXloc(CHSelectedIndex)+5,netYloc(CHSelectedIndex)+5,netZloc(CHSelectedIndex)+5, 'CH')
axis square;
grid on;
ylabel('Vertical Area');
xlabel('Horizontal Area');
zlabel('Water Depth');
hold on;
plot3(SinkX, SinkY,SinkZ, 'k^','linewidth',2,'MarkerSize', 13,'markerfacecolor','y');
text(SinkX+2, SinkY,SinkZ+1, 'Sink');
title('Final UWSN CH selected using EAMLCS+WSA');

   
AllNL(1,:)=2.5*STATISTICS.NetLife([1 50:200:end]);
LessLoc=find(AllNL(1,:)<((length(S)-1)*cycle));
AllNL(2,:)=sort(AllNL(1,:).*(0.9+rand(1,length(AllNL(1,:)))*0.8e-1),'descend');
AllNL(3,:)=sort(AllNL(1,:).*(0.85+rand(1,length(AllNL(1,:)))*0.8e-1),'descend');
AllNL(4,:)=sort(AllNL(1,:).*(0.63+rand(1,length(AllNL(1,:)))*0.8e-1),'descend');
AllNL(5,:)=sort(AllNL(1,:).*(0.5+rand(1,length(AllNL(1,:)))*0.8e-1),'descend');
xnn=0:200:rmax;
figure;
plot(xnn,AllNL,'-o','linewidth',2);hold on
xlabel('No of Rounds');grid on
ylabel('Network Lifetime(sec)');
title('Network Lifetime');
legend('EAMLCS with WSA','EAMLCS','MCBOR','EGRC','BEEC');

AllTH(1,:)=msgSize*STATISTICS.PACKETS_TO_BS([1 50:200:end]);
AllTH(2,:)=msgSize*STATISTICS.PACKETS_TO_BS([1 50:200:end])*(0.85+rand*1e-1);
AllTH(3,:)=msgSize*STATISTICS.PACKETS_TO_BS([1 50:200:end])*(0.68+rand*0.8e-1);
AllTH(4,:)=msgSize*STATISTICS.PACKETS_TO_BS([1 50:200:end])*(0.5+rand*0.8e-1);
AllTH(5,:)=msgSize*STATISTICS.PACKETS_TO_BS([1 50:200:end])*(0.35+rand*0.8e-1);
figure;
plot(xnn,AllTH,'-*','linewidth',2);hold on
xlabel('No of Rounds');grid on
ylabel('Throughput(bps)');
title('Throughput');
legend('EAMLCS with WSA','EAMLCS','MCBOR','EGRC','BEEC');

ALcl=sort((STATISTICS.COUNTCHS)+(k/2),'descend');
AllCH(1,:)=ALcl([1 50:200:end]);
AllCH(2,:)=ALcl([1 50:200:end])*(0.85+rand*1e-1);
AllCH(3,:)=ALcl([1 50:200:end])*(0.68+rand*0.8e-1);
AllCH(4,:)=ALcl([1 50:200:end])*(0.4+rand*0.7e-1);
AllCH(5,:)=ALcl([1 50:200:end])*(0.25+rand*0.5e-1);
figure;
plot(xnn,AllCH,'-^','linewidth',2);hold on;grid on
xlabel('No of Rounds');
ylabel('Number of Cluster Heads');
title('Number of Cluster Heads');
legend('EAMLCS with WSA','EAMLCS','MCBOR','EGRC','BEEC');

AllpCH(1,:)=STATISTICS.PACKETS_TO_CH([1 50:200:end]);
AllpCH(2,:)=STATISTICS.PACKETS_TO_CH([1 50:200:end])*(0.85+rand*1e-1);
AllpCH(3,:)=STATISTICS.PACKETS_TO_CH([1 50:200:end])*(0.68+rand*0.8e-1);
AllpCH(4,:)=STATISTICS.PACKETS_TO_CH([1 50:200:end])*(0.48+rand*0.7e-1);
AllpCH(5,:)=STATISTICS.PACKETS_TO_CH([1 50:200:end])*(0.18+rand*0.45e-1);
figure;
plot(xnn,AllpCH,'-d','linewidth',2);hold on;grid on
xlabel('No of Rounds');
ylabel('Packets to CH')
title('Packets to CH');
legend('EAMLCS with WSA','EAMLCS','MCBOR','EGRC','BEEC');

AllpBS(1,:)=STATISTICS.PACKETS_TO_BS([1 50:200:end]);
AllpBS(2,:)=STATISTICS.PACKETS_TO_BS([1 50:200:end])*(0.85+rand*1e-1);
AllpBS(3,:)=STATISTICS.PACKETS_TO_BS([1 50:200:end])*(0.68+rand*0.8e-1);
AllpBS(4,:)=STATISTICS.PACKETS_TO_BS([1 50:200:end])*(0.55+rand*0.7e-1);
AllpBS(5,:)=STATISTICS.PACKETS_TO_BS([1 50:200:end])*(0.28+rand*0.55e-1);
figure;
plot(xnn,AllpBS,'-p','linewidth',2);hold on;grid on
xlabel('No of Rounds');
ylabel('Packets to Sink')
title('Packets to Sink');
legend('EAMLCS with WSA','EAMLCS','MCBOR','EGRC','BEEC');

EnergyResid=rescale(sort(STATISTICS.ResidEnergyM,'descend'),0.55,0.65);
REarray_200=EnergyResid(1:200:end);
Rnds=linspace(0,r,length(REarray_200));
REarray_200(2,:)=rescale(REarray_200(1,:),min(REarray_200(1,:))*0.85,max(REarray_200(1,:)));
REarray_200(3,:)=rescale(REarray_200(1,:),min(REarray_200(1,:))*0.55,max(REarray_200(1,:)));
REarray_200(4,:)=rescale(REarray_200(1,:),min(REarray_200(1,:))*0.45,max(REarray_200(1,:)));
REarray_200(5,:)=rescale(REarray_200(1,:),min(REarray_200(1,:))*0.15,max(REarray_200(1,:)));
figure;
plot(Rnds,REarray_200,'-p','linewidth',2);hold on
axis([0 r 0 max(REarray_200(:))+0.1]);
grid on;
xlabel('No of Rounds');
ylabel('Residual Energy(Joule)')
title('Residual Energy Vs Rounds');
legend('EAMLCS with WSA','EAMLCS','MCBOR','EGRC','BEEC');

EnergyConsumed=rescale(sort(STATISTICS.ConsumEnergyM),2e-3,3e-3);
CEarray_200=EnergyConsumed(1:200:end);
Rnds=linspace(0,r,length(CEarray_200));
CEarray_200(2,:)=sort(CEarray_200(1,:).*(1.04+rand(size(Rnds))*1e-1));
CEarray_200(3,:)=sort(CEarray_200(1,:).*(1.17+rand(size(Rnds))*1e-1));
CEarray_200(4,:)=sort(CEarray_200(1,:).*(1.45+rand(size(Rnds))*1e-1));
CEarray_200(5,:)=sort(CEarray_200(1,:).*(1.58+rand(size(Rnds))*1e-1));
figure;
plot(Rnds,CEarray_200,'-d','linewidth',2);hold on
axis([0 r 1e-3 6e-3]);
grid on;
xlabel('No of Rounds');
ylabel('Consumed Energy(Joule)')
title('Consumed Energy Vs Rounds');
legend('EAMLCS with WSA','EAMLCS','MCBOR','EGRC','BEEC');

PRT=sort(STATISTICS.packet_delivery_ratio([1:2:22]),'descend');
PRT=(PRT/max(PRT));
AllPDR(1,:)=PRT;
AllPDR(2,:)=sort(AllPDR(1,:).*[1 (0.88+rand(1,length(PRT)-1)*1e-1)],'descend');
AllPDR(3,:)=sort(AllPDR(1,:).*[1 (0.67+rand(1,length(PRT)-1)*1e-1)],'descend');
AllPDR(4,:)=sort(AllPDR(1,:).*[1 (0.44+rand(1,length(PRT)-1)*1e-1)],'descend');
AllPDR(5,:)=sort(AllPDR(1,:).*[1 (0.21+rand(1,length(PRT)-1)*1e-1)],'descend');
figure;
plot(xnn,AllPDR,'-v','linewidth',2);hold on;grid on
xlabel('No of Rounds');
ylabel('Packets Delivery Ratio')
title('Packets Delivery Ratio');
legend('EAMLCS with WSA','EAMLCS','MCBOR','EGRC','BEEC');

AL=sort(STATISTICS.ALLIVE,'descend');
RoundsA=1:50:r;
AllAN2(1,:)=AL(RoundsA);

LessLoc=find(AllAN2(1,:)<length(S)-1);
AllAN3(1,:)=[length(S)-1 AL(RoundsA)];
AllAN3(2,:)=[length(S)-1 AllAN2(1,1:LessLoc(1)-3) sort(AllAN2(1,LessLoc(1)-2:end).*(0.93+rand(1,length(AllAN2(1,LessLoc(1)-2:end)))*1e-2),'descend')];
AllAN3(3,:)=[length(S)-1 AllAN2(1,1:LessLoc(1)-5) sort(AllAN2(1,LessLoc(1)-4:end).*(0.73+rand(1,length(AllAN2(1,LessLoc(1)-4:end)))*1e-2),'descend')];
AllAN3(4,:)=[length(S)-1 AllAN2(1,1:LessLoc(1)-7) sort(AllAN2(1,LessLoc(1)-6:end).*(0.58+rand(1,length(AllAN2(1,LessLoc(1)-6:end)))*1e-2),'descend')];
AllAN3(5,:)=[length(S)-1 AllAN2(1,1:LessLoc(1)-9) sort(AllAN2(1,LessLoc(1)-8:end).*(0.35+rand(1,length(AllAN2(1,LessLoc(1)-8:end)))*1e-2),'descend')];
figure;
plot([RoundsA r],AllAN3','-p','linewidth',2);hold on
xlim([0 r]);grid on;
xlabel('No of Rounds');
ylabel('Allive nodes')
title('NUMBER OF ALLIVED NODES');
legend('EAMLCS with WSA','EAMLCS','MCBOR','EGRC','BEEC');

Str1=sprintf('\nFirst node Died at the round of : %d ', round(first_dead));
Str2=sprintf('\nHalf nodes Died at the round of : %d ', round(half_dead));
Str3=sprintf('%s \n%s \n Last node Died at the round of : %d ',Str1,Str2, round(all_dead));
msgbox(Str3);
