%function [leadc] = simulation_example
    nn = 50;
    timeinit = 1000; % number of time steps to initialize probabilities.
    timereal = 2000;  % number of time steps to do real simulation.
    simu.duration = timeinit + timereal; % number of time steps
    simu.number_of_nodes = [25 25];
    simu.alpha = [0 100];   
    simu.grid_size = 3; 
    simu.zone_size = 220;
    simu.time_step = 1;
    simu.rwp_speed = [0 20] / 3.6;
    simu.rwp_pause_time = [1 5];
    simu.zone_speed = [3 10] / 3.6;
    simu.zone_time = [20 30];
    simu.radio_range = 50;
    transrange = simu.radio_range;
    tweakfloodval = 1.5;
    noofstartevents = 50;
    TTL = 12;
    %disp('----------------STARTING SIMULATIONS------------------')
    Eo = 0.5;
    Eelec=50*10^(-9);
    ETx=50*10^(-9);
    ERx=50*10^(-9); 
    Eamp=100*10^(-12); 
    EDA=5*10^(-9); 
    ks =4000;
    energy =0;
    if (1)
        steps_fig = figure;
        s = get(0,'ScreenSize');
        set(steps_fig,'OuterPosition',[s(3)*1.8/4 s(4)/5 500 500],'NextPlot','add');
        steps_plot1 = plot(0);
        steps_axes = get(steps_fig,'CurrentAxes');
        set(steps_axes,'XLimMode','manual','XLim',[0 simu.grid_size*simu.zone_size],'YLimMode','manual','YLim',[0 simu.grid_size*simu.zone_size],'XTick',0 : simu.zone_size : simu.grid_size*simu.zone_size,'YTick',0 : simu.zone_size : simu.grid_size*simu.zone_size,'XGrid','on','YGrid','on','NextPlot','add');
        steps_plot2 = plot(0);
        steps_plot3 = plot(0);
        set(steps_plot1,'LineStyle','none','Marker','o','MarkerSize',7,'MarkerFaceColor','b');
        set(steps_plot2,'LineStyle','none','Marker','.','MarkerEdgeColor','g','MarkerSize',1);
        set(steps_plot3,'LineStyle','none','Marker','.','MarkerEdgeColor','r','MarkerSize',4.7);
    end

    nodedata = struct('a',[],'b',[],  'xcc', [] , 'ycc', []);
    % creation of nodes
    nodes = Group(simu.alpha,simu.number_of_nodes,simu.grid_size,simu.rwp_speed,simu.rwp_pause_time,simu.zone_speed,simu.zone_time,simu.zone_size,simu.time_step,simu.radio_range);
    %disp(nodes)
    sample(1) = nodes;
    for ii = 1:simu.duration
        % movement of nodes
        nodes.move
        %nodes.move;
        % connections between nodes
        C = nodes.connect;
        sample(ii+1) = nodes;
        %%disp(nodes)
        if (1)
            C(eye(size(C))==1) = 0;
            [v1 v2] = find(C==1);
            links = [v1 v2];
            set(steps_plot1,'YData',nodes.coords(:,2),'XData',nodes.coords(:,1))
            [xc,yc] = circle(nodes.coords(:,1),nodes.coords(:,2),simu.radio_range);
            set(steps_plot2,'XData',xc,'YData',yc)
            nodedata(ii).xcc =  nodes.coords(:,1); %nodes.coords(n(:,1),:);
            nodedata(ii).ycc =  nodes.coords(:,2); %nodes.coords(n(:,2),:);
            if ~isempty(links)
                n = unique(sort(links,2),'rows');
                [xl,yl] = line_plot(nodes.coords(n(:,1),:),nodes.coords(n(:,2),:));
                nodedata(ii).a = nodes.coords(n(:,1),:);
                nodedata(ii).b = nodes.coords(n(:,2),:);
                set(steps_plot3,'XData',xl,'YData',yl);
            end
            %str(ii+1) = nodes;

            drawnow
        end
    end

    for i2 = 1: nn
        NN(i2).id=i2;
        NN(i2).trk = 0;
        NN(i2).E = Eo;
        NN(i2).cond=0;
        NN(i2).count=1;
        NN(i2).vect = struct('Prev_Hash',{},'Timestamp',{},'Data',{},'Hash',{}, 'Sentby', {}, 'TTL',{}, 'Dest', {},'Status', {});
        NN(i2).vect(1).Prev_Hash = 0;
        NN(i2).vect(1).Timestamp = datestr(now,'mm dd, yyyy HH:MM:SS.FFF ');
        NN(i2).vect(1).Data = i2;
        NN(i2).vect(1).Sentby = -1;
        NN(i2).vect(1).TTL = TTL;
        strr = strcat(mat2str(NN(i2).vect(1).Prev_Hash),mat2str(NN(i2).vect(1).Timestamp),mat2str(NN(i2).vect(1).Data),mat2str(NN(i2).vect(1).Sentby),mat2str(NN(i2).vect(1).TTL),mat2str(NN(i2).vect(1).Dest),mat2str(NN(i2).vect(1).Status));
        hasher = System.Security.Cryptography.SHA256Managed;
        hash0 = uint8(hasher.ComputeHash(uint8(strr)));
        NN(i2).vect(1).Hash = dec2hex(hash0);
    end
    prob = zeros(nn,nn);
    k = zeros(nn,nn);
    eventarr = randi(nn,1,noofstartevents);
    traversearr = randi(nn,1,noofstartevents);

    ageconst = 0.98;
    transconst = 0.25;
% ---------------------------------------------------------------------------------------
%
% 
% 
% 
% 
%        IMPLEMENTATION OF BLOCKCHAIN HAS BEEN REMOVED BECAUSE RESEARCH STILL GOING
%          ON. CONTACT AUTHOR IF ADDITIONAL DETAILS ARE NEEDED.
% 
% 
% 
% 
% 
% 
% ----------------------------------------------------------------------------------------
    function [x,y] = circle(a,b,r)
    t = (0:0.01:2*pi)';
    cosinus = r*cos(t);
    sinus = r*sin(t);
    n = length(t);
    m = length(a);
    x = repmat(a,n,1) + repmat(cosinus,m,1);
    y = repmat(b,n,1) + repmat(sinus,m,1);
    end

    function [x,y] = line_plot(a,b)
    t = (0:0.01:1)';
    n = length(t);
    m = size(a,1);
    x = repmat(a(:,1),n,1) + repmat(t,m,1).*repmat(b(:,1)-a(:,1),n,1);
    y = repmat(a(:,2),n,1) + repmat(t,m,1).*repmat(b(:,2)-a(:,2),n,1);
    end
%end
