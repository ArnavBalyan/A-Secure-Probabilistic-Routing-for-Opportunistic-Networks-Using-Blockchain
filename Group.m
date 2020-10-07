% Inputs (in appearance order) : 
% - alpha
% - node_config
% - grid_size
% - v_rwp
% - t_rwp
% - v_zone
% - t_zone
% - zone_size
% - time_step
% - radio_range
classdef Group < handle
    
    properties
        v_rwp 
        t_rwp
        v_zone
        t_zone
        time_step
        zone_size
        pause_time
        speed
        next_waypoint
        coords
        number_of_nodes
        grid_size
        grid_coords      
        zone_time
        alpha
        attachment_zone
        d_cdf
        d_max
        radio_range
        roaming
    end
    
    methods
        
        function self = Group(alpha, node_config, grid_size, v_rwp, t_rwp, v_zone, t_zone, zone_size, time_step, radio_range)
            self.v_rwp = v_rwp;
            self.t_rwp = t_rwp;
            self.v_zone = v_zone;
            self.t_zone = t_zone;
            self.zone_size = zone_size;
            self.time_step = time_step;
            self.number_of_nodes = sum(node_config); 
            self.grid_size = grid_size;
            self.alpha = alpha;
            self.radio_range = radio_range;
          
            % initialisation
%             grid_init = 1:2:self.grid_size;
%             self.grid_coords = [randsample(grid_init,self.number_of_nodes,1)' randsample(grid_init,self.number_of_nodes,1)']; % equi-distanced zones
            
self.grid_coords = randi(self.grid_size,self.number_of_nodes,2)-1; % uniformly
            
            self.attachment_zone = self.grid_coords;
            self.coords = self.grid_coords*self.zone_size + self.zone_size * rand(self.number_of_nodes,2);
            self.next_waypoint = self.grid_coords*self.zone_size + self.zone_size * rand(self.number_of_nodes,2);
            self.speed = (self.v_rwp(2) - self.v_rwp(1)) * rand(self.number_of_nodes,1) + self.v_rwp(1);
            self.pause_time = (self.t_rwp(2) - self.t_rwp(1)) * rand(self.number_of_nodes,1) + self.t_rwp(1);
            self.zone_time = (self.t_zone(2) - self.t_zone(1))*rand(self.number_of_nodes,1) + self.t_zone(1);
            self.roaming = zeros(self.number_of_nodes,1);
            
            % Matrix d_cdf (culmulative probability distribution of distance) serves for the calculation of new zone
            self.d_max = floor(grid_size/2);
            for i=1:length(node_config)
                beta = sum((1 + (0:self.d_max)).^(-alpha(i)))^-1;
                aux = beta*((0:self.d_max)' + 1).^(-alpha(i));
                C = sum(triu(repmat(aux,1,self.d_max+1)));
                self.d_cdf = cat(1,self.d_cdf,repmat(C,node_config(i),1));
            end
      
        end
        
        %------------------------------------------------------------------
        % Function making the movement of all the nodes in one time step
        %------------------------------------------------------------------      
        function move(self)
            % nodes staying inside zones
            in_zone = self.zone_time ~= 0;
   
            % % moving
            in_moving = in_zone & self.pause_time == 0;
            
            d_to_next_wp = ones(self.number_of_nodes,1)*-1;
            d_to_next_wp(in_moving) = sqrt(sum((self.next_waypoint(in_moving,:)-self.coords(in_moving,:)).^2,2));
            
            
            % % % reached the Waypoint
            in_arrived = in_moving & d_to_next_wp == 0;
            
            if ~isequal(in_arrived,zeros(self.number_of_nodes,1))
                new_wp = self.zone_size*rand(self.number_of_nodes,2);
                self.next_waypoint(in_arrived,:) = self.grid_coords(in_arrived,:) * self.zone_size + new_wp(in_arrived,:);                
                new_pause = (self.t_rwp(2) - self.t_rwp(1)) * rand(self.number_of_nodes,2) + self.t_rwp(1);
                self.pause_time(in_arrived) = new_pause(in_arrived);              
                new_speed = (self.v_rwp(2) - self.v_rwp(1)) * rand(self.number_of_nodes,1) + self.v_rwp(1);
                self.speed(in_arrived) = new_speed(in_arrived);
            end
            
            % % % on the way to the next waypoint
            in_way = in_moving & d_to_next_wp ~= 0;
            d_move = ones(self.number_of_nodes,1)*-1;
            d_move(in_way) = self.time_step * self.speed(in_way);
         
            self.coords(in_way & d_move >= d_to_next_wp,:) = self.next_waypoint(in_way & d_move >= d_to_next_wp,:);
            if ~isequal(in_way & d_move < d_to_next_wp, zeros(self.number_of_nodes,1))
                aux = self.next_waypoint(in_way & d_move < d_to_next_wp,:) - self.coords(in_way & d_move < d_to_next_wp,:);
                aux1 = [aux(:,1).*self.time_step.*self.speed(in_way & d_move < d_to_next_wp)./d_to_next_wp(in_way & d_move < d_to_next_wp) aux(:,2).*self.time_step.*self.speed(in_way & d_move < d_to_next_wp)./d_to_next_wp(in_way & d_move < d_to_next_wp)];
                self.coords(in_way & d_move < d_to_next_wp,:) = self.coords(in_way & d_move < d_to_next_wp,:) + aux1;
            end
            
            % % taking a pause
            pausing = in_zone & self.pause_time ~= 0 & ~d_to_next_wp==0; % node o trong zone, dang pause nhung ko phai la nhung node vua moi den noi
            self.pause_time(pausing & self.pause_time <= self.time_step) = 0;
            self.pause_time(pausing & self.pause_time > self.time_step) = self.pause_time(pausing & self.pause_time > self.time_step) - self.time_step;
        
            % update residual time of the stay
            if ~isequal(in_zone & self.zone_time <= self.time_step,zeros(self.number_of_nodes,1))
                self.zone_time(in_zone & self.zone_time <= self.time_step) = 0;
                new_zone = self.power_law(self.attachment_zone,self.grid_size,self.number_of_nodes);
                new_wp = ones(self.number_of_nodes,2) * self.zone_size .* rand(self.number_of_nodes,2);
                new_speed = (self.v_zone(2) - self.v_zone(1)) * rand(self.number_of_nodes,1) + self.v_zone(1);
                self.next_waypoint(in_zone & self.zone_time <= self.time_step,:) = new_zone(in_zone & self.zone_time <= self.time_step,:) * self.zone_size + new_wp(in_zone & self.zone_time <= self.time_step,:);
                self.grid_coords(in_zone & self.zone_time <= self.time_step,:) = new_zone(in_zone & self.zone_time <= self.time_step,:);
                self.speed(in_zone & self.zone_time <= self.time_step) = new_speed(in_zone & self.zone_time <= self.time_step);
            end
            self.zone_time(in_zone & self.zone_time > self.time_step) = self.zone_time(in_zone & self.zone_time > self.time_step) - self.time_step;
           
            
            % nodes moving to a new zone
            out_zone = self.zone_time == 0 & ~in_zone; % node het thoi gian o trong zone, nhung ko phai la nhung node vua moi het thoi gian duoc cap nhat o tren
            if ~isequal(out_zone,zeros(self.number_of_nodes,1))
                [direction d_to_next_zone] = self.shortest_distance_torus(self.coords, self.next_waypoint, out_zone);
           
                % % arrived at the new zone
                out_arrived = out_zone & d_to_next_zone == 0;
                if ~isequal(out_arrived,zeros(self.number_of_nodes,1))
                    new_zone_time = (self.t_zone(2) - self.t_zone(1))*rand(self.number_of_nodes,1) + self.t_zone(1);
                    self.zone_time(out_arrived) = new_zone_time(out_arrived);
                    new_wp = ones(self.number_of_nodes,2) * self.zone_size .* rand(self.number_of_nodes,2);
                    self.next_waypoint(out_arrived,:) = self.grid_coords(out_arrived,:)*self.zone_size + new_wp(out_arrived,:);
                    new_speed = (self.v_rwp(2) - self.v_rwp(1)) * rand(self.number_of_nodes,1) + self.v_rwp(1);
                    self.speed(out_arrived) = new_speed(out_arrived);
                    self.roaming(out_arrived) = 0;
                end
                
                % % on the way to the new zone
                out_way = out_zone & d_to_next_zone ~= 0;
                if ~isequal(out_way,zeros(self.number_of_nodes,1))
                    d_move = ones(self.number_of_nodes,1)*-1;
                    d_move(out_way) = self.time_step * self.speed(out_way);
                    
                    self.coords(out_way & d_move >= d_to_next_zone,:) = self.next_waypoint(out_way & d_move >= d_to_next_zone,:);
                    
                    out_way_muti_steps = out_way & d_move < d_to_next_zone;
                    if ~isequal(out_way_muti_steps,zeros(self.number_of_nodes,1))
                        aux = [d_move(out_way_muti_steps) .* direction(out_way_muti_steps,1) d_move(out_way_muti_steps) .* direction(out_way_muti_steps,2)];
                        self.coords(out_way_muti_steps,:) = self.coords(out_way_muti_steps,:) + aux;
                        aux1 = self.coords(out_way_muti_steps,:);
                        aux1(aux1 < 0) = aux1(aux1 < 0) + self.zone_size*self.grid_size;
                        aux1(aux1 >= self.zone_size*self.grid_size) = aux1(aux1 >= self.zone_size*self.grid_size) - self.zone_size*self.grid_size;
                        self.coords(out_way_muti_steps,:) = aux1;
                    end
                    
                    self.roaming(out_way) = 1;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Function finding the distance between all the nodes and check if
        % there exists connections between them (i.e. distance < radio
        % range)
        %  - connection_matrix : the network's adjacency matrix
        %------------------------------------------------------------------        
        function connection_matrix = connect(self)
            [X1,X2]=meshgrid(self.coords(:,1));
            [Y1,Y2]=meshgrid(self.coords(:,2));
            
            dX1 = X2 - X1;
            dY1 = Y2 - Y1;
            
            % For translated coordinates
            dX2 = abs(dX1) - self.grid_size*self.zone_size;
            dY2 = abs(dY1) - self.grid_size*self.zone_size;
            
            D1 = sqrt(dX1.^2 + dY1.^2);
            D2 = sqrt(dX2.^2 + dY2.^2);
            D3 = sqrt(dX2.^2 + dY1.^2);
            D4 = sqrt(dX1.^2 + dY2.^2);

            D = min(D1, min(D2, min(D3,D4)));
            
            connection_matrix = eye(self.number_of_nodes);
            connection_matrix(D <= self.radio_range) = 1;
          
            % excluding contacts while node is roaming
%             mask_roaming = or(repmat(self.roaming,1,self.number_of_nodes),repmat(self.roaming',self.number_of_nodes,1));
%             connection_matrix(mask_roaming & eye(size(connection_matrix))==0) = 0;
        end
        
        %------------------------------------------------------------------
        % Function generating random zone coords following the power law 
        %                 P(x) = beta*x^(-alpha)
        %------------------------------------------------------------------
        function new_zone = power_law(self, attachment_zone, grid_size, number_of_nodes)
            R = rand(number_of_nodes,1);
            RR = repmat(R,1,self.d_max+1);
            
            distance = self.d_max + 1 - sum(RR <= self.d_cdf,2);
            r = [randi(4,number_of_nodes,1);-1]; % padding -1 for the case where there is 1 node, the right hand element of equation D(r==x,:) has size (0,2) when distance(r==x) is empty
            r(distance == 0) = 0; % the following process is wrong for d=0
            
            D = zeros(number_of_nodes,2);
            D(r==1,:) = [-distance(r==1).*ones(size(distance(r==1))) floor(rand(size(distance(r==1)))*2.*distance(r==1))-distance(r==1)];
            D(r==2,:) = [floor(rand(size(distance(r==2)))*2.*distance(r==2))-distance(r==2) distance(r==2).*ones(size(distance(r==2)))];
            D(r==3,:) = [distance(r==3).*ones(size(distance(r==3))) floor(rand(size(distance(r==3)))*2.*distance(r==3))-distance(r==3)+1];
            D(r==4,:) = [floor(rand(size(distance(r==4)))*2.*distance(r==4))-distance(r==4)+1 -distance(r==4).*ones(size(distance(r==4)))];
           
%             D = floor(rand(number_of_nodes,2).*repmat(2*distance+1,1,2)) - repmat(distance,1,2);
            new_zone = attachment_zone + D;
            new_zone(new_zone < 0) = new_zone(new_zone < 0) + grid_size;
            new_zone(new_zone >= grid_size) = new_zone(new_zone >= grid_size) - grid_size;
        end
        
        %------------------------------------------------------------------
        % Function returning the shortest distance between 2 points X, Y in
        % a torus and the angle created by the line XY and axe Ox. 
        %  - X,Y : matrix containing the coordinates of points 
        %  - Z : a filter, only X(Z,:) and Y(Z,:) are useful
        %------------------------------------------------------------------
        function [direction, distance] = shortest_distance_torus(self, X, Y, Z)
            dX = Y(:,1) - X(:,1);
            dY = Y(:,2) - X(:,2);
            
            field_size = self.grid_size * self.zone_size;
           
            C1 = [X(:,1) + dX X(:,2) + dY];
            C2 = [X(:,1) + dX X(:,2) + (-dY./abs(dY)).*(field_size - abs(dY))];
            C3 = [X(:,1) + (-dX./abs(dX)).*(field_size - abs(dX)) X(:,2) + (-dY./abs(dY)).*(field_size - abs(dY))];
            C4 = [X(:,1) + (-dX./abs(dX)).*(field_size - abs(dX)) X(:,2) + dY];
            
            C = zeros(size(X));
            C(abs(dX) < field_size - abs(dX) & abs(dY) < field_size - abs(dY),:) = C1(abs(dX) < field_size - abs(dX) & abs(dY) < field_size - abs(dY),:);
            C(abs(dX) < field_size - abs(dX) & abs(dY) >= field_size - abs(dY),:) = C2(abs(dX) < field_size - abs(dX) & abs(dY) >= field_size - abs(dY),:);
            C(abs(dX) >= field_size - abs(dX) & abs(dY) >= field_size - abs(dY),:) = C3(abs(dX) >= field_size - abs(dX) & abs(dY) >= field_size - abs(dY),:);
            C(abs(dX) >= field_size - abs(dX) & abs(dY) < field_size - abs(dY),:) = C4(abs(dX) >= field_size - abs(dX) & abs(dY) < field_size - abs(dY),:);
            
            distance = ones(size(X,1),1)*-1;
            direction = ones(size(X))*-1;
            
            distance(Z) = sqrt(sum((C(Z,:)-X(Z,:)).^2,2));
            direction(Z,:) = [(C(Z,1) - X(Z,1))./distance(Z) (C(Z,2) - X(Z,2))./distance(Z)];
        end
        
    end
    
    methods (Static)
        function d = testPowerLaw
            number_of_nodes = 1000;
            v=[3 5]./3.6;
            t_rwp=[0 0];
            t_zone = [10 10];
            grid_size = 10000;
            square_size=10;
            time_step = 5;
            alpha = 1;
       
            group = RWP_group(alpha,number_of_nodes,grid_size,v,t_rwp,t_zone,square_size,time_step);
            
            [Z d] = group.power_law(zeros(number_of_nodes,2),grid_size,number_of_nodes);
            figure(2)
            hist(d,max(d)+1)
%             [f,x] = ecdf(d);
%             
%             f = 1 - f;
%             figure(1)
%             loglog(x,f,'.-')
        end
        
        function coords = checkDistribution
            clc
            close all
            % random picking
            N = 10000;
            d = 2;
            tic
            r = randi(4,N,1);
            coords = zeros(N,2);
            coords(r==1,:) = [-d*ones(size(coords(r==1,:),1),1) randsample(-d:d-1,size(coords(r==1,:),1),1)'];
            coords(r==2,:) = [randsample(-d:d-1,size(coords(r==2,:),1),1)' d*ones(size(coords(r==2,:),1),1)];
            coords(r==3,:) = [d*ones(size(coords(r==3,:),1),1) randsample(-d+1:d,size(coords(r==3,:),1),1)'];
            coords(r==4,:) = [randsample(-d+1:d,size(coords(r==4,:),1),1)' -d*ones(size(coords(r==4,:),1),1)];
            toc
            % verify the distribution
            D(:,1) = [-d:d-1 ones(1,2*d)*d -d+1:d ones(1,2*d)*-d]';
            D(:,2) = [ones(1,2*d)*d -d+1:d ones(1,2*d)*-d -d:d-1]';
            
            DD = zeros(N,1);
            for i=1:N
                DD(i) = find(D(:,1)==coords(i,1)&D(:,2)==coords(i,2));
            end
      
            hist(DD,8*d);          
        end
    end
    
end
