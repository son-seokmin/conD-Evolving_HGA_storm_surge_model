clc;    clear;  close all;
%% Model selection
name_model = 'Houston_AL_canesm_20thcal';       % CanESM, present
% name_model = 'Houston_AL_canesm_ssp585cal';     % CanESM, future
% name_model = 'Houston_AL_gfdl6_20thcal';        % GFDL, present
% name_model = 'Houston_AL_gfdl6_ssp585cal';      % GFDL, future
% name_model = 'Houston_AL_hadgem6_20thcal';      % HadGEM, present
% name_model = 'Houston_AL_hadgem6_ssp585cal';    % HadGEM, future

path_model = [pwd,'\',name_model];
load([path_model,'\',name_model,'.mat']);       % load dataset
load([path_model,'\vnet.mat']);                 % load velocity data

%% TC filtering
lim1 = 96*0.44704;
lim2 = 111*0.44704;
lim3 = 130*0.44704;
lim4 = 157*0.44704;
lon_limit = -93.87;
lat_limit = 29.55;
cat_limit = 27.5;
vnet_max = zeros(4500,1);
for j=1:4500
    k = find(latstore(j,:)>lat_limit,1);
    if isempty(k)==0
        if longstore(j,k)-360<lon_limit
        	k_lim = find(latstore(j,:)>cat_limit);
            vnet_max(j) = max(vnet(j,k_lim))*0.51444;
        end
    end
end
track_list = transpose(find(vnet_max>=lim2));

%% TC modeling
for track_num = track_list(1:end)
    spider_name = [path_model,'\Wind_grid_',num2str(track_num),'.spw'];
    idx_start = find(latstore(track_num,:)>=22,1);      % start time point
    idx_final = find(monthstore(track_num,:)==0,1);     % end time point
    
    year = double(yearstore(1,track_num));
    month = monthstore(track_num,idx_start:idx_final-1);
    day = daystore(track_num,idx_start:idx_final-1);
    hour = hourstore(track_num,idx_start:idx_final-1);
    time = datenum(year,month,day,hour,0,0)*86400;
    lat = latstore(track_num,idx_start:idx_final-1);        % latitude
    lon = longstore(track_num,idx_start:idx_final-1)-360;   % longitude
    p_c = pstore(track_num,idx_start:idx_final-1)*100;      % central pressure
    V_max = vnet(track_num,idx_start:idx_final-1)/1.944;    % maximum wind velocity
    R_max = rmstore(track_num,idx_start:idx_final-1)*1000;  % radius of maximum wind
    omega = 2*pi/24/3600;                                   % Angular velocity of the Earth
    f = 2*omega*sind(lat);                                  % Coriolis parameter
    
    % Constant parameters
    R_earth = 6371000;
    rho = 1.2;
    p_n = 101500;
    
    % Grid in polar coordinates
    r = 1000:1000:3000000;
    angle = mod(-(0:10:359)+180,360);
    
    % Velocity coefficients
    K = 1.08/1.32;      % Conversion factor (Gust factor)
    V_max = K*V_max;
    Ct = 0.5;           % Coefficient in Xie et al. (2006)
    
    %% Time marching
    for t = 2:length(time)
        % TC translation
        dx = (lon(t)-lon(t-1))*R_earth*cosd(lat(t-1))*pi/180;   % translation in x
        dy = (lat(t)-lat(t-1))*R_earth*pi/180;                  % translation in y
        dt = time(t)-time(t-1);                                 % time-step
        Vt = sqrt(dx^2 + dy^2) / dt;                            % translation velocity
        alpha = atan2(dy,dx)*180/pi;                            % translation direction
        Vt_x = Vt*cosd(alpha);                                  % x-component
        Vt_y = Vt*sind(alpha);                                  % y-component
        
        % Holland parameter (Holland et al. 2010)
        B = rho*exp(1)*V_max(t)^2/(p_n - p_c(t));
        A = R_max(t)^B;
        
        for i=1:length(r)
            % Graham and Nunn (1959)
            if r(i)<=R_max(t)
                beta = 10*r(i)/R_max(t);
            elseif (r(i)>R_max(t) && r(i)<=1.2*R_max(t))
                beta = 10 + 75*(r(i)/R_max(t)-1);
            else
                beta = 25;
            end
            
            for j=1:length(angle)
                % pressure field
                p(i,j,t) = p_c(t) + (p_n - p_c(t))*exp(-A/r(i)^B);    
                % wind velocity field
                Vc = sqrt(A*B*(p_n - p_c(t))*exp(-A/r(i)^B)/rho/r(i)^B);    % cyclostrophic wind speed
                Vg = sqrt(Vc^2 + (f(t)*r(i))^2/4) - f(t)*r(i)/2;            % geostrophic wind speed
                Vg_x = Vg*cosd(angle(j) + beta);                            % x-component
                Vg_y = Vg*sind(angle(j) + beta);                            % y-component
                % Xie et al. (2006)
                Vx = Vg_x + Ct*Vt_x;                                        
                Vy = Vg_y + Ct*Vt_y;
                V(i,j,t) = sqrt(Vx^2 + Vy^2);
                theta(i,j,t) = mod(270-atan2(Vy,Vx)*180/pi,360);            % angle in nautical convention
            end
        end    
    end
    p(find(isnan(p))) = p_n;
    V(find(isnan(V))) = 0;
    
    % Write spw file
    h1 = fopen(spider_name,'w');
    fprintf(h1,'FileVersion = 1.03\n');
    fprintf(h1,'filetype = meteo_on_spiderweb_grid\n');
    fprintf(h1,'NODATA_value = -999.000\n');
    fprintf(h1,'n_cols = %4i\n',length(angle));
    fprintf(h1,'n_rows = %4i\n',length(r));
    fprintf(h1,'grid_unit = degree\n');
    fprintf(h1,'spw_radius = %8.3f\n',max(r));
    fprintf(h1,'spw_rad_unit = m\n');
    fprintf(h1,'n_quantity = 3\n');
    fprintf(h1,'quantity1 = wind_speed\n');
    fprintf(h1,'quantity2 = wind_from_direction\n');
    fprintf(h1,'quantity3 = p_drop\n');
    fprintf(h1,'unit1 = m s-1\n');
    fprintf(h1,'unit2 = degree\n');
    fprintf(h1,'unit3 = Pa\n');
    fclose(h1);
   
    p_drop_c = p_n - p_c;   % pressure drop (central p)
    p_drop = p_n - p;       % pressure drop (p field)
    t_interval = (time - time(2))/3600;
    for t=2:length(time)
        h1 = fopen(spider_name,'a');
        fprintf(h1,'TIME = %5.1f hours since %4i-%02i-%02i %02i:%02i:00 +00:00\n',t_interval(t),year,month(2),day(2),hour(2),0);
        fprintf(h1,'x_spw_eye = %8.3f\n',lon(t));
        fprintf(h1,'y_spw_eye = %8.3f\n',lat(t));
        fprintf(h1,'p_drop_spw_eye = %8.3f\n',p_drop_c(t));
        fclose(h1);
        allwrite=[V(:,:,t); theta(:,:,t); p_drop(:,:,t)];
        dlmwrite(spider_name,allwrite,'delimiter',' ','precision',7,'-append');
    end 
end
