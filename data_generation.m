%% TRAIN DATASET
clear all;close all;clc;

% SETUP PARAMETERS
damage_sensitive_mode = 0;
% damage_sensitive_mode = 0 => raw data
% damage_sensitive_mode = 1 => processing data

damage_sensitive_type = 0;
% damage_sensitive_mode = 1 => log(RMS-W)
% damage_sensitive_mode = 2 => variance
% damage_sensitive_mode = 3 => std
% damage_sensitive_mode = 4 => skewness
% damage_sensitive_mode = 5 => kurtosis

filter = 0 ;
% filter = 0 => raw
% filter = 1 => remove down
% filter = 2 => remove up
% filter = 3 => wavenumber adaptive image filtering 

xsize = 40;
ysize = 40;

if filter == 1 || filter == 2 || filter == 4 || filter == 3
    start = 32*2 - 32;
else
    start = 32*2;
end
gap = 15*2;
stop = start + gap;

w_para = 2;
simulation_time = 400;
x_interval = 40;
y_interval = 40;
xsize = 40;
ysize = 40;
temp_mean_s = zeros(xsize,ysize,1); % Mean
frame_mean_s = zeros(xsize,ysize,simulation_time*2); % Frame mean
path_fol = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/3. Frame saved/Train/%g/%g';
%%%%%%%%%%%%%%%%%%%%%%
k_num = 20;
%%%%%%%%%%%%%%%%%%%%%%
for bolt_num = 31
    save_full =[];
    save_slice_final = [];
    save_slice_label = [];
    fprintf('====== Processing Bolt %d ======\n',bolt_num)
    for torque_value = [5 10 15 20 25]
        fprintf('With torque value %d\n',torque_value)
        path_fol_s = sprintf(path_fol,bolt_num,torque_value);
%         disp(path_fol_s);
        cd(path_fol_s);
        % LOAD FRAME
        load('frame.mat');
        
        numsample = length(frame);
        
        if filter == 1
            ft=fftn(frame);
            part1=ft(:,:,1:numsample/2);
            part1(1:xsize/2+0.5,:,:)=zeros(xsize/2,ysize,numsample/2);
            re_data=real(ifftn(part1));
            remove_down =  real(ifftn(part1));
            frame = remove_down;
            disp('Change frame to filter mode');
                 
        elseif filter == 2
            ft=fftn(frame);
            part1=ft(:,:,1:numsample/2);
            part1(1:xsize/2+0.5,:,:)=zeros(xsize/2,ysize,numsample/2);
            re_data=real(ifftn(part1));
            remove_down =  real(ifftn(part1));;
            part2=ft(:,:,1:numsample/2);
            clear ft;
            part2(xsize/2+0.5:xsize,:,:)=zeros(xsize/2,ysize,numsample/2);
            remove_up = real(ifftn(part2));
            frame = remove_up;
            disp('Change frame to filter mode');
            
        elseif filter == 3
            load('frame.mat');
            G =fftn(frame);

            Ghalf=G(:,:,1:400);
            dx = 2;
            dy = dx;

            fs= 2e+6;       % sampling frequency
            fr= fs*(0:size(G,3)/2-1)/(size(G,3)/2);        % frequency
            Kx= 2*pi*(-size(G,2)/2:size(G,2)/2-1)/(size(G,2)*dx); % wavenumber in y direction
            Ky= 2*pi*(-size(G,1)/2:size(G,1)/2-1)/(size(G,1)*dy); % wavenumber in x direction
            % CACULATE THE AVERAGED WAVEFIELD IN WAVENUMBER-WAVENUMBER DOMAINE
            MagStack = [];
            MagSum = [];
            k = 1;
            for interest_time = 1:stop
                Magntd=abs(Ghalf(:,:,interest_time)); % the spectrum at the frequency of interest
                MagStack(:,:,k) = Magntd;
                k = k +1;
            end

            MagSum = zeros(xsize,ysize);
            for sum_time = 1:stop-start+1
                MagSum = MagSum + MagStack(:,:,sum_time);
            end
            MagAvg = MagSum./(stop-1);
            %%%%%%%%%%%%%%%%%%%%%%%

            Magntd_flatten = [];
            Magntd_flatten = reshape(MagAvg.',1,[]);
            Magntd_sort = sort(Magntd_flatten);

            num_remove = round(2*length(Magntd_sort)/100,0);
            threshold = Magntd_sort(length(Magntd_sort)-num_remove);

            [xf, yf] = find(MagAvg>threshold);

            M_filter = ones(xsize,ysize);
            for j = 1:length(xf)
                M_filter(xf(j),yf(j)) = 0;
            end
            
            for time = 1:length(Ghalf)
                frame(:,:,time) = real(ifftn(M_filter.*Ghalf(:,:,time)));
            end
            disp('Change frame to filter mode');
        elseif filter ==4
            ft=fftn(frame);
            part1=ft(:,:,1:numsample/2);
            part1(1:xsize/2+0.5,:,:)=zeros(xsize/2,ysize,numsample/2);
            remove_down =  real(ifftn(part1));
            frame = remove_down;
            disp('Change frame to filter mode');

            part2=ft(:,:,1:numsample/2);
            part2(xsize/2+0.5:xsize,:,:)=zeros(xsize/2,ysize,numsample/2);
            
            remove_up = real(ifftn(part2));
            
            temp_swe = zeros(xsize,ysize,1);;
            frame_swe = zeros(xsize,ysize,simulation_time*2);
            for f = 1:simulation_time
                temp_swe = temp_swe + frame(:,:,f).^2 - ...
                    remove_down(:,:,f).^2 - remove_up(:,:,f).^2;
                swe = (temp_swe);
                frame(:,:,f) = swe;
            end
            
        end 
        
        if damage_sensitive_mode == 1
            % Mean
            for f = 1:simulation_time
                temp_mean_s = temp_mean_s + frame(:,:,f);
                mean_s = temp_mean_s / f;
                frame_mean_s(:,:,f) = mean_s;
            end
            % Variance
            temp_variance_s = zeros(xsize,ysize,1);
            frame_variance_s = zeros(xsize,ysize,simulation_time*2);
            for f = 1:simulation_time
                temp_variance_s = temp_variance_s + (frame(:,:,f) - mean_s).^2;
                variance_s = temp_variance_s / f;
                frame_variance_s(:,:,f) = variance_s;
            end
            % Standarad deviation 
            temp_std_s = zeros(xsize,ysize,1);
            frame_std_s = zeros(xsize,ysize,simulation_time*2);
            for f = 1:simulation_time
                temp_std_s = temp_std_s + (frame(:,:,f) - mean_s).^2;
                std_s = sqrt(temp_std_s / f);
                frame_std_s(:,:,f) = std_s;
            end
            % Skewness
            temp_skewness_s = zeros(xsize,ysize,1);
            frame_skewness_s = zeros(xsize,ysize,simulation_time*2);
            for f = 1:simulation_time
                temp_skewness_s = temp_skewness_s + (frame(:,:,f) - mean_s).^3;
                skewness_s = (temp_skewness_s/f)./(std_s.^3);
                frame_skewness_s(:,:,f) = skewness_s;
            end
            % Kurtosis
            temp_kurtosis_s = zeros(xsize,ysize,1);
            frame_kurtosis_s = zeros(xsize,ysize,simulation_time*2);
            for f = 1:simulation_time
                temp_kurtosis_s = temp_kurtosis_s + (frame(:,:,f) - mean_s).^4;
                kurtosis_s = (temp_kurtosis_s/f)./(std_s.^4);
                frame_kurtosis_s(:,:,f) = kurtosis_s;
            end
            % RMS
            temp_rms_s = zeros(xsize,ysize,1);
            frame_rms_s = zeros(xsize,ysize,simulation_time*2);
            for f = 1:simulation_time
                temp_rms_s = temp_rms_s + frame(:,:,f).^2;
                rms_s = sqrt(temp_rms_s/f);
                frame_rms_s(:,:,f) = rms_s;
            end
            % RMS weight
            temp_rms_w_s = zeros(xsize,ysize,1);
            frame_rms_w_s = zeros(xsize,ysize,simulation_time*2);
            for f = 1:simulation_time
                temp_rms_w_s = temp_rms_w_s + (frame(:,:,f).^2)*(f^w_para);
                rms_w_s = sqrt(temp_rms_w_s/f);
                frame_rms_w_s(:,:,f) = rms_w_s;
            end 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Change frame to damage sensitive mode');
            frame = log(frame_rms_w_s);  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        % STACK
        save_slice = [];
        for time = start:stop
            slice = frame(:,:,time);
            row = reshape(slice.',1,[]);
            save_slice = [save_slice; row];
        end
        
        label = torque_value*ones(stop-start+1,1);
        
        save_slice_label = [save_slice, label];
        save_slice_final = [save_slice_final; save_slice_label];
    end
    save_full = [save_full; save_slice_final];
    
    % SAVE .CSV
    if filter == 0 && damage_sensitive_mode == 0
        path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/1.Raw/Train';
        cd(path_save);
        file_name = sprintf('Bolt%d.csv',bolt_num);
        csvwrite(file_name,save_full)
    elseif filter == 1 && damage_sensitive_mode == 0
        path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/2.Filter1/Train';
        cd(path_save);
        file_name = sprintf('Bolt%d.csv',bolt_num);
        csvwrite(file_name,save_full)       
    elseif filter == 2 && damage_sensitive_mode == 0
        path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/3.Filter2/Train';
        cd(path_save);
        file_name = sprintf('Bolt%d.csv',bolt_num);
        csvwrite(file_name,save_full)
    elseif filter == 3 && damage_sensitive_mode == 0
        path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/4.Adaptive/Train';
        cd(path_save);
        file_name = sprintf('Bolt%d.csv',bolt_num);
        csvwrite(file_name,save_full)
    elseif filter == 4 && damage_sensitive_mode == 0
        path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/9.StandingWaves/Train';
        cd(path_save);
        file_name = sprintf('Bolt%d.csv',bolt_num);
        csvwrite(file_name,save_full)
        
    elseif damage_sensitive_type == 1 && filter == 0
        path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/5.RawLogRMSW/Train';
        cd(path_save);
        file_name = sprintf('Bolt%d.csv',bolt_num);
        csvwrite(file_name,save_full)
    elseif damage_sensitive_type == 1 && filter == 1
        path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/6.Filter1LogRMSW/Train';
        cd(path_save);
        file_name = sprintf('Bolt%d.csv',bolt_num);
        csvwrite(file_name,save_full)
    elseif damage_sensitive_type == 1 && filter == 2
        path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/7.Filter2LogRMSW/Train';
        cd(path_save);
        file_name = sprintf('Bolt%d.csv',bolt_num);
        csvwrite(file_name,save_full) 
    elseif damage_sensitive_type == 1 && filter == 3 
        path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/8.AdaptiveLogRMSW/Train';
        cd(path_save);
        file_name = sprintf('Bolt%d.csv',bolt_num);
        csvwrite(file_name,save_full) 
        elseif damage_sensitive_type == 1 && filter == 4
        path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/10.StandingWavesLogRMSW/Train';
        cd(path_save);
        file_name = sprintf('Bolt%d.csv',bolt_num);
        csvwrite(file_name,save_full)  
    end
end
disp('Done');
% %% MERGE TRAIN - DEV - TEST SET
if filter == 0 && damage_sensitive_mode == 0
    path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/1.Raw/Train';
    cd(path_save);

elseif filter == 1 && damage_sensitive_mode == 0
    path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/2.Filter1/Train';
    cd(path_save);
    
elseif filter == 2 && damage_sensitive_mode == 0
    path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/3.Filter2/Train';
    cd(path_save);

elseif filter == 3 && damage_sensitive_mode == 0
    path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/4.Adaptive/Train';
    cd(path_save);
elseif filter == 4 && damage_sensitive_mode == 0
    path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/9.StandingWaves/Train';
    cd(path_save);
elseif damage_sensitive_type == 1 && filter == 0
    path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/5.RawLogRMSW/Train';
    cd(path_save);

elseif damage_sensitive_type == 1 && filter == 1
    path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/6.Filter1LogRMSW/Train';
    cd(path_save);

elseif damage_sensitive_type == 1 && filter == 2
    path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/7.Filter2LogRMSW/Train';
    cd(path_save);

elseif damage_sensitive_type == 1 && filter == 3 
    path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/8.AdaptiveLogRMSW/Train';
    cd(path_save);
elseif damage_sensitive_type == 1 && filter == 4
    path_save = 'E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/4. CSV saved/2. Regression/10.StandingWavesLogRMSW/Train';
    cd(path_save);
end

train = [];
dev = [];
test = [];
for kfol = 1:k_num
    train = [];
    dev = [];
    fprintf('====== Processing kfold %d ======\n',kfol)
    for bolt_num = 1:k_num
       fprintf('====== Processing Bolt %d ======\n',bolt_num)
       file_name = sprintf('Bolt%d.csv',bolt_num);
       data_name = sprintf('Bolt%d',bolt_num);
       load(file_name);
       if  bolt_num == kfol
           dev = [dev; eval(data_name)];
           fprintf('Bolt %d in the DEV set\n',bolt_num);
       else 
           train = [train; eval(data_name)]; 
           fprintf('Bolt %d in the TRAIN set\n',bolt_num);
       end
    end
    disp(size(train))
    disp(size(dev))
    csvwrite(sprintf('train%d.csv',kfol), train);
    csvwrite(sprintf('dev%d.csv',kfol), dev);
    disp('Done');
end
disp('Done');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%