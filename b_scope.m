clearvars;close all;clc;
path_color = 'C:/Users/windows/Desktop/Research/3. Code/Old source code/color map';
cd(path_color);
load('MyColormaps.mat');
etime = 70;
dx = 2;
for torque_value = [999]
    for bolt_num = [18]
        folder = "E:/1. LASER-ULTRASONIC DATA/0. Official 11 bolt (May 2020)/3. Frame saved/Official data (26.5)/%g/%g/";
        path_fol = sprintf(folder,bolt_num,torque_value);
        cd(path_fol);
        load('frame.mat');
        frame = frame(1:50,:,:);

        for yj=1:size(frame,2)
        xi=1; xf=41; yi=yj;yf=yj;
        xx= xi:xf;
        yy =yi:yf;

        for ff= 10:etime*2  
          %for ff= 1:etime % only etime for downpropagating
          for M=xf-xi+1
              for N=yf-yi+1
                G2(1:M,1:N,ff)= frame(xx,yy,ff);%*(ff^0.5);
              end
          end
        end
        
        time=0.5:0.5:ff*1000000/2000000; 
        G22_1 =reshape(G2,M,ff);
        distance = dx*(0:size(xx,2)-1);
        [T,D] = meshgrid(time,distance);
        
        figure 
        h=surf(T,D,G22_1); 
        shading interp; 
        colormap(mycmap); 
        xlabel('Time(\mus)'); 
        ylabel('Distance(mm)'); 
        axis tight; view(0,90)
        
        end
    end
end
