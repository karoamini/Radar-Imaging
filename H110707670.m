clc 
clear all
%% configuration

c = 3*10^8; % light speed
f = 9.6*10^9; % carrier frequency
lambda = c / f;
number_of_tx = 3; 
number_of_rx = 17;
max_theta = atan(0.35);
dx = lambda/(2*sin(max_theta)); 
dx_rx = dx ; % distance between tow consecutive receivers
dx_tx =number_of_rx*dx_rx ; % distance between two consecutive transmitters

%% part 3.1.1 : place the target in a random position
xBox= [-35,-35,35,35,-35];
yBox= [100,150,150,100,100];
figure(1)
plot(xBox,yBox)
hold on

target_x1 = -35+70*rand; %In general, you can generate N random numbers in the interval (a,b) with the formula r = a + (b-a).*rand(N,1).
target_y1 = 100+50*rand;
plot (target_x1, target_y1, '^' )
xlim ([-80 80])
ylim ([-10 160])



%% part 3.1.2 placing TX and RX elements in 2D space


mode_tx = mod(number_of_tx,2);
mode_rx = mod(number_of_rx,2);

if mode_rx == 0
    xn = dx_rx*((number_of_rx+1)/2)
    rx_position = (-xn+dx_rx : dx_rx : xn-dx_rx );
     
else
    xn = ((number_of_rx-1)/2)*dx_rx;  
    rx_position = (-xn : dx_rx : xn);
   % rx_position = rx_position + dx_rx/2;
    
    
end


if mode_tx == 0
    xn = dx_tx*((number_of_tx+1)/2);
    tx_position = (-xn+dx_tx : dx_tx : xn-dx_tx );
     
else
   xn = ((number_of_tx-1)/2)*dx_tx;  
    tx_position = (-xn : dx_tx : xn);

    
    
end

tx_position;
rx_position;

 plot(tx_position,0,'*')
 plot(rx_position,0,'.')
hold off



%% part 3.1.3 : theta1 calculation

theta1 =  rad2deg(atan(target_x1 / target_y1))

%% 3.1.4 : propagating the signal

R0_1 = zeros(1,number_of_tx);
R1_1= zeros(1,number_of_rx);
for m = 1:1:number_of_tx
   R0_1(1,m) = sqrt((target_x1 -(tx_position(1,m)))^2 + target_y1^2); % calculating R1 and store it in an array
end

for n = 1:1:number_of_rx
     R1_1(1,n) = sqrt((target_x1 -rx_position(1,n))^2 + target_y1^2);% calculating R2 and store it in an array
end
w=rand;
phi=rand*2*pi;
row_p = w*exp(-j*phi) ; ;
Smn1 = zeros(number_of_tx, number_of_rx);
R1 = zeros(1,number_of_tx*number_of_rx);
temp = 1;
for m = 1:1:number_of_tx
   for n = 1:1: number_of_rx
       R1(1,temp) =  R0_1(1,m)+R1_1(1,n);
       temp = temp + 1;
   end
end
received_vector1 = (row_p./R1).*exp(-(j*2*pi.*R1)/lambda); %% received signal stored in an array


%% part 3.1.5 fft calculation


max_theta =(atan(0.35));
fs = 4 * sin(max_theta) / lambda ; % sampling frequency
dx = 1/fs ; 
sample_number = 2000*length(received_vector1);% zero padding
ff1 = fftshift(fft(received_vector1,sample_number)); % do the fft shift it on the center

power1 = abs(ff1).^2 / sample_number; % calculating the PSD
x = (-fs/2 : fs/sample_number : fs/2-fs/sample_number); 
figure(2)
plot(x , power1) % plot the PSD

%% part 3.1.6 angle estimation
[M1,I1] = max(power1);
fx1 = -fs/2 + I1*fs/sample_number;% calculating the fx of the maxmumim
sine_theta1 = fx1 * lambda / (2) ;
estimated_theta1 = rad2deg(asin(sine_theta1))% estimate the theta


%% part 3.2 tow object detection

xBox= [-35,-35,35,35,-35];
yBox= [100,150,150,100,100];
figure(3)
 plot(xBox,yBox)
 hold on
a=rand;
target_x2 = -35+70*rand; %In general, you can generate N random numbers in the interval (a,b) with the formula r = a + (b-a).*rand(N,1).
target_y2 = 100+50*rand;
plot (target_x1, target_y1, '+' )
plot (target_x2, target_y2, '*' )
xlim ([- 80 80])
ylim ([-10 150])

hold off

theta2 =  rad2deg(atan(target_x2 / target_y2));


R0_2 = zeros(1,number_of_tx);
R1_2= zeros(1,number_of_rx);
for m = 1:1:number_of_tx
   R0_2(1,m) = sqrt((target_x2 -(tx_position(1,m)))^2 + target_y2^2); % calculate the R0
end

for n = 1:1:number_of_rx
     R1_2(1,n) = sqrt((target_x2 -rx_position(1,n))^2 + target_y2^2); % calculate the R2
end
R2 = zeros(1,number_of_tx*number_of_rx);
temp = 1;
for m = 1:1:number_of_tx
   for n = 1:1: number_of_rx
       R2(1,temp) =  R0_2(1,m)+R1_2(1,n);
       temp = temp + 1;
   end
end

received_vector2 = (row_p./R2).*exp(-(j*2*pi.*R2)/lambda);% 

ff2 = fftshift(fft(received_vector2,sample_number));
power2 = abs(ff2.^2)/sample_number;
x = (-fs/2 : fs/sample_number : fs/2-fs/sample_number);
figure(4)
plot(x,power2)

received_vector = received_vector1 + received_vector2; % received signal of 2 object
ff = fftshift(fft(received_vector,sample_number));
power = abs(ff.^2)/sample_number;
figure(5)
plot(x,power)

% finding 2 peaks and their cooresponding fx 

[local_peaks location] = findpeaks(power);

[max1 , index1] = max(local_peaks);
location1 = location(index1);
local_peaks(index1) = [];
location(index1) = [];

[max2 , index2] = max(local_peaks);

location2 = location(index2);
fx_1 = -fs/2 + location1*fs/sample_number;
sine_theta_1 = fx_1 * lambda / (2) ;
estimated_theta_1 = rad2deg(asin(sine_theta_1));

fx_2 = -fs/2 + location2*fs/sample_number;
sine_theta_2 = fx_2 * lambda / (2) ;
estimated_theta_2 = rad2deg(asin(sine_theta_2));
two_obj_theta = [theta1 theta2]
two_object_estimated_theta = [  estimated_theta_1 estimated_theta_2] 
%%

x_mesh = -35 :  dx : 35;
y_mesh = 100 : dx : 150;
[x,y] = meshgrid(x_mesh,y_mesh);% creat the pixel
Ntot = number_of_rx * number_of_tx;

if mode(Ntot,2) == 0
    xn = dx*((Ntot+1)/2)
    virtual_antenna_position = (-xn+dx : dx : xn-dx );
     
else
    xn = ((Ntot-1)/2)*dx;  
    virtual_antenna_position = (-xn : dx : xn);   
    
end



Sp = zeros(length(y_mesh),length(x_mesh));

for l = 1:1:length(received_vector1)
    D=sqrt((x-virtual_antenna_position(l)).^2+(y).^2);
    S = received_vector1(1,l)*D.*exp(j*4*pi*D/lambda); % propagate the signal of each virtual antena
    Sp = Sp + S ; % sum the signals propagated by each virtual antennas

end
figure(6)
imagesc(x_mesh,y_mesh,real(Sp))


K = abs(Sp);
maximum = max(max(K));
[alfa,beta]=find(K==maximum);
rad2deg(atan(x(alfa,beta)/y(alfa,beta))) % estimate the theta