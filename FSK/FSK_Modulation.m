% <<<<<<<<<<<<<<<<<<<< FSK Modulation and Demodulation >>>>>>>>>>>>>>>>>>>>

clc, clear all, close all;
% ******************* Digital/Binary input information ********************
  
% Binary information as stream of bits (binary signal 0 or 1)
x = randi(2, [1,30], 'int32') - 1; % auto generate the binary sequence
N = length(x);
Tb = 0.0001;   %Data rate = 1MHz i.e., bit period (second)
disp('Binary Input Information at Transmitter: ');
disp(x);

% ************* Represent input information as digital signal *************

% Tb: Thời gian truyền 1 bit
% nb: 1 bit được biểu diễn bởi nb signal (số mẫu lấy trong 1 chu kì truyền
% 1 bit)
% Tb/nb: thời gian truyền một signal

nb = 100;   % Digital signal per bit
digit = []; 
for n = 1:1:N
    if x(n) == 1;
       sig = ones(1,nb);
    else x(n) == 0;
        sig = zeros(1,nb);
    end
     digit = [digit sig];
end

t1 = Tb/nb : Tb/nb : nb*N*(Tb/nb);   % Time period
figure('Name','FSK Modulation','NumberTitle','off');
subplot(4,1,1);
plot(t1,digit,'LineWidth',2.5);
grid on;
axis([0 Tb*N -0.5 1.5]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('Digital Input Signal');

% *************************** FSK Modulation *****************************
Ac = 10;      % Carrier amplitude for binary input
br = 1/Tb;    % Bit rate
Fc1 = br*4;      % Carrier phase for binary input '1'
Fc2 = br*6;     % Carrier phase for binary input '0' 

t2 = Tb/nb:Tb/nb:Tb;   % Signal time
mod = [];

for (i = 1:1:N)
    if (x(i) == 1)
        y = Ac*cos(2*pi*Fc1*t2);   % Modulation signal with carrier signal 1
    else
        y = Ac*cos(2*pi*Fc2*t2);   % Modulation signal with carrier signal 2
    end
    mod = [mod y];
end

t3 = Tb/nb:Tb/nb:Tb*N; % Time period

y1 = Ac*cos(2*pi*Fc1*t3); % carrier function 1
y2 = Ac*cos(2*pi*Fc2*t3); % carrier function 2

% plot carrier 1
subplot(4,1,2);
plot(t3,y1);
grid on;
ylim([-20,20]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('2-FSK Carrier 1');

% plot carrier 2
subplot(4,1,3);
plot(t3,y2);
grid on;
ylim([-20,20]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('2-FSK Carrier 2');

subplot(4,1,4);
plot(t3,mod);
grid on;
ylim([-20,20]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('FSK Modulated Signal');

% ************************** End of the program ***************************