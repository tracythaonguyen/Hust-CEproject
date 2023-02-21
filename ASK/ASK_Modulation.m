% <<<<<<<<<<<<<<<<<<<< ASK Modulation and Demodulation >>>>>>>>>>>>>>>>>>>>

clc, clear all, close all;

% ******************* Digital/Binary input information ********************
num_bit = 10;
x = randi(2, [1,num_bit], 'int32') - 1; % Randomly generate a Binary information as stream of bits (binary signal 0 or 1)
x1 = x; % x1 stores the input binary signal 
N = length(x); % length of the binary information
Tb = 0.001;   % one bit lasts in Tb second or bit rate = 1/Tb
Rb = 1/Tb; % Rb = bit rate 


% ************* Represent input information as digital signal *************

% Tb: Thời gian truyền 1 bit
% nb: 1 bit được biểu diễn bởi nb signal (số mẫu lấy trong 1 chu kì truyền 1 bit)
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

% plot the input binary information
t1 = Tb/nb : Tb/nb : nb*N*(Tb/nb);   % time axis
figure('Name','2-ASK Modulation and Demodulation With Noise','NumberTitle','off');
subplot(7,1,1);
plot(t1,digit,'LineWidth',2.5);
grid on;
axis([0 Tb*N -0.5 1.5]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('Digital Input Signal');

% *************************** ASK Modulation *****************************
Eb = 0.5; % energy per bit (V^2)

% 100% depth BASK unipolar modulation
Ac1 = 2*sqrt(Eb/Tb);     % Carrier amplitude for binary input '1'
Ac2 = 0;      % Carrier amplitude for binary input '0'

% 50% depth BASK 
% A = 4* sqrt(Eb/ (5*Tb));
% Ac1 = A;    % Carrier amplitude for binary input '1'
% Ac2 = A/2;      % Carrier amplitude for binary input '0'

Fc = Rb*10;   % Carrier frequency
mod = [];

y1 = Ac1*cos(2*pi*Fc*t1); % carrier function 1
y2 = Ac2*cos(2*pi*Fc*t1); % carrier function 2

% plot carrier 1
subplot(7,1,2);
plot(t1,y1);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('2-ASK Carrier 1');

% plot carrier 2
subplot(7,1,3);
plot(t1,y2);
ylim([-20,20]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('2-ASK Carrier 2');

t2 = Tb/nb:Tb/nb:Tb;   % Signal time
for (i = 1:1:N)
    if (x(i) == 1)
        y = Ac1*cos(2*pi*Fc*t2);   % Modulation signal with carrier signal 1
    else
        y = Ac2*cos(2*pi*Fc*t2);   % Modulation signal with carrier signal 2
    end
    mod = [mod y];
end

% plot modulated signal
subplot(7,1,4);
plot(t1,mod);
grid on;
ylim([-20,20]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('ASK Modulated Signal');

% ************************** End of the program ***************************
