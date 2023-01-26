% <<<<<<<<<<<<<<<<<<<< FSK Modulation and Demodulation >>>>>>>>>>>>>>>>>>>>

clc, clear all, close all;
% ******************* Digital/Binary input information ********************

% x = input('Enter Digital Input Information = ');   % Binary information as stream of bits (binary signal 0 or 1)
x = randi(2, [1,10], 'int32') - 1; % auto generate the binary sequence
x1 = x
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
figure('Name','ASK Modulation and Demodulation With Noise','NumberTitle','off');
subplot(6,1,1);
plot(t1,digit,'LineWidth',2.5);
grid on;
axis([0 Tb*N -0.5 1.5]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('Digital Input Signal');

% *************************** FSK Modulation *****************************
Ac = 1;      % Carrier amplitude for binary input
br = 1/Tb;    % Bit rate
Fc1 = br;      % Carrier phase for binary input '1'
Fc2 = br*10;     % Carrier phase for binary input '0' 
N0 = sqrt(10);
t2 = Tb/nb:Tb/nb:Tb;   % Signal time
mod = [];

x_ = Tb/nb:Tb/nb:Tb*N; % Time period
y1 = Ac*cos(2*pi*Fc1*x_); % carrier function 1
y2 = Ac*cos(2*pi*Fc2*x_); % carrier function 2

% plot carrier 1
subplot(6,1,2);
plot(x_,y1);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('2-FSK Carrier 1');

% plot carrier 2
subplot(6,1,3);
plot(x_,y2);
ylim([-20,20]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('2-FSK Carrier 2');

for (i = 1:1:N)
    if (x(i) == 1)
        y = Ac*cos(2*pi*Fc1*t2);   % Modulation signal with carrier signal 1
    else
        y = Ac*cos(2*pi*Fc2*t2);   % Modulation signal with carrier signal 2
    end
    mod = [mod y];
end

t3 = Tb/nb:Tb/nb:Tb*N;   % Time period
subplot(6,1,4);
plot(t3,mod);
grid on;
ylim([-20,20]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('FSK Modulated Signal');

% ********************* Transmitted signal x ******************************
x = mod;

% ********************* Channel model h and w *****************************
h = 1;   % Signal fading 
w = sqrt(N0/2) * rand(1, N*nb);   % Noise

% ********************* Received signal y *********************************
y = h.*x + w;   % Convolution

t3=Tb/nb:Tb/nb:Tb*N;   % Time period
subplot(6,1,5);
plot(t3,y);
grid on;
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('FSK Demodulated Signal');

% *************************** FSK Demodulation ****************************
s = length(t2); % số cặp giá trị (t,y) để biểu diễn 1 tín hiệu
demod = [];
Fc = br*10;
for n = s:s:length(y) 
    t4 = Tb/nb:Tb/nb:Tb;    % Time period
    c = cos(2*pi*Fc*t4);    % Carrier signal % basis of M
    mm = c.*y((n-(s-1)):n); % Convolution 
    t5 = Tb/nb:Tb/nb:Tb;
    z = trapz(t5,mm);       % Intregation 
    rz = round(2*z/Tb); % khong hieu sao phải *2/Tb
    if (rz > Ac/2)
        a = 0;
    else   
        a = 1;
    end
    demod = [demod a];
end
disp('Demodulated Binary Information at Receiver: ');
disp(demod);

% ********** Represent demodulated information as digital signal **********


digit = [];
for n = 1:length(demod);
    if demod(n) == 1;
       sig = ones(1,nb);
    else demod(n) == 0;
        sig = zeros(1,nb);
    end
     digit = [digit sig];
end

t5 = Tb/nb:Tb/nb:nb*length(demod)*(Tb/nb);   % Time period
subplot(6,1,6)
plot(t5,digit,'LineWidth',2.5);grid on;
axis([0 Tb*length(demod) -0.5 1.5]);
grid on;
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('FSK Demodulated Binary Data');
num = xor(demod,x1);
disp(sum(num));
% ************************** End of the program ***************************


