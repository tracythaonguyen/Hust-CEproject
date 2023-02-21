% <<<<<<<<<<<<<<<<<<<< ASK Modulation and Demodulation >>>>>>>>>>>>>>>>>>>>

clc, clear all, close all;
% ******************* Digital/Binary input information ********************
num_bit = 1000;
x = randi(2, [1,num_bit], 'int32') - 1; % Randomly generate a Binary information as stream of bits (binary signal 0 or 1)
x1 = x; % x1 stores the input binary signal 
N = length(x); % length of the binary information
Tb = 0.001;   % one bit lasts in Tb second or bit rate = 1/Tb
Rb = 1/Tb; % Rb = bit rate 
% disp('Binary Input Information at Transmitter: ');
% disp(x);

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
subplot(6,1,1);
axis tight;
plot(t1,digit,'LineWidth',2.5);
grid on;
axis([0 Tb*N -0.5 1.5]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('Digital Input Signal');

% *************************** ASK Modulation *****************************
dB = 4;
Eb = 5; % energy per bit (V^2)

% 100% depth BASK unipolar modulation
% Ac1 = 2*sqrt(Eb/Tb);     % Carrier amplitude for binary input '1'
% Ac2 = 0;      % Carrier amplitude for binary input '0'

% 50% depth BASK 
A = 4* sqrt(Eb/ (5*Tb));
Ac1 = A;    % Carrier amplitude for binary input '1'
Ac2 = A/2;      % Carrier amplitude for binary input '0'

Fc = Rb*10;   % Carrier frequency
mod = [];

y1 = Ac1*cos(2*pi*Fc*t1); % carrier function 1
y2 = Ac2*cos(2*pi*Fc*t1); % carrier function 2

% plot carrier 1
axis tight;
subplot(6,1,2);
plot(t1,y1);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('2-ASK Carrier 1');

% plot carrier 2
axis tight;
subplot(6,1,3);
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
axis tight;
subplot(6,1,4);
plot(t1,mod);
grid on;
ylim([-20,20]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('ASK Modulated Signal');

% ********************* Transmitted signal x ******************************
x = mod;

% ********************* Channel model h and noise *****************************
%h = 1;   % Signal fading 

% ********************* Received signal y *********************************

% Received signal without noise
y = x;

% plot Received signal without noise
% axis tight;
% subplot(6,1,5);
% plot(t1,y);
% grid on;
% xlabel('Time(Sec)');
% ylabel('Amplitude(Volts)');
% title('2-ASK Received Signal Without Noise');

% -------------------------------------------------------------------------
% Received signal with white noise
N0 = Eb/(10.^(dB/10));
mean = 0;
sigma = sqrt(N0/2);
noise = randn(size(x)); % Generate random numbers from a standard normal distribution
noise = mean + sigma*noise; % Scale the random numbers to obtain Gaussian distribution with mean and sigma
% display(noise);
y = x + noise;   % add noise

% plot Received signal with white noise
axis tight;
subplot(6,1,5);
plot(t1,y);
grid on;
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('2-ASK Received Signal With White Noise');
% axis tight;

% *************************** ASK Demodulation ****************************
s = length(t2);
demod = [];
for n = s:s:length(y) 
    t4 = Tb/nb:Tb/nb:Tb;    % Time period
    c = cos(2*pi*Fc*t4);    % basis of M
    mm = c.*y((n-(s-1)):n); % Convolution 
    t5 = Tb/nb:Tb/nb:Tb;
    z = trapz(t5,mm);       % Integration 
    rz = 2*z/Tb;
    d1 = abs(rz-Ac1);
    d2 = abs(rz-Ac2);
    if (d1 < d2)
        a = 1;
    else   
        a = 0;
    end
    demod = [demod a];
end

%disp('Demodulated Binary Information at Receiver: ');
%disp(demod);

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
axis tight;
plot(t5,digit,'LineWidth',2.5);grid on;
axis([0 Tb*length(demod) -0.5 1.5]);
grid on;
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('2-ASK Demodulated Binary Data');

% *************************** Calculate BER *******************************
num = xor(demod, x1);
fprintf('*** Compare input signal and received signal ***\n');
fprintf('Number of binary information bit: %d\n', num_bit);
fprintf('Number of incorrect bit = %i\n', sum(num));
fprintf('BER = %.4f \n', sum(num)/num_bit);

fprintf('\n*** BER Theory ***\n')
fprintf('Eb/N0 = %f\n',Eb/N0);
fprintf('10*lg(Eb/N0) = %f\n',10*log10(Eb/N0));
fprintf('P(e) = %f', 0.5 * erfc(sqrt(Eb/(4*N0))));

% ************************** End of the program ***************************


