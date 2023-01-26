%<<<<<<<<<<<<<<<<<<< BPSK Modulation and Demodulation >>>>>>>>>>>>>>>>>>>
clc, clear all, close all;
%******************* Digital/Binary input information ********************
disp('Enter Digital Input Information = '); % Binary information as stream of bits (binary signal 0 or 1)
x = randi(2, [1,10], 'int32') - 1;
x1 = x;
N = length(x);
Tb = 0.0001;   %Data rate = 1MHz i.e., bit period (second)
disp('Binary Input Information at Transmitter: ');
disp(x);
%************* Represent input information as digital signal *************
nb = 20;   % Digital signal per bit
digit = []; 
for n = 1:1:N    
    if x(n) == 1;   
       sig = ones(1,nb);
    else x(n) == 0;
        sig = zeros(1,nb);
    end
     digit = [digit sig];
end
t1=Tb/nb:Tb/nb:nb*N*(Tb/nb);   % Time period 
figure('Name','BPSK Modulation and Demodulation','NumberTitle','off');
subplot(4,1,1);
plot(t1,digit,'lineWidth',2.5);
grid on;
axis([0 Tb*N -0.5 1.5]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('Digital Input Signal');
%**************************** BPSK Modulation **************************** 
Ac = 1;      % Carrier amplitude for binary input  
br = 1/Tb;    % Bit rate
N0 = sqrt(10);
Fc = br;   % Carrier frequency 
Pc1 = 0;      % Carrier phase for binary input '1'
Pc2 = pi;     % Carrier phase for binary input '0'
t2 = Tb/nb:Tb/nb:Tb;   % Signal time                 
mod = [];
for (i = 1:1:N)
    if (x(i)==1)
        y = Ac*cos(2*pi*Fc*t2+Pc1);   % Modulation signal with carrier signal 1
    else
        y = Ac*cos(2*pi*Fc*t2+Pc2);   % Modulation signal with carrier signal 2
    end
    mod=[mod y];
end
t3=Tb/nb:Tb/nb:Tb*N;   % Time period
subplot(4,1,2);
plot(t3,mod);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('BPSK Modulated Signal');
%********************* Transmitted signal x ******************************
x = mod;
%********************* Channel model h and w *****************************
h = 1;   % Signal fading 
w = 0;   % Noise
%********************* Received signal y *********************************
y = h.*x + w;   % Convolution

t3=Tb/nb:Tb/nb:Tb*N;   % Time period
subplot(4,1,3);
plot(t3,y);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('BPSK Demodulated Signal');

%*************************** BPSK Demodulation ***************************
s = length(t2);
demod = []; 
for n = s:s:length(y)
  t4 = Tb/nb:Tb/nb:Tb;
  c = cos(2*pi*Fc*t4);      % carrier signal 
  mm = c.*y((n-(s-1)):n);   % Convolution
  t5 = Tb/nb:Tb/nb:Tb;
  z = trapz(t5,mm);         % integral 
  rz = round((2*z/Tb));                                     
  if(rz > Ac/2)             % Logical condition 
    a = 1;
  else
    a = 0;
  end
  demod = [demod a];
end 
disp('Demodulated Binary Information at Receiver: ');
disp(demod);
%********** Represent demodulated information as digital signal **********
digit = [];
for n = 1:length(demod);
    if demod(n) == 1;
       sig = ones(1,nb);
    else demod(n) == 0;
        sig = zeros(1,nb);
    end
     digit = [digit sig];
end
t5=Tb/nb:Tb/nb:nb*length(demod)*(Tb/nb);   % Time period
subplot(4,1,4)
plot(t5,digit,'LineWidth',2.5);grid on;
axis([0 Tb*length(demod) -0.5 1.5]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('BPSK Demodulated Binary Data');
num = xor(demod,x1);
disp(sum(num));
%************************** End of the program ***************************





