% %*******************************************************************************
% %********************************* TIME DOMAIN *******************************
% %*******************************************************************************
% 
clc ;
clear all;
% 
% % ================================
% % ======== reading audios ========
% % ================================
% 
[x1,fs1]=audioread('audio1.wav');
[x2,fs2]=audioread('audio2.wav');
[x3,fs3]=audioread('audio3.wav');

x1=interp(x1,12);
x2=interp(x2,12);
x3=interp(x3,12);

fs1=fs1*12;
TS1 = 1/fs1;                            %Time sampling for Second signal
N1 = [-length(x1)/2:(length(x1)/2-1)];
L1 = length(x1);                        %length of second signal
f1 = [-L1/2:L1/2-1]*(fs1/L1);           %vector for frequancy sampling for the second

fs2=fs2*12;
TS2 = 1/fs2;                            %Time sampling for Second signal
N2 = [-length(x2)/2:(length(x2)/2-1)];
L2 = length(x2);                        %length of second signal
f2 = [-L2/2:L2/2-1]*(fs2/L2);           %vector for frequancy sampling for the second

fs3=fs3*12;
TS3 = 1/fs3;                            %Time sampling for Second signal
N3 = [-length(x3)/2:(length(x3)/2-1)];
L3 = length(x3);                        %length of second signal
f3 = [-L3/2:L3/2-1]*(fs3/L3);           %vector for frequancy sampling for the second

%figure(1);
%subplot(3,1,1);
%plot(f1,abs(fftshift(fft(x1))));
%subplot(3,1,2);
%plot(f2,abs(fftshift(fft(x2))));
%subplot(3,1,3);
%plot(f3,abs(fftshift(fft(x3))));

% % ================================
% % ========== Modulating ==========
% % ================================
Fc1=10000;
carrier1 = cos(2*pi*Fc1*TS1*N1);
x1_modulated = x1'.* carrier1;
Fc2=100000;
carrier2 = cos(2*pi*Fc2*TS2*N2);
x2_modulated = x2'.*carrier2;
Fc3=100000;
carrier3 = sin(2*pi*Fc3*TS3*N3);
x3_modulated = x3'.*carrier3;

figure(2);
subplot(3,1,1);
plot(f1,abs(fftshift(fft(x1_modulated))));
subplot(3,1,2);
plot(f2,abs(fftshift(fft(x2_modulated))));
subplot(3,1,3);
plot(f3,abs(fftshift(fft(x3_modulated))));



 
% %*******************************************************************************
% %******************************* COMBINED SIGNAL *****************************
% %******************************************************************************* 

x1_modulated=[x1_modulated, zeros(1,L2-L1)];
x2_modulated=[x2_modulated, zeros(1,L2-L2)];
x3_modulated=[x3_modulated, zeros(1,L2-L3)];
signal = x1_modulated + x2_modulated + x3_modulated ;
figure();
plot(signal);
title('Time Domain')
figure();
plot(f2,abs(fftshift(fft(signal))));
title('Frequency Domain')

[C1,R1]=size(x1);
[C2,R2]=size(x2);
[C3,R3]=size(x3);

% %*******************************************************************************
% %******************************* Demodulate ************************************
% %*******************************************************************************

x1_new = signal(1:R1, 1:C1) .* carrier1;
lp = designfilt('lowpassfir', 'FilterOrder',64, 'CutoffFrequency',5*10^3, 'SampleRate', fs1);%lowpassfilter
x1_new=filter(lp,x1_new);
% [b,a] = butter(4,0.9,'low');
% x1_new = filter(b,a,x1_new);
%figure;
%plot(f1,abs(fftshift(fft(x1_new))));

OUTPUT1=downsample(x1_new,12);
%sound(OUTPUT1,fs1/12);
audiowrite('audio1Out.wav',OUTPUT1,fs1/12);



x2_new = signal(1:R2, 1:C2) .* carrier2;
lp = designfilt('lowpassfir', 'FilterOrder',128, 'CutoffFrequency',10*10^3, 'SampleRate', fs2);%lowpassfilter
x2_new=filter(lp,x2_new);
% [b,a] = butter(4,0.9,'low');
% x2_new = filter(b,a,x2_new);
%figure;
%plot(f2,x2_new);%abs(fftshift(fft(x2_new))));
OUTPUT2=downsample(x2_new,12);
%sound(OUTPUT2,fs2/12);
audiowrite('audio2Out.wav',OUTPUT2,fs2/12);


x3_new = signal(1:R3, 1:C3) .* carrier3;
lp = designfilt('lowpassfir', 'FilterOrder',2048, 'CutoffFrequency',10*10^3, 'SampleRate', fs3);%lowpassfilter
x3_new=filter(lp,x3_new);
% [b,a] = butter(4,0.9,'low');
% x2_new = filter(b,a,x2_new);
%figure;
%plot(f3,x3_new);%abs(fftshift(fft(x2_new))));
OUTPUT3=downsample(x3_new,12);
%sound(OUTPUT3,fs3/12);
audiowrite('audio3Out.wav',OUTPUT3,fs3/12);




% %*******************************************************************************************
% %******************************* Demodulate Phase Shift ************************************
% %*******************************************************************************************

%90

carrier1 = cos(2*pi*Fc1*TS1*N1 + pi/2);
x1_new = signal(1:R1, 1:C1) .* carrier1;
lp = designfilt('lowpassfir', 'FilterOrder',64, 'CutoffFrequency',5*10^3, 'SampleRate', fs1);%lowpassfilter
x1_new=filter(lp,x1_new);
% [b,a] = butter(4,0.9,'low');
% x1_new = filter(b,a,x1_new);
%figure;
%plot(f1,abs(fftshift(fft(x1_new))));

OUTPUT1=downsample(x1_new,12);
%sound(OUTPUT1,fs1/12);
audiowrite('audio1Out90.wav',OUTPUT1,fs1/12);


carrier2 = cos(2*pi*Fc2*TS2*N2 + pi/2);
x2_new = signal(1:R2, 1:C2) .* carrier2;
lp = designfilt('lowpassfir', 'FilterOrder',128, 'CutoffFrequency',10*10^3, 'SampleRate', fs2);%lowpassfilter
x2_new=filter(lp,x2_new);
% [b,a] = butter(4,0.9,'low');
% x2_new = filter(b,a,x2_new);
%figure;
%plot(f2,x2_new);%abs(fftshift(fft(x2_new))));
OUTPUT2=downsample(x2_new,12);
%sound(OUTPUT2,fs2/12);
audiowrite('audio2Out90.wav',OUTPUT2,fs2/12);

carrier3 = sin(2*pi*Fc3*TS3*N3 + pi/2);
x3_new = signal(1:R3, 1:C3) .* carrier3;
lp = designfilt('lowpassfir', 'FilterOrder',2048, 'CutoffFrequency',10*10^3, 'SampleRate', fs3);%lowpassfilter
x3_new=filter(lp,x3_new);
% [b,a] = butter(4,0.9,'low');
% x2_new = filter(b,a,x2_new);
%figure;
%plot(f3,x3_new);%abs(fftshift(fft(x2_new))));
OUTPUT3=downsample(x3_new,12);
%sound(OUTPUT3,fs3/12);
audiowrite('audio3Out90.wav',OUTPUT3,fs3/12);

%30

carrier1 = cos(2*pi*Fc1*TS1*N1 + pi/6);
x1_new = signal(1:R1, 1:C1) .* carrier1;
lp = designfilt('lowpassfir', 'FilterOrder',64, 'CutoffFrequency',5*10^3, 'SampleRate', fs1);%lowpassfilter
x1_new=filter(lp,x1_new);
% [b,a] = butter(4,0.9,'low');
% x1_new = filter(b,a,x1_new);
%figure;
%plot(f1,abs(fftshift(fft(x1_new))));

OUTPUT1=downsample(x1_new,12);
%sound(OUTPUT1,fs1/12);
audiowrite('audio1Out30.wav',OUTPUT1,fs1/12);


carrier2 = cos(2*pi*Fc2*TS2*N2 + pi/6);
x2_new = signal(1:R2, 1:C2) .* carrier2;
lp = designfilt('lowpassfir', 'FilterOrder',128, 'CutoffFrequency',10*10^3, 'SampleRate', fs2);%lowpassfilter
x2_new=filter(lp,x2_new);
% [b,a] = butter(4,0.9,'low');
% x2_new = filter(b,a,x2_new);
%figure;
%plot(f2,x2_new);%abs(fftshift(fft(x2_new))));
OUTPUT2=downsample(x2_new,12);
%sound(OUTPUT2,fs2/12);
audiowrite('audio2Out30.wav',OUTPUT2,fs2/12);

carrier3 = sin(2*pi*Fc3*TS3*N3 + pi/6);
x3_new = signal(1:R3, 1:C3) .* carrier3;
lp = designfilt('lowpassfir', 'FilterOrder',2048, 'CutoffFrequency',10*10^3, 'SampleRate', fs3);%lowpassfilter
x3_new=filter(lp,x3_new);
% [b,a] = butter(4,0.9,'low');
% x2_new = filter(b,a,x2_new);
%figure;
%plot(f3,x3_new);%abs(fftshift(fft(x2_new))));
OUTPUT3=downsample(x3_new,12);
%sound(OUTPUT3,fs3/12);
audiowrite('audio3Out30.wav',OUTPUT3,fs3/12);

%10

carrier1 = cos(2*pi*Fc1*TS1*N1 + pi/18);
x1_new = signal(1:R1, 1:C1) .* carrier1;
lp = designfilt('lowpassfir', 'FilterOrder',64, 'CutoffFrequency',5*10^3, 'SampleRate', fs1);%lowpassfilter
x1_new=filter(lp,x1_new);
% [b,a] = butter(4,0.9,'low');
% x1_new = filter(b,a,x1_new);
%figure;
%plot(f1,abs(fftshift(fft(x1_new))));

OUTPUT1=downsample(x1_new,12);
%sound(OUTPUT1,fs1/12);
audiowrite('audio1Out10.wav',OUTPUT1,fs1/12);


carrier2 = cos(2*pi*Fc2*TS2*N2 + pi/18);
x2_new = signal(1:R2, 1:C2) .* carrier2;
lp = designfilt('lowpassfir', 'FilterOrder',128, 'CutoffFrequency',10*10^3, 'SampleRate', fs2);%lowpassfilter
x2_new=filter(lp,x2_new);
% [b,a] = butter(4,0.9,'low');
% x2_new = filter(b,a,x2_new);
%figure;
%plot(f2,x2_new);%abs(fftshift(fft(x2_new))));
OUTPUT2=downsample(x2_new,12);
%sound(OUTPUT2,fs2/12);
audiowrite('audio2Out10.wav',OUTPUT2,fs2/12);

carrier3 = sin(2*pi*Fc3*TS3*N3 + pi/18);
x3_new = signal(1:R3, 1:C3) .* carrier3;
lp = designfilt('lowpassfir', 'FilterOrder',2048, 'CutoffFrequency',10*10^3, 'SampleRate', fs3);%lowpassfilter
x3_new=filter(lp,x3_new);
% [b,a] = butter(4,0.9,'low');
% x2_new = filter(b,a,x2_new);
%figure;
%plot(f3,x3_new);%abs(fftshift(fft(x2_new))));
OUTPUT3=downsample(x3_new,12);
%sound(OUTPUT3,fs3/12);
audiowrite('audio3Out10.wav',OUTPUT3,fs3/12);