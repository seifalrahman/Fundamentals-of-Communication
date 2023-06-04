pkg load communications; %m7tagha 3shan al Eye Diagram
%%%%%%%%%% PART 1 %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Transmitter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bit_burst_stream=randi([0,1],1,10000); % dah al generated bits bta3tna
%%%%%%%%%%%%%%%%%%%%%%%%UNI-POLAR-NRZ%%%%%%%%%%%%%%%%%%%%%%%%

Vcc=1.2; %hna 3ndna al Variables bta3na -Vcc
lvl=[0,Vcc]; %- lvl dah bn3mlo 3shan 3ndna leveleen 7rfyan xD
bit_Coded=lvl(bit_burst_stream+1); % - bit coded dah bn3mlo 3shan n5ly al bits bta3na bdl 1 wa 0 tb2a Vcc wa 0
N=1000000;           %dah al N aly hn3ml byha sample
pulse_Duration=0.01; %dah al Rs(sample rate)
 global t ; %3rfna al time global 3shan hnst3mlo kteer
 t = 0:pulse_Duration:(N-1)*pulse_Duration; %hna bnmla al time bta3na bkol pulse
TxUPNRZ=zeros(1,N);  %hna bn3ml al signal kolha esfar
fs=100;
T=10000;      %3shan 10000 bit
E=ceil(T/pulse_Duration); %bn7sb  al size
df=1/T;
if (rem(E,2)==0)
  f= -(0.5*fs):df:(0.5*fs-df); %Even
else
  f= -(0.5*fs-0.5*df):df:(0.5*fs-0.5*df); %odd
end


%%%%%%%%%

h=100; %l wa h dool 3shan hytkkro kteer mn al a5er bn3mlha kda 3shan n5lyy al bit kolha mn 1 le 100 (kol al samples y3ni)
l=1;   %fa n5ly al kol al bits value wa7ed
for i=1:length(bit_Coded) %hna ant bta5od al bits aly 3mlnlha generate
  TxUPNRZ(l:h)=bit_Coded(i); % wa bt3mlha sample 100 mra y3ni mn al a5er bta5od al bit al 1 btkrrha 100 mra
  h=h+100; l=l+100;
end

figure
plot(t,TxUPNRZ);
ylim([-3 3]); grid on;
xlabel('Time'); ylabel('Magnitude');
title('TxUPNRZ');
eyediagram(TxUPNRZ,200,1) %function al eye diagram 200 dyh 3shan al bnsample fy 100 fa 5lnah yrsm 2 mn al a5er
title('TxUPNRZ');
POTATOxUPNRZ=fftshift(fft(TxUPNRZ))*pulse_Duration; %wa dyh 3shan n3ml al spectral domain
figure              %(potato asm al spectral domain bta3na xD)
plot(f,abs(POTATOxUPNRZ).^2);grid on;%power 2 hna 3shan hoa 3awz power xD
xlabel('Frequency'); ylabel('Magnitude');
title('TxUPNRZ');

%%%%%%%%%%%%%%%%%%%%%%%%POLAR-NRZ%%%%%%%%%%%%%%%%%%%%%%%%
lvl2=[-Vcc,Vcc]; %hna 3mlna level tani 3shan dah polar fa mn -Vcc le Vcc b3dha nfs al steps mtkrra
bit_Coded2=lvl2(bit_burst_stream+1);
%%%%%%%%%
TxPNRZ=zeros(1,N);
h=100;
l=1;
for i=1:length(bit_Coded2)
  TxPNRZ(l:h)=bit_Coded2(i);
  h=h+100; l=l+100;
end
figure
plot(t,TxPNRZ);
ylim([-3 3]); grid on;
xlabel('Time'); ylabel('Magnitude');
title('TxPNRZ');
eyediagram(TxPNRZ,200,1)
title('TxPNRZ');
POTATOxPNRZ=fftshift(fft(TxPNRZ))*pulse_Duration;
figure
plot(f,abs(POTATOxPNRZ).^2);grid on;
xlabel('Frequency'); ylabel('Magnitude');
title('TxPNRZ');

%%%%%%%%%%%%%%%%%%%%%%%%UNI-POLAR-RZ%%%%%%%%%%%%%%%%%%%%%%%%
TxUPRZ=zeros(1,N);

for i=1:length(bit_Coded) %hna asmna al 100 bta3tna noseen mn al a5er 3shan n3ml nos Vcc wa al nos al tani zero
  l=(i-1)*100+1;
  h=l+49;
  TxUPRZ(l:h)=bit_Coded(i);
  l=h+1; h=l+49;
  TxUPRZ(l:h)=0;
end
figure
plot(t,TxUPRZ);
ylim([-3 3]); grid on;
xlabel('Time'); ylabel('Magnitude');
title('TxUPRZ');
eyediagram(TxUPRZ,200,1)
title('TxUPRZ');
POTATOxUPRZ=fftshift(fft(TxUPRZ))*pulse_Duration;
figure
plot(f,abs(POTATOxUPRZ).^2);grid on;
xlabel('Frequency'); ylabel('Magnitude');
title('TxUPRZ');

%%%%%%%%%%%%%%%%%%%%%%%%BI-POLAR-RZ%%%%%%%%%%%%%%%%%%%%%%%%
TxBIPRZ=zeros(1,N);
flag=1.2;
for i=1:length(bit_Coded)   %3mlna b2a mix le kol aly fat
  l=(i-1)*100+1;
  h=l+49;
  if (bit_Coded(i) ~=0)
  TxBIPRZ(l:h)=flag;
  flag=-flag;
  l=h+1; h=l+49;
  TxBIPRZ(l:h)=0;
  else
  TxBIPRZ(l:h)=0;
  l=h+1; h=l+49;
  TxBIPRZ(l:h)=0;
  end
end
figure
plot(t,TxBIPRZ);
ylim([-3 3]); grid on;
xlabel('Time'); ylabel('Magnitude');
title('TxBIPRZ');
eyediagram(TxBIPRZ,200,1)
title('TxBIPRZ');
POTATOxBIPRZ=fftshift(fft(TxBIPRZ))*pulse_Duration;
figure
plot(f,abs(POTATOxBIPRZ).^2);grid on;
xlabel('Frequency'); ylabel('Magnitude');
title('TxBIPRZ');

%%%%%%%%%%%%%%%%%%%%%%%%MANCH%%%%%%%%%%%%%%%%%%%%%%%%
TxMANCH=zeros(1,N);
flag=1.2;
for i=1:length(bit_Coded) %3mlna condition le al 0 wa condition lel 1
  l=(i-1)*100+1;
  h=l+49;
  if (bit_Coded(i) ~=0)
  TxMANCH(l:h)=flag;

  l=h+1; h=l+49;
  TxMANCH(l:h)=-flag;
  else
  TxMANCH(l:h)=-flag;
  l=h+1; h=l+49;
  TxMANCH(l:h)=flag;
  end
end
figure
plot(t,TxMANCH);
ylim([-3 3]); grid on;   %lel bys2l eh -3 wa 3 dah by5ly fyh max fo2 wa t7t lel graph awl ma yzhr fy wshk 3shan al mzhr ykon latef
xlabel('Time'); ylabel('Magnitude');
title('TxMANCH');
eyediagram(TxMANCH,200,1)
title('TxMANCH');
POTATOxMANCH=fftshift(fft(TxMANCH))*pulse_Duration;
figure
plot(f,abs(POTATOxMANCH).^2);grid on;
xlabel('Frequency'); ylabel('Magnitude');
title('TxMANCH');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%n5osh fy al receive b2a

%%
sigma=linspace(0,Vcc,10); %dyh al sigma aly hoa tlabha mnna linspace dyh bt3ml an minimum 0 wa max vcc wa 10 7agat
BERxUPNRZ =zeros(1,10); %b3ml Bit error rate vector le kol no3
BERxPNRZ =zeros(1,10);
BERxUPRZ =zeros(1,10);
BERxBIPRZ =zeros(1,10);
BERxMANCH =zeros(1,10);
%%%%
function [Rx]=Receiver(signal,type) %3mlna al receive function wa hnst3mlha t7t (dyh zyada mn 3ndna mlhash msh lazem function y3ni)
Rx=zeros(1,1000000);
Vcc=1.2;
Rlvl=[-Vcc/2,0,Vcc/2]; %3 levels 3shan ne-detect byhom
%%%%%%%%%%%%%%%%%%
if (type==1)
l=1;
h=100;
for i=50:100:length(signal) %hna bna5od nos al signal dyh aly hnshof 3ndha al value bta3 al signal 3shan a7na 3amlyn sample b 100
if(signal(i)>=Rlvl(3))      % wa bnmshy 100 100 3shan ant 3ndk kol bit wa5da 100 lw7dha
Rx(l:h)=Vcc;      %3amlyn hna lw al signla akbr mn Vcc/2 n3mlha Vcc else 0
else
Rx(l:h)=0;
end
l=l+100;
h=h+100;
end
end
%%%%%%%%%%%%%%%%%%
if (type==2) %type 2 bta3 al PNRZ
l=1;
h=100;
for i=50:100:length(signal)
if(signal(i)>Rlvl(2))
Rx(l:h)=Vcc;
else
Rx(l:h)=-Vcc;
end
l=l+100;
h=h+100;
end
end
%%%%%%%%%%%%%%%%%%
if (type==3) %type 3 bta3 al UPRZ
l=1;
h=50;
for i=50:100:length(signal)
if(signal(i)>=Rlvl(3))
Rx(l:h)=Vcc;
Rx((h+1):(h+50))=0;
else
Rx(l:(h+50))=0;
end
l=l+100;
h=h+100;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(type==4) %type 4 bta3 al BIPRZ
l=1;
h=50;
for i=50:100:length(signal)
if(signal(i)>Rlvl(3) && signal(i)~=0 )
Rx(l:h)=Vcc;
Rx(h+1:h+50)=0;
elseif(signal(i)<Rlvl(1) && signal(i)~=0)
Rx(l:h)=-Vcc;
Rx(h+1:h+50)=0;
else
signal(l:h+50)=0;
end
l=l+100;
h=h+100;

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(type==5) %type 5 MANCHESTER
l=1;
h=50;

for i=50:100:length(signal) %2 conditions lel al 0 wa al 1
if(signal(i)>Rlvl(2))
Rx(l:h)=Vcc;
Rx(h+1:h+50)=-Vcc;
else
Rx(l:h)=-Vcc;
Rx(h+1:h+50)=Vcc;
end
l=l+100;
h=h+100;
end
end
end


%%%%%


for i=1:10
n = sigma(i) * randn(1,length(t) );% hna bncreate al noise bl fucntion aly Dr mdyhalna `
signal_With_Noise=n+TxUPNRZ; %b7oto 3la al signal aly transmitted
RxUPNRZ = Receiver(signal_With_Noise,1); %bnst3ml al function aly 3mlnaha fo2 dyh (mfesh prototype fy Octave 3shan kda 7atenha fo2)
noerror = sum(RxUPNRZ~=TxUPNRZ); %bn3d al errors 3shan nst3mlha fy Bit error rate
ber=noerror/N; %al 3dd 3la al size 3shan ntl3 al rate
BERxUPNRZ(i)=ber; %bngm3 al error fy Vector
disp(['BERxUPNRZ: ' num2str(ber)]);
end
figure
plot(t,RxUPNRZ);
ylim([-3 3]); grid on;
xlabel('Time'); ylabel('Magnitude');
title('RxUPNRZ');


%%
for i=1:10  %nfs al klam aly fo2 bnkrro bs 3la 7sb no3 al signal
n = sigma(i) * randn(1,length(t) );
signal_With_Noise=n+TxPNRZ;
RxPNRZ = Receiver(signal_With_Noise,2);
noerror = sum(RxPNRZ~=TxPNRZ);
ber=noerror/N;
BERxPNRZ(i)=ber;
disp(['BERxPNRZ: ' num2str(ber)]);
end
figure
plot(t,RxPNRZ);
ylim([-3 3]); grid on;
xlabel('Time'); ylabel('Magnitude');
title('RxPNRZ');


%%


for i=1:10
n = sigma(i) * randn(1,length(t) );
signal_With_Noise=n+TxUPRZ;
RxUPRZ = Receiver(signal_With_Noise,3);
noerror = sum(RxUPRZ~=TxUPRZ);
ber=2*noerror/N;
BERxUPRZ(i)=ber;
disp(['BERxUPRZ: ' num2str(ber)]);
end
figure
plot(t,RxUPRZ);
ylim([-3 3]); grid on;
xlabel('Time'); ylabel('Magnitude');
title('RxUPRZ');

%%
Errors=zeros(1,10); % Vector hshyl fyh al errors
Eflag=0; %%current bit mn al a5er %hna dyh goz2yt al bonus bta3 al bi polar
next=0; %next bit
previous=0; %previous bit
Ecount=0; % Counter lel errors dah aly e7na 3awzeeno mn kol dah

for i=1:10
n = sigma(i) * randn(1,length(t) );
signal_With_Noise=n+TxBIPRZ;
RxBIPRZ = Receiver(signal_With_Noise,4);
noerror = sum(RxBIPRZ~=TxBIPRZ);
ber=2*noerror/N;
BERxBIPRZ(i)=ber;
disp(['BERxBIPRZ: ' num2str(ber)]);
Ecount=0;
for ii=50:100:length(RxBIPRZ)-100; %bsm allah al logic behind al klam dah fy voice 3shan mt3bsh nfsi xD
  if(RxBIPRZ(ii)>= Vcc/2 && RxBIPRZ(ii)~=0)
  Eflag=Vcc;
  next=RxBIPRZ(ii+100);
  if(ii~=50)
  previous=RxBIPRZ(ii-100);
  end
elseif(RxBIPRZ(ii)<= -Vcc/2 && RxBIPRZ(ii)~=0)
  Eflag= -Vcc;
  next=RxBIPRZ(ii+100);
  if(ii~=50)
  previous=RxBIPRZ(ii-100);
  end
end

if (Eflag==next && Eflag~=0)
  Ecount=Ecount+1;
elseif(Eflag~=next && Eflag==previous && Eflag~=0)
 Ecount=Ecount+1;
end
end


Errors(i)=Ecount;
end


figure
plot(t,RxBIPRZ);
ylim([-3 3]); grid on;
xlabel('Time'); ylabel('Magnitude');
title('RxBIPRZ');


%%
for i=1:10 % r3gna lel 3adi hn3ml al reciever bta3 al manchester
n = sigma(i) * randn(1,length(t));
signal_With_Noise=n+TxMANCH;
RxMANCH = Receiver(signal_With_Noise,5);
noerror = sum(RxMANCH~=TxMANCH);
ber=noerror/N;
BERxMANCH(i)=ber;
disp(['BER: ' num2str(ber)]);
end
figure
plot(t,RxMANCH);
ylim([-3 3]); grid on;
xlabel('Time'); ylabel('Magnitude');
title('RxMANCH');

%%

figure %wa dyh aly hya log graph bta3 al Bit Error Rate le kol wa7da fy graph wa7ed
semilogy(sigma,BERxUPNRZ);
hold on;
semilogy(sigma,BERxPNRZ);
semilogy(sigma,BERxUPRZ);
semilogy(sigma,BERxBIPRZ);
semilogy(sigma,BERxMANCH);
legend UPNRZ PNRZ UPRZ BIPRZ MANCH;
title('BER vs Sigma')
hold off;
%bs kda tm B7md allah %S7s n Yosry <3
%%%%%%%PART2%%%%%%%%%%%%%%%%%
bit__stream2=randi([0,1],1,100);%%stream of random bits
lvl__2=[-1,1]; %hna 3mlna level tani 3shan dah polar fa mn -Vcc le Vcc b3dha nfs al steps mtkrra
bit__Coded__2=lvl__2(bit__stream2+1);
%%%%%%%%%
N__2=10000 ; %%total bits sampling
ts = 1e-11 ; 
tb = 1e-9;
tt = 0:ts:(N__2-1)*ts;%%1e-7 step 1e-11
N__2__2=100;
tc = 0:ts:(N__2__2-1)*ts;%%1e-9 step 1e-11
TxPNRZ__2=zeros(1,N__2);
hh=100;
ll=1;

for j=1:length(bit__Coded__2)
  TxPNRZ__2(ll:hh) =  bit__Coded__2(j);%%%%%1e4
  hh=hh+100; ll=ll+100;
end

fss=100;
TT=100;      %3shan 100 bit
pulse_Duration_2=0.01;
EE=ceil(TT/pulse_Duration_2); #bn7sb  al size
dff=1/TT;
if (rem(EE,2)==0)
  ff= -(0.5*fss):dff:(0.5*fss-dff); %Even
else
  ff= -(0.5*fss-0.5*dff):dff:(0.5*fss-0.5*dff); %odd
end

%%%%%TXPNRZ%%%%%%%
%%TIME%%%
figure
plot(tt,TxPNRZ__2);
ylim([-3 3]); grid on;
xlabel('Time'); ylabel('Magnitude');
title('TxPNRZ__part2');
%%%FREQ%%%%
POTATOxPNRZ__2=fftshift(fft(TxPNRZ__2))*pulse_Duration_2;
figure
plot(ff,abs(POTATOxPNRZ__2));grid on;
xlabel('Frequency'); ylabel('Magnitude');
title('TxPNRZ part2 freq');
%%%%%%%%POWER%%%%%
POTATOxPNRZ__2=fftshift(fft(TxPNRZ__2))*pulse_Duration_2;
figure
plot(ff,abs(POTATOxPNRZ__2).^2);grid on;
xlabel('Frequency'); ylabel('Magnitude');
title('TxPNRZ part2 power');

%%%%%%%%MODULATION%%%%%%

fc= 1000000000;
phi = cos(2* pi *fc*tc) ;%%%%eb=tb/2
mod =[];
for i=1:length(bit__Coded__2)
  if (bit__Coded__2(i) > 0)
  modulated_BPSK=phi;
  else
   modulated_BPSK=-phi;
  end
  mod = [mod modulated_BPSK];
end

%%%%%TIME%%%%%%%%
figure
plot(tt,mod);
ylim([-3 3]); grid on;
xlabel('Time'); ylabel('Magnitude');
title('modulated BPSK');
axis ([0 1e-8 -2 2]); %%%%ZOOM
%%%%%%%FREQ%%%%
figure
freq_modulated_BPSK=fftshift(fft(mod))*pulse_Duration_2;
plot(ff,abs(freq_modulated_BPSK));grid on;
xlabel('Frequency'); ylabel('Magnitude');
title('freq modulated BPSK') ;
%%%%%%POWER%%%%
figure
freq_modulated_BPSK=fftshift(fft(mod))*pulse_Duration_2;
plot(ff,abs(freq_modulated_BPSK).^2);grid on;
xlabel('Frequency'); ylabel('Magnitude');
title('Power modulated BPSK') ;

%%%%%%%%%%%%%%%%%%%DEMODULATION%%%%%%%%%%%%%%%%%%%%%%%%%%
Demod = [];
for i=100:100:length(mod)
  de_mod=phi.*mod((i-99):i); %%convolution
  x = trapz(tc,de_mod); %%LPF
  rd = round (2*x/tb);
  if (rd > 0)
    y = 1;
  else
    y = -1;
  end
  Demod = [Demod y]; %%100bit
end

Demod__2=zeros(1,N__2);
hhh=100;
lll=1;
for j=1:length(bit__Coded__2)
  Demod__2(lll:hhh) =  bit__Coded__2(j);%%%%%1e4
  hhh=hhh+100; lll=lll+100;
end
%%%%TIME%%%%%
figure
plot(tt,Demod__2);
ylim([-3 3]); grid on;
xlabel('Time'); ylabel('Magnitude');
title('demodulated BPSK');
%%%%FREQ%%%%
figure
freq_demodulated_BPSK=fftshift(fft(Demod__2))*pulse_Duration_2;
plot(ff,abs(freq_demodulated_BPSK));grid on;
xlabel('Frequency'); ylabel('Magnitude');
title('freq demodulated BPSK');
%%%%%%%%%BER%%%%%%%%%
noerror = sum(bit__Coded__2~=Demod);
ber=2*noerror/100;
disp(['BER: ' num2str(ber)]);
