clear all
close all
clc
load('data(1028505).mat')
n=1:5000;
% ------------------------------- QRS detection---------------------------

m=zeros(1,length(n));

surf(CWT_db3_qrs.coefs);
title('CWT lead V4 (QRS detect) with scale factor a=[1:100]');
figure
%% finding scale factor
[hh,nn]=max(CWT_db3_qrs.coefs);
[kk,index]=max(hh);
scale_factor0=nn(index);
plot(CWT_db3_qrs.coefs(scale_factor0,:));
title(['CWT lead V4 (QRS detect) with scale factor a0=',num2str(scale_factor0)]);
figure
abs_C_qrs=abs(CWT_db3_qrs.coefs(scale_factor0,:));
plot(abs_C_qrs);
title(['absolute lead V4 (QRS detect) with scale factor a0=',num2str(scale_factor0)]);
[x,y] = hist(abs_C_qrs,9);
h=sum(x.*y)/sum(y);
for i=1:length(n)
    if abs_C_qrs(1,i) >= h
     m(1,i)=1   ;
    end
end
mQRS=m;
for i=1:4990
   
    if mQRS(1,i)== 1 && mQRS(1,i+5)== 0 && mQRS(1,i+10)== 1
        for j=1:10
      mQRS(1,i+j) = 1;
        end
    end
end

figure
stairs(mQRS);
title('mask mQRS');
axis([0 5000 0 2])
Vecg=VarName11;
Vqrs=mQRS.*Vecg';
V1ecg=Vecg'-Vqrs;
figure
subplot(3,1,1)
plot(Vecg)
title('ECG lead I');
subplot(3,1,2)
plot(Vqrs)
title('QRS detect');
subplot(3,1,3)
plot(V1ecg)
title('ECG - QRS');
%------------------------------- T detection---------------------------------
mt=zeros(1,length(n));
figure
surf(CWT_db3_T.coefs);
%% finding scale factor
[hh,nn]=max(CWT_db3_T.coefs);
[kk,index]=max(hh);
scale_factor1=nn(index);
title('CWT lead I (T detect) with scale factor a=[1:100]');
figure
plot(CWT_db3_T.coefs(scale_factor1,:));
title(['CWT lead I (T detect) with scale factor a1=',num2str(scale_factor1)]);
figure
abs_C_T=abs(CWT_db3_T.coefs(scale_factor1,:));
plot(abs_C_T);
title(['absolute lead I (T detect) with scale factor a1=',num2str(scale_factor1)]);
[x1,y1] = hist(abs_C_T,24);
h1=sum(x1.*y1)/sum(y1);

for i=1:length(n)
   
    if abs_C_T(1,i) >= h1
     mt(1,i)=1   ;
    end
end
mT=mt;
for i=1:length(n)
   
    if mT(1,i)== 1 && mT(1,i+18)== 0 && mT(1,i+36)== 1
        for j=1:36
      mT(1,i+j) = 1;
        end
    end
end
figure
plot(V1ecg)
hold on
stairs(60*mT);
axis([0 5000 -150 250])
title('mask mT');
VT=mT.*V1ecg;
V2ecg=V1ecg-VT;
figure
subplot(3,1,1)
plot(V1ecg)
title('V1ecg ');
subplot(3,1,2)
plot(VT)
title('T detect');
subplot(3,1,3)
plot(V2ecg)
title('P=V1ecg - T');
%-------------------------------- P detection-----------------------------
mp=zeros(1,length(n));
figure
surf(CWT_db3_P.coefs);
%% finding scale factor
[hh,nn]=max(CWT_db3_P.coefs);
[kk,index]=max(hh);
scale_factor2=nn(index);
title('CWT lead V4 (P detect) with scale factor a=[1:100]');
figure
plot(CWT_db3_P.coefs(scale_factor2,:));
title(['CWT lead V4 (P detect) with scale factor a2=',num2str(scale_factor2)]);
figure
abs_C_P=abs(CWT_db3_P.coefs(scale_factor2,:));
plot(abs_C_P);
title(['absolute lead V4 (P detect) with scale factor a2=',num2str(scale_factor2)]);
[x2,y2] = hist(abs_C_P,100);
h2=sum(x2.*y2)/sum(y2);

for i=1:4088
   
    if abs_C_P(1,i) >= h2
     mp(1,i)=1;
    end
end
mP=mp;

for i=1:length(n)
   
    if mP(1,i)== 1  && mP(1,i+70)== 1
        for j=1:70
      mP(1,i+j) = 1;
        end
    end
end

figure
plot(V2ecg)
hold on
stairs(20*mP);
axis([0 5000 -30 30])
title('mask mP');
VP=mP.*V2ecg;
V3ecg=V2ecg-VP;
figure
subplot(3,1,1)
plot(V2ecg)
title('V2ecg ');
subplot(3,1,2)
plot(VP)
title('P detect');
subplot(3,1,3)
plot(V3ecg)
title('V2ecg - P');

start=170;
T=Vqrs(start:start+60);
size_signal = 5000;
t=0:0.001:(size_signal);
T=T-mean(T);
for k=1:(size_signal - 61)
    s=Vqrs(k:k+60)-mean(Vqrs(k:k+60));
    rxy(k)=sum(s.*T)/sqrt(sum(s.^2)*sum(T.^2));
    if rxy(k)> 0.6
        R(k) =1;
    else R(k)=0;
    end 
end
[C,D] = findpeaks(R);
[A,B,delay]=pan_tompkin(Vqrs,250,1);
B = B(1:12);
RMSE = sqrt(((D(1)-B(1))^2+(D(2)-B(2))^2+(D(3)-B(3))^2+(D(4)-B(4))^2+(D(5)-B(5))^2+(D(6)-B(6))^2+(D(7)-B(7))^2+(D(8)-B(8))^2+(D(9)-B(9))^2+(D(10)-B(10))^2+(D(11)-B(11))^2+(D(12)-B(12))^2)./12)