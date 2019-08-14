%Code_for_Fig3_and_Fig4.m
%ULA structure
% Estimate the DOA based on the real array data from Yang Feng using 
% Root MUSIC algorithm
% +noise！！！ +SR！！！

clear all;
close all;
clc;

fc=20e6;
Monte_Carlo_time=200;
Parallel_SR_Number=100;
lambda=physconst('lightspeed') / fc;

%% 以下考虑将ULA阵列天线个数设为4

hula=phased.ULA('NumElements',4,'ElementSpacing',lambda/2);
hula.Element.FrequencyRange=[19e6 21e6];

%% 相位差0_90_180_270 -50dbm 
signal=load('data_file_30.txt');

%每次循环时改变sigma_n的值（噪声信号标准差）
for sigma_iter=1:11
    sigma_n(sigma_iter)=0.3*(sigma_iter-1);

for i=1:Monte_Carlo_time
    
  (i+(sigma_iter-1)*Monte_Carlo_time)/(3*Monte_Carlo_time*11),
    
  start_point1=unidrnd(990000);
  rxsig1=hilbert(signal(start_point1:(start_point1+999),3:6));
  start_point2=unidrnd(990000);
  rxsig2=hilbert(signal(start_point2:(start_point2+999),3:6));
  hroot=phased.RootMUSICEstimator('SensorArray',hula,...
      'OperatingFrequency',fc,'NumSignalsSource','Property',...
  'NumSignals',1,'ForwardBackwardAveraging',true);

  
  rxsig=rxsig1;
    DoA_RootMUSIC=step(hroot,rxsig);
  
  % +noise！！！ 
  Rxsig1=rxsig1+sigma_n(sigma_iter)*(randn(size(rxsig1))+1i*(randn(size(rxsig1))));
  Rxsig=Rxsig1;
  DoA_RootMUSIC2=step(hroot,Rxsig);

  
  % +SR！！！
  delta_t=0.0195;
  [row_number,column_number]=size(Rxsig);
  Rxsig2=zeros(size(Rxsig));
  for j=1:column_number
      xn_real=real(Rxsig(:,j));
      std_xn_real=std(xn_real);
      for k=1:row_number
          xn2_real(k)=xn_real(k)./std_xn_real;
      end
      x_real(1)=2.0*rand-1.0;
      for m=2:row_number
          x_real(m)=x_real(m-1)+delta_t*(x_real(m-1)-power(x_real(m-1),3)+60.0*xn2_real(m-1));
      end
      
      xn_imag=imag(Rxsig(:,j));
      std_xn_imag=std(xn_imag);
      for k=1:row_number
          xn2_imag(k)=xn_imag(k)./std_xn_imag;
      end
      x_imag(1)=2.0*rand-1.0;
      for m=2:row_number
          x_imag(m)=x_imag(m-1)+delta_t*(x_imag(m-1)-power(x_imag(m-1),3)+60.0*xn2_imag(m-1));
      end
      
      Rxsig2(:,j)=x_real+1i*x_imag;
  end
  DoA_RootMUSIC3=step(hroot,Rxsig2);
  
  %parallel SR for each receiving antenna!
  Rxsig3=zeros(size(Rxsig));
  for j=1:column_number
      xn_real=real(Rxsig(:,j));
      std_xn_real=std(xn_real);
      for k=1:row_number
          xn2_real(k)=xn_real(k)./std_xn_real;
      end
      
      xn_imag=imag(Rxsig(:,j));
      std_xn_imag=std(xn_imag);
      for k=1:row_number
          xn2_imag(k)=xn_imag(k)./std_xn_imag;
      end
      
      for p=1:Parallel_SR_Number
          x_real_PSR(p,1)=2.0*rand-1.0;
          for m=2:row_number
              x_real_PSR(p,m)=x_real_PSR(p,m-1)+delta_t*(x_real_PSR(p,m-1)-power(x_real_PSR(p,m-1),3)+60.0*xn2_real(m-1));
          end

          x_imag_PSR(p,1)=2.0*rand-1.0;
          for m=2:row_number
              x_imag_PSR(p,m)=x_imag_PSR(p,m-1)+delta_t*(x_imag_PSR(p,m-1)-power(x_imag_PSR(p,m-1),3)+60.0*xn2_imag(m-1));
          end
      
      end
      
          temp_real3=mean(x_real_PSR);
          temp_imag3=mean(x_imag_PSR);
          x_real3=temp_real3';
          x_imag3=temp_imag3';
      %end
      
      Rxsig3(:,j)=x_real3+1i*x_imag3;
  end
  DoA_RootMUSIC4=step(hroot,Rxsig3);
  
  
  Error_DoA_RootMUSIC(i)=30-DoA_RootMUSIC;
  Error_DoA_RootMUSIC2(i)=30-DoA_RootMUSIC2;
  Error_DoA_RootMUSIC3(i)=30-DoA_RootMUSIC3;
  Error_DoA_RootMUSIC4(i)=30-DoA_RootMUSIC4;
  

end


Average_Error_DoA=mean(Error_DoA_RootMUSIC),
Std_Dev_DoA=std(Error_DoA_RootMUSIC);

Average_Error_DoA2=mean(Error_DoA_RootMUSIC2),
Std_Dev_DoA2=std(Error_DoA_RootMUSIC2);

Average_Error_DoA3=mean(Error_DoA_RootMUSIC3),
Std_Dev_DoA3=std(Error_DoA_RootMUSIC3);

Average_Error_DoA4=mean(Error_DoA_RootMUSIC4),
Std_Dev_DoA4=std(Error_DoA_RootMUSIC4);



% 画曲线的数据
no_noise_no_SR_4array(sigma_iter)=Average_Error_DoA;
has_noise_no_SR_4array(sigma_iter)=Average_Error_DoA2;
has_noise_has_single_SR_4array(sigma_iter)=Average_Error_DoA3;
has_noise_has_parallel_SR_4array(sigma_iter)=Average_Error_DoA4;

end

%% 以下考虑将ULA阵列天线个数变为24

hula=phased.ULA('NumElements',24,'ElementSpacing',lambda/2);
hula.Element.FrequencyRange=[19e6 21e6];

%% 相位差0_90_180_270 -50dbm 
signal=load('data_file_30.txt');

%每次循环时改变sigma_n的值（噪声信号标准差）
for sigma_iter=1:11
    sigma_n(sigma_iter)=0.3*(sigma_iter-1);

for i=1:Monte_Carlo_time
    
  (1.0/3.0)+(i+(sigma_iter-1)*Monte_Carlo_time)/(3*Monte_Carlo_time*11),
    
  rxsig1=hilbert(signal(start_point1:(start_point1+999),3:6));
  rxsig2=hilbert(signal(start_point2:(start_point2+999),3:6));
  hroot=phased.RootMUSICEstimator('SensorArray',hula,...
      'OperatingFrequency',fc,'NumSignalsSource','Property',...
  'NumSignals',1,'ForwardBackwardAveraging',true);

  rxsig=rxsig1;
  rxsig=[rxsig1,rxsig1,rxsig1,rxsig1,rxsig1,rxsig1];
  DoA_RootMUSIC5=step(hroot,rxsig);
  
  % +noise！！！ 
  Rxsig1=rxsig1+sigma_n(sigma_iter)*(randn(size(rxsig1))+1i*(randn(size(rxsig1))));
  Rxsig=Rxsig1;
  Rxsig=[Rxsig1,Rxsig1,Rxsig1,Rxsig1,Rxsig1,Rxsig1];
  DoA_RootMUSIC6=step(hroot,Rxsig);

  
  % +SR！！！
  delta_t=0.0195;
  [row_number,column_number]=size(Rxsig);
  Rxsig2=zeros(size(Rxsig));
  for j=1:column_number
      xn_real=real(Rxsig(:,j));
      std_xn_real=std(xn_real);
      for k=1:row_number
          xn2_real(k)=xn_real(k)./std_xn_real;
      end
      x_real(1)=2.0*rand-1.0;
      for m=2:row_number
          x_real(m)=x_real(m-1)+delta_t*(x_real(m-1)-power(x_real(m-1),3)+60.0*xn2_real(m-1));
      end
      
      xn_imag=imag(Rxsig(:,j));
      std_xn_imag=std(xn_imag);
      for k=1:row_number
          xn2_imag(k)=xn_imag(k)./std_xn_imag;
      end
      x_imag(1)=2.0*rand-1.0;
      for m=2:row_number
          x_imag(m)=x_imag(m-1)+delta_t*(x_imag(m-1)-power(x_imag(m-1),3)+60.0*xn2_imag(m-1));
      end
      
      Rxsig2(:,j)=x_real+1i*x_imag;
  end
  DoA_RootMUSIC7=step(hroot,Rxsig2);
  
  %parallel SR for each receiving antenna!
  Rxsig3=zeros(size(Rxsig));
  for j=1:column_number
      xn_real=real(Rxsig(:,j));
      std_xn_real=std(xn_real);
      for k=1:row_number
          xn2_real(k)=xn_real(k)./std_xn_real;
      end
      
      xn_imag=imag(Rxsig(:,j));
      std_xn_imag=std(xn_imag);
      for k=1:row_number
          xn2_imag(k)=xn_imag(k)./std_xn_imag;
      end
      
      for p=1:Parallel_SR_Number
          x_real_PSR(p,1)=2.0*rand-1.0;
          for m=2:row_number
              x_real_PSR(p,m)=x_real_PSR(p,m-1)+delta_t*(x_real_PSR(p,m-1)-power(x_real_PSR(p,m-1),3)+60.0*xn2_real(m-1));
          end

          x_imag_PSR(p,1)=2.0*rand-1.0;
          for m=2:row_number
              x_imag_PSR(p,m)=x_imag_PSR(p,m-1)+delta_t*(x_imag_PSR(p,m-1)-power(x_imag_PSR(p,m-1),3)+60.0*xn2_imag(m-1));
          end
      
      end
      
          temp_real3=mean(x_real_PSR);
          temp_imag3=mean(x_imag_PSR);
          x_real3=temp_real3';
          x_imag3=temp_imag3';
      %end
      
      Rxsig3(:,j)=x_real3+1i*x_imag3;
  end
  DoA_RootMUSIC8=step(hroot,Rxsig3);
  
  
  Error_DoA_RootMUSIC5(i)=30-DoA_RootMUSIC5;
  Error_DoA_RootMUSIC6(i)=30-DoA_RootMUSIC6;
  Error_DoA_RootMUSIC7(i)=30-DoA_RootMUSIC7;
  Error_DoA_RootMUSIC8(i)=30-DoA_RootMUSIC8;
  

end


Average_Error_DoA5=mean(Error_DoA_RootMUSIC5),
Std_Dev_DoA5=std(Error_DoA_RootMUSIC5);

Average_Error_DoA6=mean(Error_DoA_RootMUSIC6),
Std_Dev_DoA6=std(Error_DoA_RootMUSIC6);

Average_Error_DoA7=mean(Error_DoA_RootMUSIC7),
Std_Dev_DoA7=std(Error_DoA_RootMUSIC7);

Average_Error_DoA8=mean(Error_DoA_RootMUSIC8),
Std_Dev_DoA8=std(Error_DoA_RootMUSIC8);

%distance error in 30km case
30000*abs(sin(Average_Error_DoA5*pi/180));
30000*abs(sin(Average_Error_DoA6*pi/180));
30000*abs(sin(Average_Error_DoA7*pi/180));
30000*abs(sin(Average_Error_DoA8*pi/180));


% 画曲线的数据
no_noise_no_SR_24array(sigma_iter)=Average_Error_DoA5;
has_noise_no_SR_24array(sigma_iter)=Average_Error_DoA6;
has_noise_has_single_SR_24array(sigma_iter)=Average_Error_DoA7;
has_noise_has_parallel_SR_24array(sigma_iter)=Average_Error_DoA8;

end

%% 以下考虑将ULA阵列天线个数变为48

hula=phased.ULA('NumElements',48,'ElementSpacing',lambda/2);
hula.Element.FrequencyRange=[19e6 21e6];

%% 相位差0_90_180_270 -50dbm 
signal=load('data_file_30.txt');

%每次循环时改变sigma_n的值（噪声信号标准差）
for sigma_iter=1:11
    sigma_n(sigma_iter)=0.3*(sigma_iter-1);

for i=1:Monte_Carlo_time
    
  (2.0/3.0)+(i+(sigma_iter-1)*Monte_Carlo_time)/(3*Monte_Carlo_time*11),
    
  rxsig1=hilbert(signal(start_point1:(start_point1+999),3:6));
  rxsig2=hilbert(signal(start_point2:(start_point2+999),3:6));
  hroot=phased.RootMUSICEstimator('SensorArray',hula,...
      'OperatingFrequency',fc,'NumSignalsSource','Property',...
  'NumSignals',1,'ForwardBackwardAveraging',true);

  rxsig=rxsig1;
  rxsig=[rxsig1,rxsig1,rxsig1,rxsig1,rxsig1,rxsig1,rxsig1,rxsig1,rxsig1,rxsig1,rxsig1,rxsig1];
  DoA_RootMUSIC9=step(hroot,rxsig);
  
  % +noise！！！ 
  Rxsig1=rxsig1+sigma_n(sigma_iter)*(randn(size(rxsig1))+1i*(randn(size(rxsig1))));
  Rxsig=Rxsig1;
  Rxsig=[Rxsig1,Rxsig1,Rxsig1,Rxsig1,Rxsig1,Rxsig1,Rxsig1,Rxsig1,Rxsig1,Rxsig1,Rxsig1,Rxsig1];
  DoA_RootMUSIC10=step(hroot,Rxsig);

  
  % +SR！！！
  delta_t=0.0195;
  [row_number,column_number]=size(Rxsig);
  Rxsig2=zeros(size(Rxsig));
  for j=1:column_number
      xn_real=real(Rxsig(:,j));
      std_xn_real=std(xn_real);
      for k=1:row_number
          xn2_real(k)=xn_real(k)./std_xn_real;
      end
      x_real(1)=2.0*rand-1.0;
      for m=2:row_number
          x_real(m)=x_real(m-1)+delta_t*(x_real(m-1)-power(x_real(m-1),3)+60.0*xn2_real(m-1));
      end
      
      xn_imag=imag(Rxsig(:,j));
      std_xn_imag=std(xn_imag);
      for k=1:row_number
          xn2_imag(k)=xn_imag(k)./std_xn_imag;
      end
      x_imag(1)=2.0*rand-1.0;
      for m=2:row_number
          x_imag(m)=x_imag(m-1)+delta_t*(x_imag(m-1)-power(x_imag(m-1),3)+60.0*xn2_imag(m-1));
      end
      
      Rxsig2(:,j)=x_real+1i*x_imag;
  end
  DoA_RootMUSIC11=step(hroot,Rxsig2);
  
  %parallel SR for each receiving antenna!
  Rxsig3=zeros(size(Rxsig));
  for j=1:column_number
      xn_real=real(Rxsig(:,j));
      std_xn_real=std(xn_real);
      for k=1:row_number
          xn2_real(k)=xn_real(k)./std_xn_real;
      end
      
      xn_imag=imag(Rxsig(:,j));
      std_xn_imag=std(xn_imag);
      for k=1:row_number
          xn2_imag(k)=xn_imag(k)./std_xn_imag;
      end
      
      for p=1:Parallel_SR_Number
          x_real_PSR(p,1)=2.0*rand-1.0;
          for m=2:row_number
              x_real_PSR(p,m)=x_real_PSR(p,m-1)+delta_t*(x_real_PSR(p,m-1)-power(x_real_PSR(p,m-1),3)+60.0*xn2_real(m-1));
          end

          x_imag_PSR(p,1)=2.0*rand-1.0;
          for m=2:row_number
              x_imag_PSR(p,m)=x_imag_PSR(p,m-1)+delta_t*(x_imag_PSR(p,m-1)-power(x_imag_PSR(p,m-1),3)+60.0*xn2_imag(m-1));
          end
      
      end
      
          temp_real3=mean(x_real_PSR);
          temp_imag3=mean(x_imag_PSR);
          x_real3=temp_real3';
          x_imag3=temp_imag3';
      %end
      
      Rxsig3(:,j)=x_real3+1i*x_imag3;
  end
  DoA_RootMUSIC12=step(hroot,Rxsig3);
  
  
  Error_DoA_RootMUSIC9(i)=30-DoA_RootMUSIC9;
  Error_DoA_RootMUSIC10(i)=30-DoA_RootMUSIC10;
  Error_DoA_RootMUSIC11(i)=30-DoA_RootMUSIC11;
  Error_DoA_RootMUSIC12(i)=30-DoA_RootMUSIC12;
  

end


Average_Error_DoA9=mean(Error_DoA_RootMUSIC9),
Std_Dev_DoA9=std(Error_DoA_RootMUSIC9);

Average_Error_DoA10=mean(Error_DoA_RootMUSIC10),
Std_Dev_DoA10=std(Error_DoA_RootMUSIC10);

Average_Error_DoA11=mean(Error_DoA_RootMUSIC11),
Std_Dev_DoA11=std(Error_DoA_RootMUSIC11);

Average_Error_DoA12=mean(Error_DoA_RootMUSIC12),
Std_Dev_DoA12=std(Error_DoA_RootMUSIC12);


% 画曲线的数据
no_noise_no_SR_48array(sigma_iter)=Average_Error_DoA9;
has_noise_no_SR_48array(sigma_iter)=Average_Error_DoA10;
has_noise_has_single_SR_48array(sigma_iter)=Average_Error_DoA11;
has_noise_has_parallel_SR_48array(sigma_iter)=Average_Error_DoA12;

end


%plot figures
figure(1);
plot(sigma_n,abs(no_noise_no_SR_4array),'bs-','LineWidth',2,'MarkerSize',10);
hold on;
plot(sigma_n,abs(has_noise_no_SR_4array),'g^-','LineWidth',2,'MarkerSize',10);
hold on;
plot(sigma_n,abs(has_noise_has_single_SR_4array),'r*-','LineWidth',2,'MarkerSize',10);
hold on;
plot(sigma_n,abs(has_noise_has_parallel_SR_4array),'ko-','LineWidth',2,'MarkerSize',10);
hold on;
plot(sigma_n,abs(no_noise_no_SR_24array),'bs--','LineWidth',2,'MarkerSize',10);
hold on;
plot(sigma_n,abs(has_noise_no_SR_24array),'g^--','LineWidth',2,'MarkerSize',10);
hold on;
plot(sigma_n,abs(has_noise_has_single_SR_24array),'r*--','LineWidth',2,'MarkerSize',10);
hold on;
plot(sigma_n,abs(has_noise_has_parallel_SR_24array),'ko--','LineWidth',2,'MarkerSize',10);
hold on;
plot(sigma_n,abs(no_noise_no_SR_48array),'bs:','LineWidth',2,'MarkerSize',10);
hold on;
plot(sigma_n,abs(has_noise_no_SR_48array),'g^:','LineWidth',2,'MarkerSize',10);
hold on;
plot(sigma_n,abs(has_noise_has_single_SR_48array),'r*:','LineWidth',2,'MarkerSize',10);
hold on;
plot(sigma_n,abs(has_noise_has_parallel_SR_48array),'ko:','LineWidth',2,'MarkerSize',10);
hold on;
title('multiple-antenna ULA');
legend('no noise no SR (4-antenna)','has noise no SR (4-antenna)','has noise has single SR (4-antenna)','has noise has parallel SR (4-antenna)', ...
       'no noise no SR (24-antenna)','has noise no SR (24-antenna)','has noise has single SR (24-antenna)','has noise has parallel SR (24-antenna)', ...
       'no noise no SR (48-antenna)','has noise no SR (48-antenna)','has noise has single SR (48-antenna)','has noise has parallel SR (48-antenna)','location','northeastoutside');
xlabel('Noise Standard Deviation \sigma_n');
ylabel('DoA Estimation Error (degree)');
grid on;
hold on;
