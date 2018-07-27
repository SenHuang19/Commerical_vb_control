% clear all
% clc

function [mm_now,T_next,P_signal_now,P_now]=Commercial_Control(T_now,T_out,time)


if isrow(T_now)
    T_now=T_now';
end

% A=xlsread('3820.xlsx','sheet1','A24:K40');
% B=xlsread('3820.xlsx','sheet1','O24:S40');
% m_min=xlsread('3820.xlsx','sheet1','Y24:Y40');

A=[1	2	3	18	0.52793091	0.610067	1.00E-02	0	0.00499903	0.373154	0.005
2	1	0	18	0.56617236	0.905277	0.01	2.81E-02	0.170145	8.02E-02	0
3	1	7	18	0.21547656	0.902911	0.01	0	0.47428	0.0620888	0.005
4	3	5	18	0.03989817	0.965	1.00E-02	1.87E-15	0.412673	4.92E-16	0.005
5	3	4	18	0.04085172	0.5	1.00E-02	0.112961	0.451178	3.64E-17	0.47
6	4	5	18	0.03714633	0.929514	0.0108144	1	1	0.0366324	0.005
7	8	3	18	0.26322426	0.5	0.01	0.0305154	0.424941	0.174063	0.295937
8	7	3	18	0.27065025	0.654008	1.00E-02	0.0362279	0.449949	2.62E-16	0.315992
9	10	0	18	0.74943531	0.97	0.01	0.0369616	0.557404	2.50E-16	0
10	9	0	18	1.287	0.510189	0.0349414	0.0314956	0.107384	0.459811	0
11	9	12	18	0.13857129	0.965	0.01	0.0830073	0.479562	0	5.00E-03
12	9	11	18	0.05653206	0.5	0.01	0.338011	0.447386	3.07E-16	4.70E-01
13	9	10	18	0.12399543	0.965	0.01	0.24761	0.587417	0	0.005
14	9	0	18	0.15456519	0.872358	0.0244433	3.73E-15	2.34E-15	0.0976422	0
15	16	10	18	0.85861269	0.847754	0.0223683	0	0	8.99E-15	1.22E-01
16	15	10	18	0.16213509	0.725983	0.0227574	4.15E-15	0	0	0.244017
17	10	0	18	0.17521218	0.97	0.01	6.24E-17	0.443958	0.00E+00	0];
B=[2	3	0	0	0
1	0	0	0	0
1	4	5	7	8
5	6	0	0	0
4	6	0	0	0
0	0	0	0	0
3	8	0	0	0
7	0	0	0	0
10	11	12	13	14
9	13	15	16	17
12	0	0	0	0
11	0	0	0	0
0	0	0	0	0
0	0	0	0	0
16	0	0	0	0
15	0	0	0	0
0	0	0	0	0
];

Zone=A(:,1);
A_Z1=A(:,2);
A_Z2=A(:,3);
Tc=13;

m_max=[0.52793091;0.56617236;0.21547656;0.03989817;0.04085172;0.03714633;0.26322426;0.27065025;0.74943531;1.287;0.13857129;0.05653206;0.12399543;0.15456519;0.85861269;0.16213509;0.17521218;];
m_min=m_max*0.3;
a1=A(:,6);
a2=A(:,7);
a3=A(:,8);
a4=A(:,9);
a5=A(:,10);
a6=A(:,11);

A_Z1_Opt=A_Z1;
for i=1:1:17;
    if A_Z2(i)==0;
        A_Z2_Opt(i,1)=1;
    else
        A_Z2_Opt(i,1)=A_Z2(i);
    end
end

%Fan power factors
% RR=[0.0308354498677761,0.241847658860779,-1.42193698864903,2.99284650687545];
RR=[0.0264629913582345,-0.200260041627581,0.763315238378467,-0.687871384050606];
c1=RR(1); %1;
c2=RR(2); %1;
c3=RR(3); %1;
c4=RR(4); %1;

%Chiller power factors
lambda=0.9;
fra=0.7;
Tc=13;
COP=2.5;
cp=6.5;
% cp=2;
% Prate=5.6;
T_set=22.5;
deadband=0.625;

% T=22.5*ones(17,1);
%status=[1;0;0;0;0;1;0;0;1;0;0;0;1;0;0;1;1];

% P_signal(1)=0;
% T_wp=22.5*ones(17,1);
% % P_signal=45*ones(1,L);
% % load('AGC.mat')
%load('Ramping.mat')
% % load('Loadshift.mat');
% Signal=15*Ramping(1:95)';
load('P_signal.mat')
% P_signal=P+Signal;
P_signal=P;
m_target=0;
mm_max=zeros(17,1);
T_wp=zeros(17,1);

% for i=2:1:L

%%
%Control
    T_sort2(:,1)=T_now;
    T_sort2(:,2)=1:17;
    T_sorted2=sortrows(T_sort2);

% if i==2
% a=c1;b=c2;c=(c3+cp*(fra*(sum(T_wp(:,i-1))/17)+(1-fra)*T_out(i-1)-Tc)/(lambda*COP));d=c4-P_signal(i-1);
% 
% else
  a=c1;b=c2;c=(c3+cp*(fra*(sum(T_sorted2(7:17,1))/11)+(1-fra)*T_out-Tc)/(lambda*COP));d=c4-P_signal(time);
%     a=c1;b=c2;c=(c3+cp*(fra*(sum(T_wp(3:17,i-1))/15)+(1-fra)*T_out(i-1)-Tc)/(lambda*COP));d=c4-P_signal(i-1);
  
%   Tm2_wp(i-1)=fra*(sum(T_sorted(17,1))/1)+(1-fra)*T_out(i-1);
% a=c1;b=c2;c=(c3+cp*(fra*(Tm_wp(i-2))+(1-fra)*T_out(i-1)-Tc)/(lambda*COP));d=c4-P_signal(i-1);   
% % a=c1;b=c2;c=(c3+cp*(fra*(Tm_wp(i-2)+0.3*(T_out(i-1)-T_out(i-2)))+(1-fra)*T_out(i-1)-Tc)/(lambda*COP));d=c4-P_signal(i-1);    
% end

Q=((2*b^3-9*a*b*c+27*a*a*d)^2-4*(b^2-3*a*c)^3)^(1/2);
C=(0.5*(Q+2*b^3-9*a*b*c+27*a*a*d))^(1/3);

m_target=-b/(3*a)-C/(3*a)-(b^2-3*a*c)/(3*a*C);
m_ref=m_target;
flag=17;
% x11(i-1)=-b/(3*a)+C*(1+sqrt(-3))/(6*a)+(1-sqrt(-3))*(b^2-3*a*c)/(6*a*C);
% x111(i-1)=-b/(3*a)+C*(1-sqrt(-3))/(6*a)+(1+sqrt(-3))*(b^2-3*a*c)/(6*a*C);

    for j=1:1:17
       if A_Z2(j,1)==0
       mm_max(j,1)=(a1(j,1)*T_now(j,1)+a2(j,1)*T_out+a4(j,1)+a5(j,1)*T_now(A_Z1(j,1),1)-(T_set-deadband))/(a3(j,1)*(T_now(j,1)-Tc));    
       else
       mm_max(j,1)=(a1(j,1)*T_now(j,1)+a2(j,1)*T_out+a4(j,1)+a5(j,1)*T_now(A_Z1(j,1),1)+a6(j,1)*T_now(A_Z2(j,1),1)-(T_set-deadband))/(a3(j,1)*(T_now(j,1)-Tc));
       end
    end
    
    
    T_sort(:,1)=zeros(17,1);
    T_sort(:,2)=T_now(:,1);
    T_sort(:,3)=1:17;
    
    for ff2=1:1:17
           if A_Z2(T_sort(ff2,3),1)==0
   TT_mid_est=a1(T_sort(ff2,3),1)*T_now(T_sort(ff2,3),1)+a2(T_sort(ff2,3),1)*T_out+a3(T_sort(ff2,3),1)*m_min(T_sort(ff2,3),1)*(Tc-T_now(T_sort(ff2,3),1))+a4(T_sort(ff2,3),1)+a5(T_sort(ff2,3),1)*T_now(A_Z1(T_sort(ff2,3),1),1);
       else
   TT_mid_est=a1(T_sort(ff2,3),1)*T_now(T_sort(ff2,3),1)+a2(T_sort(ff2,3),1)*T_out+a3(T_sort(ff2,3),1)*m_min(T_sort(ff2,3),1)*(Tc-T_now(T_sort(ff2,3),1))+a4(T_sort(ff2,3),1)+a5(T_sort(ff2,3),1)*T_now(A_Z1(T_sort(ff2,3),1),1)+a6(T_sort(ff2,3),1)*T_now(A_Z2(T_sort(ff2,3),1),1);        
       end
     
       if TT_mid_est>T_set+deadband
           T_sort(ff2,1)=1;
       else
           T_sort(ff2,1)=0;
       end
    end  
     
    for j=1:1:17;
    mm_max_mid(j,1)=min(mm_max(j,1),m_max(j,1));
    if mm_max_mid(j,1)<m_min(j,1)
        mm_max_mid(j,1)=m_min(j,1);
    end
    end
    mm_max_mid_rec(:,1)=mm_max_mid;
    
    T_sorted=sortrows(T_sort);
    
    m_ref=m_target-sum(m_min);
    mm(:,1)=m_min;
    if m_ref<0
        m_ref=0;
    end
  
    %%
    %Arrange m_ref
    while m_ref~=0
        mm_mid=min(mm_max(T_sorted(flag,3),1),m_max(T_sorted(flag,3),1));
         flag2=0;
%       if T_wp(T_sorted(flag,2),i-1)<T_set-deadband;
   for ff=1:1:length(find(B(T_sorted(flag,3),:)>0)) 
     TT_mid=T_now(B(T_sorted(flag,3),ff),1);
     if TT_mid<T_set-deadband
         flag2=1;
     end
   end
   
       
   if flag2==0||T_sorted(flag,1)==1
       
        if mm_mid<=m_min(T_sorted(flag,3),1);
        mm(T_sorted(flag,3),1)=m_min(T_sorted(flag,3),1); 
        m_ref=m_ref;
        elseif mm_mid>m_ref+m_min(T_sorted(flag,3),1)
        mm(T_sorted(flag,3),1)=m_ref+m_min(T_sorted(flag,3),1);
        m_ref=0;
        else
        mm(T_sorted(flag,3),1)=mm_mid;
        m_ref=m_ref-mm_mid+m_min(T_sorted(flag,3),1);
        end
        flag=flag-1;
    
   else
      flag=flag-1; 
   end
        if flag==0
            m_ref=0;
        end
   
   end
   %%
   
   %%
   %Get Temperature
    
   for j=1:1:17
       if A_Z2(j,1)==0
   T_wp(j,1)=a1(j,1)*T_now(j,1)+a2(j,1)*T_out+a3(j,1)*mm(j,1)*(Tc-T_now(j,1))+a4(j,1)+a5(j,1)*T_now(A_Z1(j,1),1);
       else
   T_wp(j,1)=a1(j,1)*T_now(j,1)+a2(j,1)*T_out+a3(j,1)*mm(j,1)*(Tc-T_now(j,1))+a4(j,1)+a5(j,1)*T_now(A_Z1(j,1),1)+a6(j,1)*T_now(A_Z2(j,1),1);        
       end
   end
     

%    flag3(:,i)=zeros(17,1);
%    
%    for flag4=1:1:17
%    for ff=1:1:length(find(B(T_sorted2(flag4,2),:)>0)) 
%      TT_mid2=T_wp(B(T_sorted2(flag4,2),ff),i);
%      if TT_mid2<T_set-deadband
%          flag3(T_sorted2(flag4,2),i)=1;
%      end
%    end
%    end
    
       %%
     mm_sum(1)=sum(mm(:,1));
   P_fan_wp=c1*(mm_sum(1))^3+c2*(mm_sum(1))^2+c3*(mm_sum(1))+c4;
   E_mix_wp=0;
%    E_mix_wp2=0;
   for j=1:1:17
       E_mix_wp=E_mix_wp+mm(j,1)*T_wp(j,1);
%        E_mix_wp2=E_mix_wp2(i-1)+mm(j,i-1)*T_wp(j,i);
   end
   
   Tm_wp=fra*(E_mix_wp/mm_sum(1))+(1-fra)*T_out;
%    Tm_wp2(i-1)=fra*(E_mix_wp2(i-1)/mm_sum(i-1))+(1-fra)*T_out(i);
   
   
%    Tm2(tt)=fra*(sum(T_wp(:,tt))/17)+(1-fra)*T_out(tt);
   
   P_chiller_wp=cp*(mm_sum(1))*(Tm_wp(1)-Tc)/(lambda*COP);
   
   P_wp=P_fan_wp+P_chiller_wp;
   
   T_next=T_wp;
   mm_now=mm;
   P_signal_now=P_signal(time);
   P_now=P_wp;
%    Midd=Tm_wp;
     %Mid=P_fan_wp;
%    P_now=P_chiller_wp;
   %%
    
% end


%  plot(m_sum);hold on;plot(x1)
% plot(P_signal(2:end));hold on;plot(P_wp(2:end))
% plot(P_signal(2:end)-P_wp(2:end))
% plot(P(2:end));hold on;plot(P_wp(2:end))
% for i=1:1:17;
%     figure(i)
%     set(gcf,'DefaultAxesFontSize',16);
%     plot(T(i,:),'b','Linewidth',2);
%     hold on
%     plot(T_wp(i,:),'r','Linewidth',2);
%     plot((T_set-deadband)*ones(1,96),'Linewidth',2);
%     plot((T_set+deadband)*ones(1,96),'Linewidth',2);
%     legend('Without control','With control','Location','Best');
%     xticks([0 24 48 72 96])
% xticklabels({'0','6','12','18','24'})
% xlabel(['Time [hour] Zone' num2str(i)])
% ylabel('Temperature')
% xlim([0 96]);
% end
% plot(P)
% hold on
% plot(P_wp(2:end))

% load('Error1.mat')
% figure(1);
% set(gcf,'DefaultAxesFontSize',16);
% plot(Signal,'b','Linewidth',2);hold on;
% plot(P_wp-P,'r','Linewidth',2)
% legend('Target Signal','Tracked Signal','Location','Best');
% xticks([0 24 48 72 96])
% xticklabels({'0','6','12','18','24'})
% xlabel('Time [hour]')
% ylabel('Power [kW]')
% xlim([0 96]);
% 
% figure(2);
% set(gcf,'DefaultAxesFontSize',16);
% % plot(Error1);hold on;
% plot(Signal-P_wp+P,'b','Linewidth',2)
% legend('Signal Error','Location','Best');
% xticks([0 24 48 72 96])
% xticklabels({'0','6','12','18','24'})
% xlabel('Time [hour]')
% ylabel('Power [kW]')
% xlim([0 96]);

% maa=mm-m_min.*ones(17,95);
% for i=1:1:95;maa_z(i)=length(find(maa(:,i)~=0));end
% plotyy(2:95,P_signal(2:end)-P_wp(2:end),1:95,maa_z)

% figure(1);plot(Tm);hold on;plot(Tm_wp);
% figure(2);plot(m_sum);hold on;plot(mm_sum)
% figure(3);plot(P);hold on;plot(P_wp)
% figure(4);plot(P_chiller);hold on;plot(P_chiller_wp)
% figure(5);plot(P_fan);hold on;plot(P_fan_wp)

% figure(1)
% set(gcf,'DefaultAxesFontSize',16)
% plot(Signal,'b','Linewidth',2);
% legend('Ramping Signal 2','Location','Best');
% xticks([0 24 48 72 96])
% xticklabels({'0','6','12','18','24'})
% xlabel('Time [hour]')
% ylabel('Signal')
% xlim([0 96]);

% figure(2);
% set(gcf,'DefaultAxesFontSize',16);
% fill([1:95 fliplr(1:95)],[P_wp_up,fliplr(P_wp_low)] ,'r');alpha(0.25);hold on;
% plot(P_wp,'-bo','Linewidth',2,'MarkerFaceColor','b');hold on;plot(P_signal,'-ro','Linewidth',2,'MarkerFaceColor','r');plot(P_wp_base,'-ko','Linewidth',2,'MarkerFaceColor','k');hold on;plot(P_wp_up,'Linewidth',2);hold on;plot(P_wp_low,'Linewidth',2)
% legend('Power Bound','Actual Power','Signal Power','Power Baseline','Power Upper Bound','Power Lower Bound','Location','Best');
% xticks([0 24 48 72 96])
% xticklabels({'0','6','12','18','24'})
% xlabel('Time [hour]')
% ylabel('Power [kW]')
% xlim([0 96]);
% 
Z=12;
% figure(1)
% set(gcf,'DefaultAxesFontSize',16);
% fill([1:96 fliplr(1:96)],[(T_set-deadband)*ones(1,96),fliplr((T_set+deadband)*ones(1,96))] ,'r');alpha(0.25);hold on;
% plot(T_wp(Z,:),'-ro','Linewidth',2,'MarkerFaceColor','r');hold on;
% plot((T_set-deadband)*ones(1,96),'Linewidth',2);
% plot((T_set+deadband)*ones(1,96),'Linewidth',2);
% legend('Temperature Bound','Temperature')
% xticks([0 24 48 72 96]);xticklabels({'0','6','12','18','24'})1/22
% xlabel(['Time [hour] Zone 12'])
% ylabel('Temperature')
% xlim([0 96]);


% figure(1)
% set(gcf,'DefaultAxesFontSize',16);
% fill([1:95 fliplr(1:95)],[(m_min(Z))*ones(1,95),fliplr((m_max(Z))*ones(1,95))] ,'r');alpha(0.25);hold on;
% plot(mm(Z,:),'-ro','Linewidth',2,'MarkerFaceColor','r');hold on;
% plot((m_min(Z))*ones(1,95),'Linewidth',2);
% plot((m_max(Z))*ones(1,95),'Linewidth',2);
% legend('Airflow Bound','Airflow')
% xticks([0 24 48 72 96]);xticklabels({'0','6','12','18','24'})
% xlabel(['Time [hour] Zone 12'])
% ylabel('Airflow Rate')

% for Z=1:1:17;figure(Z)
% set(gcf,'DefaultAxesFontSize',16);
% fill([1:95 fliplr(1:95)],[(m_min(Z))*ones(1,95),fliplr((mm_max_mid_rec(Z,:)))] ,'r');alpha(0.25);hold on;
% plot(mm(Z,:),'-ro','Linewidth',2,'MarkerFaceColor','r');hold on;
% plot((m_min(Z))*ones(1,95),'Linewidth',2);
% plot((mm_max_mid_rec(Z,:)),'Linewidth',2);
% legend('Airflow Bound','Airflow')
% xticks([0 24 48 72 96]);xticklabels({'0','6','12','18','24'})
% xlabel(['Time [hour] Zone 12'])
% ylabel('Airflow Rate');end