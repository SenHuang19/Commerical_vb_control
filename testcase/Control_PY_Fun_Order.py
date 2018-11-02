# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 19:50:27 2018

@author: jwang49
"""

from __future__ import division
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from sklearn import datasets, linear_model
import math
import json
import sys

def Control_Commercial(T_now,T_out,time):
    #
    
###Input from EnergyPlus
    #T_now=[21.5]*17
    #T_out=26
    #time=30
    
    ###
    #'''
    df22 = pd.DataFrame(T_now)
    df22['Tt_now']=T_now
    df22.to_csv("T_now.csv")
    
    temps2=''
    for i in range(17):
        temps2=temps2+str(T_now[i])+','
    
    f=open('T_now2.csv','a')
    f.writelines(temps2+'\n')
    f.close()
    
    df23 = pd.DataFrame()
    df23['Tt_out']=ones(1)
    df23['Tt_out'][0]=T_out
    df23.to_csv("T_out.csv")
    
    df24 = pd.DataFrame()
    df24['Tt_time']=ones(1)
    df24['Tt_time'][0]=time
    df24.to_csv("Time.csv")
    
    if time==21:
        Tm=15.5
    elif time==22:     
        #Tm_pd=array(pd.read_csv('MixAir.csv',header=None))
        #Tm_end=Tm_pd[-2]
        #Tm=Tm_end[0]
        Tm=21
    elif time==23:     
        Tm=21             
    elif time>23:
        Tm_pd=array(pd.read_csv('MixAir.csv',header=None))
        Tm_end=Tm_pd[-1]
        Tm=Tm_end[0]
    else:
        Tm=0
   
    #'''
    #Read parameters for thermal model
    tab=pd.read_csv('ParameterT.csv')
    tab2=pd.read_csv('AdZone.csv')
    #Read base power and normalized power signal
    tab3=pd.read_csv('ParameterP.csv')
    
    Zone=tab['Zone']
    A_z1=tab['Ad1']
    A_z2=tab['Ad2']
    a1=tab['a1']
    a2=tab['a2']
    a3=tab['a3']
    a4=tab['a4']
    a5=tab['a5']
    a6=tab['a6']
    Tc=13
    
    #Maximum airflow volume for each thermal zone
    m_max=pd.DataFrame()
    m_max['m_max']=[0.943450560000000,0.619580520000000,0.335658960000000,0.0816578100000000,0.0826897500000000,0.0366034500000000,0.411887970000000,0.421921890000000,0.939084120000000,1.39230000000000,0.244164960000000,0.0996161400000000,0.212036760000000,0.158256540000000,0.875851470000000,0.166255830000000,0.221093730000000]
    m_max=m_max['m_max']
    #Minimum airflow volume for each thermal zone
    m_min=m_max*0.3
    #Step size when changing airflow volume for each thermal zone
    m_delta=(m_max-m_min)/10
    
    #Parameter for fan model
    Rr=pd.DataFrame()
    Rr['FanCo']=[0.00746555686277102,-0.00921145141615396,0.130563408784124,0.00157993820518559]
    Rr=Rr['FanCo']
    c1=Rr[0]
    c2=Rr[1]
    c3=Rr[2]
    c4=Rr[3]
    
    #Parameters for chiller
    lam=0.9;
    fra=0.7;
    Tc=13;
    COP=6.16;
    cp=1.006;
    #Temperature setpoint and comfortable range
    T_set=21.1;
    deadband=0.8;
    P=tab3['Pbase']
    Ramping=tab3['Ramping']
    m_out=tab3['m_out']
    
    ## Input m_tar_rec.mat and mm_last.mat here, line140-143 from Commercial_control_method2.m
    
    #Get the operation point of previous time step
    if time>1:
        m_information1=pd.read_csv('m_tar_rec.csv')
        m_information2=pd.read_csv('mm_last.csv')
        m_tar_rec=m_information1['m_tar_rec']
        mm_last=m_information2['mm_last']
    else:
        mm_last=m_min-m_min
        
    #Ramping[1:20]=zeros(20)
    Zeroele=pd.DataFrame()
    Zeroele['1']=zeros(20)
    Zeroele=Zeroele['1']
    #Set the signal to be 0 when the office is close
    Ramping[0:20]=Zeroele
    Zeroele=pd.DataFrame()
    Zeroele['1']=zeros(28)
    Zeroele=Zeroele['1']
    Ramping[68:96]=Zeroele
    
    #Get the power which HVAC need to follow, base power plus signal power
    P_signal=P
    #Initial variables
    m_target=0;
    mm_max=[0]*17
    mm_toset=[0]*17
    mm_min=[0]*17
    mm_ref_up=[0]*17
    mm_ref_down=[0]*17
    T_wp=[0]*17
    
    #Control will work when the target power is not 0
    if P_signal[time-1]!=0:
        #T_sort2=T_now
        #Calculate the temperature of returned air
        #Tm_est=sum(m_min*T_now)/sum(m_min)
        Tm_est=Tm
        a=c1;b=c2;c=(c3+cp*((Tm_est)-Tc)/(lam*COP));
        d=c4-P_signal[time-1]
        #+cp*m_out[time-1]*(T_out-Tm_est)/(lam*COP);
        
        Q=((2*b**3-9*a*b*c+27*a*a*d)**2-4*(b**2-3*a*c)**3)**(1/2);
        C=(0.5*(Q+2*b**3-9*a*b*c+27*a*a*d))**(1/3);
        #Base on thermal model, calculate the target total airflow volume
        m_target=-b/(3*a)-C/(3*a)-(b**2-3*a*c)/(3*a*C);
        
        #For each thermal zone, calculate the airflow volume that will change its temperature to lower bound
        for j in range(0,17):
            if A_z2[j]==0:
                mm_max[j]=(a1[j]*T_now[j]+a2[j]*T_out+a4[j]+a5[j]*T_now[A_z1[j]-1]-(T_set-deadband))/(a3[j]*(T_now[j]-Tc));
            else:
                mm_max[j]=(a1[j]*T_now[j]+a2[j]*T_out+a4[j]+a5[j]*T_now[A_z1[j]-1]+a6[j]*T_now[A_z2[j]-1]-(T_set-deadband))/(a3[j]*(T_now[j]-Tc));
        
        #For each thermal zone, calculate the airflow volume that will change its temperature to setpoint               
        for j in range(0,17):
            if A_z2[j]==0:
                mm_toset[j]=(a1[j]*T_now[j]+a2[j]*T_out+a4[j]+a5[j]*T_now[A_z1[j]-1]-(T_set))/(a3[j]*(T_now[j]-Tc));
            else:
                mm_toset[j]=(a1[j]*T_now[j]+a2[j]*T_out+a4[j]+a5[j]*T_now[A_z1[j]-1]+a6[j]*T_now[A_z2[j]-1]-(T_set))/(a3[j]*(T_now[j]-Tc));

        #For each thermal zone, calculate the airflow volume that will change its temperature to upper bound                
        for j in range(0,17):
            if A_z2[j]==0:
                mm_min[j]=(a1[j]*T_now[j]+a2[j]*T_out+a4[j]+a5[j]*T_now[A_z1[j]-1]-(T_set+deadband))/(a3[j]*(T_now[j]-Tc));
            else:
                mm_min[j]=(a1[j]*T_now[j]+a2[j]*T_out+a4[j]+a5[j]*T_now[A_z1[j]-1]+a6[j]*T_now[A_z2[j]-1]-(T_set+deadband))/(a3[j]*(T_now[j]-Tc));
        
        for j in range(0,17):
            if A_z2[j]==0:
                mm_ref_up[j]=(a1[j]*T_now[j]+a2[j]*T_out+a4[j]+a5[j]*T_now[A_z1[j]-1]-(T_set-0.7*deadband))/(a3[j]*(T_now[j]-Tc));
            else:
                mm_ref_up[j]=(a1[j]*T_now[j]+a2[j]*T_out+a4[j]+a5[j]*T_now[A_z1[j]-1]+a6[j]*T_now[A_z2[j]-1]-(T_set-0.7*deadband))/(a3[j]*(T_now[j]-Tc));

        for j in range(0,17):
            if A_z2[j]==0:
                mm_ref_down[j]=(a1[j]*T_now[j]+a2[j]*T_out+a4[j]+a5[j]*T_now[A_z1[j]-1]-(T_set+0.7*deadband))/(a3[j]*(T_now[j]-Tc));
            else:
                mm_ref_down[j]=(a1[j]*T_now[j]+a2[j]*T_out+a4[j]+a5[j]*T_now[A_z1[j]-1]+a6[j]*T_now[A_z2[j]-1]-(T_set+0.7*deadband))/(a3[j]*(T_now[j]-Tc));
        
        #Get the temperature at current time
        T_sort=np.zeros((17,3))
        T_sort[:,0]=array([0]*17)
        T_sort[:,1]=array(T_now)       
        T_sort[:,2]=array(range(0,17))
        
        #Get the estimated temperature when operating at minimum airflow volume
        for j in range(0,17):
            if A_z2[T_sort[j,2]]==0:
                TT_mid_est=a1[T_sort[j,2]]*T_now[int(T_sort[j,2])]+a2[T_sort[j,2]]*T_out+a3[T_sort[j,2]]*m_min[T_sort[j,2]]*(Tc-T_now[int(T_sort[j,2])])+a4[T_sort[j,2]]+a5[T_sort[j,2]]*T_now[A_z1[T_sort[j,2]]-1];
            else:
                TT_mid_est=a1[T_sort[j,2]]*T_now[int(T_sort[j,2])]+a2[T_sort[j,2]]*T_out+a3[T_sort[j,2]]*m_min[T_sort[j,2]]*(Tc-T_now[int(T_sort[j,2])])+a4[T_sort[j,2]]+a5[T_sort[j,2]]*T_now[A_z1[T_sort[j,2]]-1]+a6[T_sort[j,2]]*T_now[A_z2[T_sort[j,2]]-1];
        #If the estimated temprature is higher than upper bound, give these zone a pirority to increase airflow volume        
            if TT_mid_est>T_set+deadband:
                T_sort[j,0]=1;
            else:
                T_sort[j,0]=0;
         
        #Sort the temperatures of all thermal zones    
        T_sorted=array(sorted(T_sort, key=lambda x: (x[0], x[1], x[2])))
        
        #This part needs to be revised as m_tar_rec[time]=m_target
        
        #Record the target total airflow volume at this time
        if time>1:
#            m_tar_rec.set_value(max(m_tar_rec.index) + 1,m_target)
             m_tar_rec.loc[len(m_tar_rec)]=m_target
             temp=pd.DataFrame()
             temp['m_tar_rec']=m_tar_rec
             m_tar_rec=temp
        else:
             m_tar_rec=pd.DataFrame()
             m_tar_rec['m_tar_rec']=ones(1)
             m_tar_rec['m_tar_rec'][0]=m_target

        #m_tar_rec=m_target
        #m_tar_rec=m_tar_rec+0
        
        Check_Cap=np.zeros((17,1))
        
        #This part needs to be revised as mm=mm_last
        #Set the initial airflow volume at this time the same as last time step
        mm=mm_last
        mm=mm+0
        
        '''
        for j in range(0,17):
            if T_now[j]<=T_set-deadband:
                mm[j]=max(min(mm_max[j],mm[j]-m_delta[j]),max(m_min[j],mm_toset[j]));
            elif T_set-deadband<T_now[j] and T_now[j]<=T_set:
                mm[j]=min(mm[j]-m_delta[j],max(m_min[j],mm_toset[j]));
            elif T_set<T_now[j] and T_now[j]<=T_set+deadband:
                mm[j]=max(mm[j]+m_delta[j],min(m_max[j],mm_toset[j]));
            elif T_now[j]>T_set+deadband:
                mm[j]=min(max(mm_min[j],mm[j]+m_delta[j]),min(m_max[j],mm_toset[j]));
        '''        
        
        for j in range(0,17):
            if mm_min[j]<m_min[j]:
                mm_min[j]=m_min[j]
                
            if mm_min[j]>m_max[j]:
                mm_min[j]=m_max[j]
                
            if mm_toset[j]<m_min[j]:
                mm_toset[j]=m_min[j]
                
            if mm_toset[j]>m_max[j]:
                mm_toset[j]=m_max[j]
                    
            if mm_max[j]<m_min[j]:
                mm_max[j]=m_min[j]
                
            if mm_max[j]>m_max[j]:
                mm_max[j]=m_max[j]
     
        #If the airflow volume is out of its range, set it to the boundary value
        for j in range(17):
            mm[j]=min(max(max(mm_min[j],m_min[j]),mm[j]),min(m_max[j],mm_max[j]))
        
        for j in range(17):
            if T_now[j]>T_set+0.7*deadband:
                mm[j]=min(max(mm[j]+3*m_delta[j],m_min[j]),m_max[j])
            elif T_now[j]<T_set-0.7*deadband:
                mm[j]=min(max(mm[j]-3*m_delta[j],m_min[j]),m_max[j])
                    
            
        
        #Get the total airflow volume
        m_total=sum(mm)
        
        #If it is already the target, no more actions are needed
        if m_total==m_target:
            Check_Cap=np.ones((17,1))
        
        #If the target is smaller than the minimum allowable range, set all to be min
        if m_target<sum(m_min):
            mm=m_min
            m_total=m_target
            m_total=m_total+0
            Check_Cap=np.ones((17,1))
            Fi=0
        #
        T_group=np.zeros((17,2))
        T_group[:,0]=array(T_now)       
        T_group[:,1]=array(range(0,17))
        
        Group_high=pd.DataFrame(columns=['0','1'])
        Group_low=pd.DataFrame(columns=['0','1'])
        ck1=0
        ck2=0
        for j in range(17):
            if T_now[j]>=T_set+0.5*deadband:
                Group_high.loc[ck1]=T_group[j]
                ck1=ck1+1
            elif T_now[j]<=T_set-0.5*deadband:
                Group_low.loc[ck2]=T_group[j]
                ck2=ck2+1
                
        
              
        #If the target is larger than the total, change the airflow volume of each thermal zone from zone with higher temperature to zone with lower temperature
        #Control Group_high first
        Fi=1
        if m_total<m_target:
            
            if len(Group_high)>0:
                fg_1=0  #len(Group_high)-1
                Check_high=zeros(len(Group_high))
                round_num=0
                while sum(Check_high)<len(Group_high):
                    if Check_high[fg_1]!=1:
                        mm[int(Group_high['1'][fg_1])]=min(min(mm[int(Group_high['1'][fg_1])]+m_delta[int(Group_high['1'][fg_1])],m_target-sum(mm)+mm[int(Group_high['1'][fg_1])]),max(min(m_max[int(Group_high['1'][fg_1])],mm_toset[int(Group_high['1'][fg_1])]),m_min[int(Group_high['1'][fg_1])]));
                    if mm[int(Group_high['1'][fg_1])]==m_max[int(Group_high['1'][fg_1])] or mm[int(Group_high['1'][fg_1])]==mm_toset[int(Group_high['1'][fg_1])]: 
                        Check_high[fg_1]=1;
                        
                    fg_1=fg_1+1
                    if fg_1==len(Group_high):
                        fg_1=0
                        round_num=round_num+1
                    if round_num==3:
                        Check_high=np.ones(len(Group_high))
                        
                    if abs(sum(mm)-m_target)<0.000001:
                        Check_high=np.ones(len(Group_high))
                        Check_Cap=np.ones((17,1))
                              
            flag=16
            while sum(Check_Cap)<17:
#                print flag
                if Check_Cap[int(T_sorted[flag,2])]!=1:
                    mm[int(T_sorted[flag,2])]=min(min(mm[int(T_sorted[flag,2])]+m_delta[int(T_sorted[flag,2])],m_target-sum(mm)+mm[int(T_sorted[flag,2])]),max(min(m_max[int(T_sorted[flag,2])],mm_max[int(T_sorted[flag,2])]),m_min[int(T_sorted[flag,2])]));
                
                if mm[int(T_sorted[flag,2])]==m_max[int(T_sorted[flag,2])] or mm[int(T_sorted[flag,2])]==mm_max[int(T_sorted[flag,2])]: 
                #or (mm[int(T_sorted[flag,2])]==m_min[int(T_sorted[flag,2])] and m_min[int(T_sorted[flag,2])]>mm_max[int(T_sorted[flag,2])]):
                    Check_Cap[int(T_sorted[flag,2])]=1;
                    
                flag=flag-1
                
                if flag==-1:
                    flag=16
                                        
                if abs(sum(mm)-m_target)<0.000001:
                    Check_Cap=np.ones((17,1))
                
            if sum(mm)!=m_target:
                Fi=0
         #If the target is larger than the total, change the airflow volume of each thermal zone from zone with lower temperature to zone with higher temperature        
        if m_total>m_target:
            if len(Group_low)>0:
                fg_1=0  #len(Group_high)-1
                Check_low=zeros(len(Group_low))
                round_num=0
                while sum(Check_low)<len(Group_low):
                    if Check_low[fg_1]!=1:
                        mm[int(Group_low['1'][fg_1])]=max(max(mm[int(Group_low['1'][fg_1])]-m_delta[int(Group_low['1'][fg_1])],m_target-sum(mm)+mm[int(Group_low['1'][fg_1])]),min(max(m_min[int(Group_low['1'][fg_1])],mm_toset[int(Group_low['1'][fg_1])]),m_max[int(Group_low['1'][fg_1])]));
                    if mm[int(Group_low['1'][fg_1])]==m_min[int(Group_low['1'][fg_1])] or mm[int(Group_low['1'][fg_1])]==mm_toset[int(Group_low['1'][fg_1])]: 
                        Check_low[fg_1]=1;
                        
                    fg_1=fg_1+1
                    if fg_1==len(Group_low):
                        fg_1=0
                        round_num=round_num+1
                    if round_num==3:
                        Check_low=np.ones(len(Group_low))
                    if abs(sum(mm)-m_target)<0.000001:
                        Check_low=np.ones(len(Group_low))
                        Check_Cap=np.ones((17,1))
            
            
            flag=0
            while sum(Check_Cap)<17:
                mm[int(T_sorted[flag,2])]=max(max(mm[int(T_sorted[flag,2])]-m_delta[int(T_sorted[flag,2])],m_target-sum(mm)+mm[int(T_sorted[flag,2])]),min(max(m_min[int(T_sorted[flag,2])],mm_min[int(T_sorted[flag,2])]),m_max[int(T_sorted[flag,2])]));
                
                if mm[int(T_sorted[flag,2])]==m_min[int(T_sorted[flag,2])] or mm[int(T_sorted[flag,2])]==mm_min[int(T_sorted[flag,2])]: 
                #or (mm[int(T_sorted[flag,2])]==m_max[int(T_sorted[flag,2])] and m_max[int(T_sorted[flag,2])]<mm_min[int(T_sorted[flag,2])]):
                    Check_Cap[int(T_sorted[flag,2])]=1;
                    
                flag=flag+1
                
                if flag==17:
                    flag=0;
                
                if abs(sum(mm)-m_target)<0.000001:
                    Check_Cap=np.ones((17,1))
                    
            if sum(mm)!=m_target:
                Fi=0
        
        #Calculate the estimate temperature with the desired airflow volume input        
        for j in range(0,17):
            if A_z2[j]==0:
                T_wp[j]=a1[j]*T_now[j]+a2[j]*T_out+a3[j]*mm[j]*(Tc-T_now[j])+a4[j]+a5[j]*T_now[A_z1[j]-1];
            else:
                T_wp[j]=a1[j]*T_now[j]+a2[j]*T_out+a3[j]*mm[j]*(Tc-T_now[j])+a4[j]+a5[j]*T_now[A_z1[j]-1]+a6[j]*T_now[A_z2[j]-1];
                
        mm_sum=sum(mm)
        
        #Calculate the estimated fan power and chiller power
        
        P_fan_wp=c1*(mm_sum)**3+c2*(mm_sum)**2+c3*(mm_sum)+c4;
        E_mix_wp=0;
        for j in range(0,17):
            E_mix_wp=E_mix_wp+mm[j]*T_now[j]
            
        #Tm_wp=((mm_sum-m_out[time])*(E_mix_wp/mm_sum)+m_out[time]*T_out)/(mm_sum);
        Tm_wp=Tm
        

        
        P_chiller_wp=cp*(mm_sum)*(Tm_wp-Tc)/(lam*COP);
        P_wp=P_fan_wp+P_chiller_wp
        T_next=T_wp;
        mm_now=mm+0;
        P_signal_now=P_signal[time]+0;
        P_now=P_wp+0;
        Tm_real=E_mix_wp/mm_sum
        
    else:
        #If signal is zero, get the estimated temperature and set all power to be 0
        mm=m_min-m_min
        Tm_wp=0
        
        for j in range(0,17):
            if A_z2[j]==0:
                T_wp[j]=a1[j]*T_now[j]+a2[j]*T_out+a3[j]*mm[j]*(Tc-T_now[j])+a4[j]+a5[j]*T_now[A_z1[j]-1];
            else:
                T_wp[j]=a1[j]*T_now[j]+a2[j]*T_out+a3[j]*mm[j]*(Tc-T_now[j])+a4[j]+a5[j]*T_now[A_z1[j]-1]+a6[j]*T_now[A_z2[j]-1];
                
        #This part needs to be revised as m_tar_rec[time]=m_target
        if time>1:
            m_tar_rec.loc[len(m_tar_rec)]=m_target
            temp=pd.DataFrame()
            temp['m_tar_rec']=m_tar_rec
            m_tar_rec=temp
        else:
            m_tar_rec=pd.DataFrame()
            m_tar_rec['m_tar_rec']=ones(1)
            m_tar_rec['m_tar_rec'][0]=m_target
        
        P_chiller_wp=0
        P_fan_wp=0
       
        P_wp=P_fan_wp+P_chiller_wp
        T_next=T_wp
        mm_now=mm+0
        P_signal_now=P_signal[time]+0
        P_now=P_wp+0
        FI=1
        Tm_est=0
        Tm_real=0
        
    mm_last=mm
    #print mm_last
    Tm_wp2=''
    Tm_wp2=Tm_wp2+str(Tm_wp)+','
    f=open('Tm.csv','a')
    f.writelines(Tm_wp2+'\n')
    f.close()
    
    print time
    m_tar_rec.to_csv('m_tar_rec.csv')
    df = pd.DataFrame(mm_last)
    df['mm_last']=mm_last
    df.to_csv("mm_last.csv")
    
    df = pd.DataFrame(T_next)
    df['T_next']=T_next
    df.to_csv("T_next.csv")
    
    return mm_now, T_next, P_signal_now, P_now, P_fan_wp, P_chiller_wp
