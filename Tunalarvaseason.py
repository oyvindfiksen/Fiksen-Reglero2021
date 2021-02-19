# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 18:43:22 2019

Code for development, feeding and growth of Atlantic Bluefin Tuna. 

@author: Øyvind Fiksen (oyvind.fiksen@uib.no)
"""


import numpy as np
import matplotlib.pyplot as plt

plotfig_env = True     #Plot fig of environmental variables
plotfig_egg = False     #Plot fig of egg survival
plotfig1larva = True   #Plot larval results 1
plotfigLightLim = False

NoLightLim=False       #Turn off light limit and larva will not be food limited
NoFoodLim = False       #Superabundant food (10000/m3)
MoreFoodLim = False     #Reduce cladoceran abundance to observed in cruises

Year03 = False      #Get temperature from year 2003 heatwave early
Year04 = False      #Get temperature from year 2004
Year11 = False
Year06 = False       #Temperature 2006, heatwave late

#Some global forcing data
TempDay = np.zeros(366)
CladoDay = np.zeros(366)
DayLength = np.zeros(366)
LightLim = np.zeros((25,366)) 
#Temperature and zooplankton measures, midmonth values.
MonthDays = np.array([30,30,31,30,31,30,31,30,31,30,31,30])
TempMidMonth = [14.8049,14.1126,14.1915, 15.4093,17.9667,21.6733,24.958,
                26.6576,25.8626, 23.2371,20.0049,16.7536]
CladoceraMidMonth= [8,16,32,77,115,162,215,140,81,57,13,12] #From de Puelles& al 2007

#Interpolation to get continuous values from day 1 to day 365. Day 0 is 0. 
for i in range(0,13):
        if i==0:
            for Jday in range(0,16):
                Jan1 = (TempMidMonth[11] + TempMidMonth[0])/2.
                Jan1C = (CladoceraMidMonth[11] + CladoceraMidMonth[0])/2.  
                DayIncr = (TempMidMonth[0] - Jan1)/15.
                DayIncrC = (CladoceraMidMonth[0] - Jan1C)/15.
                TempDay[Jday] = Jan1 + Jday*DayIncr
                CladoDay[Jday] = Jan1C + Jday*DayIncrC
        elif i < 12:
            DayIncr = (TempMidMonth[i] - TempMidMonth[i-1])/MonthDays[i]
            DayIncrC = (CladoceraMidMonth[i] - CladoceraMidMonth[i-1])/MonthDays[i]
            for j in range(0,MonthDays[i]):
                TempDay[Jday] = TempMidMonth[i-1] + j*DayIncr
                CladoDay[Jday] = CladoceraMidMonth[i-1] + j*DayIncrC
                Jday += 1
        elif i ==12:
            DayIncr = (Jan1 - TempMidMonth[11])/15.
            DayIncrC = (Jan1C - CladoceraMidMonth[11])/15.
            for j in range(0,16):
                TempDay[Jday] = TempMidMonth[i-1] + j*DayIncr
                CladoDay[Jday] = CladoceraMidMonth[i-1] + j*DayIncrC
                Jday += 1

'''Read daylength with light for every day from file'''
with open('HoursofLight.txt', 'r') as reader:
        Jday=0
        for line in reader:
            DayLength[Jday] = line
            Jday += 1
''' Read temperature data from file'''
with open('TempDataPRB.txt', 'r') as reader:
       Jday=0
       for line in reader: 
           TempDay[Jday] = line
           Jday += 1
        
#Change temperature? Food? Daylength?
TempDay[:] += 0.

if(Year03==True): TempDay[:]=np.fromfile("temp2003.npy")
if(Year04==True): TempDay[:]=np.fromfile("temp2004.npy")
if(Year11==True): TempDay[:]=np.fromfile("temp2011.npy")
if(Year06==True): TempDay[:]=np.fromfile("temp2006.npy")



'''Light limitation of the feeding process is between 1 and 0,
1 if light is above ** and 0 else. In summer days are longer. 
n the morning and the afternoon a fraction of the hour can be dark, 
and the rest illuminated.'''
for i in range(0,366):        
         dawn = int(12 - DayLength[i]//2)  #int(12 - (int(DayLength[i])/2))
         dusk = int(12 + DayLength[i]//2)  #int(12 + (int(DayLength[i])/2))
         hourfraction = DayLength[i] - 2*(DayLength[i]//2)
         LightLim[dawn:dusk,i] = 1. #no light lim in day
         LightLim[dawn-1,i] = hourfraction/2.
         LightLim[dusk,i] = hourfraction/2.
         
         '''a	220,7472	b	41,3084	c	1,2717 x0	190,8144
         Fitted modified Gaussian curve to monthly data of Cladocerans (ref)'''
         CladoDay[i] = 220.74*np.exp(-0.5*np.abs((i-190.81)/41.31)**1.2717)

if(NoLightLim==True): LightLim[:,:]   = 1.  #Can feed all night and day 
if(NoFoodLim==True): CladoDay[:]   = 10000. #Superabundant food 
if(MoreFoodLim==True): CladoDay[:] *= 0.1   #0.1 here shifts total cladocera abndance to field sampling values

def egg_dev_time(Temperature :float):
    """Returns the number of hours it takes for a bluefin tuna 
       egg to hatch at a given temperture. From Reglero & al 2018
       Minimum viable Temperature is 19 degrees C"""
    return (8787.5*(Temperature**-1.701))

def egg_hatching_prob(Temperature :float):
    """Returns the hatching probability (not %) for a bluefin tuna 
       egg at a given temperature. From Reglero & al 2018. 
       Minimum viable Temperature is 19 degrees C"""
    return max(0.,(-1.27*Temperature**2. + 63.78*Temperature -727.98)*0.01)

def egg_surv(Temp):
    """Returns egg survival chance from mortality rate and development time
    """
    dh = 1./24.
    #egg_m = 0.5*dh #Egg mortality rate per hour (Daies & al 1991 = 0.5)
    egg_m = 0.5*dh*0.00022*(0.001*0.018)**-0.85 # From McGurk 1986
    egg_pred_surv = np.exp(-egg_m*egg_dev_time(Temp))
      
    return egg_pred_surv

def SGR_T(Temp):
    """Returns the specific per hour growth rate of ABFT larvae as a function of 
    temperature. From Reglero & al 2018, their fig 1. 
    NB - in the lab exp they only feed during 15 hours of light, and eat nothing during night.
    In the lab, their guts run with empty within a couple of hours after darkness (Reglero, pers com), 
    and growth in mass must then logically cease. This means we have to distrbute 
    the average daily growth potential over ca 15 hours, when they have food in their gut, and not 24.
    If the larvae have to long days in the lab, they do not thrive (Reglero, pers com) """
    return max(0.,0.0418*Temp - 0.8355)/15.

def larva_surv(wgt,dh):
    """ Returns larval chance of surviving dh days as a function of 
    body mass in mg DW (McGurk 1986). Converts to gram! Note error in Reglero&al2018, Fig. 1d legend text"""
    m = 0.5*0.00022*(0.001*wgt)**-0.85
    #m = 0.1             #Or fixed mortality
    m = np.exp(-m*dh)
    return m

def larva_feeding(day,hour,wgt):
    """The feeding process of larva - clearance from visual acuity. Return mass of prey encontered in mg dw per hour.
    Include light, prey density, larval size.. Handling time can limit feeding on nauplii, but is not a 
    factor for the larger cladocera (rare encounters, large prey)"""
    #hour = 12;     day = 180;     wgt = 0.018
    
    length_m = np.log(wgt/0.0008)/905.2
    handling = 2./3600          #s ->h     handling time, hours per prey
    DwClad = 3.15*1E-03         #mgDW DW_prey_micg->mg  (ref?)
    Clad_length = 0.8*1E-3       #Cladoceran length in mm
    swim_vel=3.*length_m*3600   #BL/s ->m/h swimming velocity (m/s) scale with BL
    msa = 4.699*((length_m*1000.)**-1.129)  #From Hilder PE, Battaglene SC, Hart NS, Collin SP, Cobcroft JM (2019) Retinal adaptations of southern bluefin tuna larvae: Implications for culture. Aquaculture 507:222-232
    behav_anatomical_ratio = 0.5 #The ratio between anatomical and behavioural measures of visual detection
    vis_r=0.5*Clad_length/np.tan(msa*0.5*2*np.pi/360) #See eg Job & Bellwood 1996 JFB 48:952-963. vis_r in m
    vis_r=behav_anatomical_ratio*vis_r         #m   How far away the larva detect prey
    Light = LightLim[hour,day]  #Day, night, dusk or dawn
    Clado_dens = CladoDay[day] * 0.1   # #/m3 Density of cladceran - reduced in open ocean by some factor
    Pr_capt_clad = min(1., max(0.,(wgt - 0.018)/(0.77-0.018))) #Capture probability of cladocerans
    
    naup_dens = 400. # 10.*Clado_dens  #nauplii/m3 - smaller larvae only eat nauplii Density is from field surveys
    naup_len =0.3*1E-3          #typical size of nauplii, Oithona 
    DwNaup = 0.5*1E-03  #Dw in mg from Hay &al 1988, 0.3 mm long naplia
    vis_r_naup = 0.5*naup_len/np.tan(msa*0.5*2*np.pi/360)
    vis_r_naup = behav_anatomical_ratio*vis_r_naup         #m   How far away the larva detect prey
    clearance_vol_naup= 0.5*np.pi*(vis_r_naup**2.)*swim_vel*Light #0.5 - larva only look up
    naup_enc = DwNaup*clearance_vol_naup*naup_dens/(1.+handling*clearance_vol_naup*naup_dens)
    
    clearance_vol_Clad = 0.5*np.pi*(vis_r**2.)*swim_vel*Light #m3/h, daylight is on or off
    food_enc = Pr_capt_clad*DwClad*clearance_vol_Clad*Clado_dens/(1.+handling*clearance_vol_Clad*Clado_dens) #Holling disk
    food_enc += naup_enc # Larvae eat both cladcerans and nauplii, but capture of Clado is size-dependent
    return food_enc

def metab_rate(wgt,Temp):
     """Find metabolic rate in mgDW per hour from body mass. 
     Weight in is in mg DW """
     q10 = 2.    #Q10 - if q10 is higher then the stomach-based model will start to deviate
                 #from the temperature dependency in oberved potential growth rate (Reglero & al. 2018)
     metabol_rate = 0.404*(wgt**0.994) # in umol O2/h at 26 degr (Reglero & al in revision)
     metabol_rate *= 32.*0.88 # -> ug O2 * RQprotein = 0.88  
     metabol_rate *= 12./34. #->ug C
     metabol_rate *= 100./(1000.*45.) #-> ug C -> mg DW (mgC=45% of DW)
     metabol_rate *= q10**((Temp-26.)/10.)  #Measured at 26 degrees only. Q10 fitted from growth 
     return metabol_rate
 
def yolk_duration(Temp):
    """Return the number of days it takes to develop through the yolk sac stage
    assuming q10 and 2.5 days at 25 degrees"""
    q10=2.
    T_measured = 25.    # The temperature where development time has been measured to 2.5 days
    Yolk_stage_days= 2.5*q10**((T_measured-Temp)/10.)
    return Yolk_stage_days
  

'''Find fitness (=survival from egg to flexion stage) of all spawning dates of the year. 
     We start one tuna egg each day of the year and follow it forward in time 
     until it reach the flexion stage. The egg has to hatch successfully - which require 
     temperatures above 18 degrees. Then the yolk sac period, and the first feeding, 
     from 0.018 mgDW to 0.77 mgDW. We also include a size-dependent mortality, 
     which gives a penalty for longer stage duration. Feeding larvae include stomach dynamics.
     '''
Year = np.arange(0,366,1,dtype = 'int')
diel_res = np.arange(1,25,1, dtype ='int')
dh = 1./24.                         #Diel resolution = 1hour
age =  np.array([0.]*366)           #age since spawned
larval_age = np.array([0.]*366)     #age since hatch
stomach   = np.array([0.]*366)      #stomach content from birthday
gut = np.zeros((366,366,24))        #follow gut of birthday 24/7
gut2 = np.zeros((366,366*24))       #Follow gut of larva continuously
bodym = np.zeros((366,366*24))
L_ageH = np.zeros((366,366*24))
larval_avgT = np.array([0.]*366)    #average temperature during larval stage
wgt_mgDW= np.array([0.018]*366)     #Larvae hatch with body mass 0.018 mg DW
surv :float = np.ones(366)          #Total survival chance
Psurv = np.ones((366,366*24))
Ps_yolk = np.ones(366)
eggsurv = []                        #Survival through egg stag
pr_hatch, timetohatch = [], []      #Hatching success, hours to hatch
#Yolk_duration = 2.5                 #Days after hatching before feeding starts
AE =    0.7                            #assimilation efficiency (ref?)
EndSize = 0.77                      #ug DW 
MaxGut = 0.12                         #Gut size, in fractions of body mass


#Release one egg per day - at midnight. Then loop through until 
for day in Year:
     eggsurv.append(egg_surv(TempDay[day])) #Egg survival from predation
     pr_hatch.append(egg_hatching_prob(TempDay[day])) #Hatching survival 
     timetohatch.append(egg_dev_time(TempDay[day]))  #Adding egg dev time to age - hours
     a_count = 0

     surv[day] = (egg_hatching_prob(TempDay[day])*egg_surv(TempDay[day]))
     
     age[day]=egg_dev_time(TempDay[day])*dh    #Hatching time in days
     Yolk_days = yolk_duration(TempDay[day])   #Duration of yolk stage

     age[day] += Yolk_days          #Time between hathing and first feeding
     surv[day] *= larva_surv(wgt_mgDW[day],Yolk_days) #Survival durng the non feeding period
     Ps_yolk[day] = surv[day]
     #For each potential birthday (day), move to hatching time and follow survival until given size
     hatch_day :int = day +  int(Yolk_days + egg_dev_time(TempDay[day])//24)    #Given timesteps of 24
     hatch_hour :int = int(Yolk_days + egg_dev_time(TempDay[day]) % 24) #Assume all hatch at midnight - NB the adaptive value of spawning at night time..
          
     #Loop over the remaining fraction of the day after hatching
     for hour in range(hatch_hour,24):
         age[day] += dh
         larval_age[day] += dh
         larval_avgT[day] += TempDay[min(hatch_day,365)]
         growthT = wgt_mgDW[day]*SGR_T(TempDay[min(hatch_day,365)])
         GutLim = MaxGut*wgt_mgDW[day]
         
         #The gut is a container which determines if the larva can grow at max rate, limited by temperature or food
         Food_enc = larva_feeding(day,hour,wgt_mgDW[day]) #food encountered, in mgDW
         if(growthT == 0.): Food_enc = 0.       #No feeding if to cold to grow
         
         #the digestion is driven by growth, metabolism and assimilation efficiency
         resp = metab_rate(wgt_mgDW[day],TempDay[min(hatch_day,365)])
         digestion = (resp + growthT)/AE   
         
         # if the growth is temperature limited; enough food in the stomach 
         if(stomach[day] >= digestion):
             wgt_mgDW[day] += growthT
             stomach[day]  += (Food_enc - digestion)     #Add prey remove digestion
             stomach[day]   = min(stomach[day],GutLim) #Constrain gutcontent
             gut[day,min(hatch_day,365),hour]  = stomach[day]/wgt_mgDW[day]
         else:
             wgt_mgDW[day] += (stomach[day]*AE - resp) #What is left in stomach is used to grow or cover respiration 
             stomach[day]   = Food_enc                 #Empty stomach, add feeding
             stomach[day]   = min(max(0.,stomach[day]),GutLim)   #and it may fill up the stomach by feeding
             gut[day,min(hatch_day,365),hour]  = stomach[day]/wgt_mgDW[day]
         
         wgt_mgDW[day] = max(wgt_mgDW[day],0.018)   #Do not shrink below hatch size
         surv[day] *= larva_surv(wgt_mgDW[day],dh)  #survival chance from size
         if(hatch_day>360): surv[day]= 0.           #If the larva do not reach flexion by the end (avoid terminal effects)

         a_count += 1                           #Variables to look at results and plot
         gut2[day,a_count] = stomach[day]
         L_ageH[day,a_count] = age[day]
         bodym[day,a_count] = wgt_mgDW[day]
         Psurv[day,a_count] = surv[day]


     #Then loop over the rest of the year, or until larva reach 0.77 mg DW
     for j in range(hatch_day+1,365):
         #Loop over 24 hours diel cycle
         for h_time in range(24): 
             if wgt_mgDW[day]>=EndSize: continue    #If flexion, done
             if(j>360): surv[day]= 0.    #If the larva do not reach flexion by the end

             GutLim = MaxGut*wgt_mgDW[day]

             age[day] += dh
             larval_age[day] += dh
             larval_avgT[day] += TempDay[j]     
             growthT = wgt_mgDW[day]*SGR_T(TempDay[j])
              
             #The gut is a container which determines if the larva can grow at max rate, limited by temperature or food
             Food_enc = larva_feeding(j,h_time,wgt_mgDW[day]) #food encountered, in mgDW
             if(growthT == 0.): Food_enc = 0.       #No feeding if to cold to grow

             #the digestion is driven by growth, metabolism and assimilation efficiency
             resp = metab_rate(wgt_mgDW[day],TempDay[j])
             digestion = (resp + growthT)/AE   
         
             # if the growth is temperature limited; enough food in the stomach 
             if(stomach[day] >= digestion):
                 wgt_mgDW[day]  += growthT
                 stomach[day]   += (Food_enc - digestion)     #Add prey remove digestion
                 stomach[day]    = min(stomach[day],GutLim) #Constrain gutcontent  #digestion is driven by growth
                 gut[day,j,h_time] = stomach[day]/wgt_mgDW[day]
             else:
                 wgt_mgDW[day]  += (stomach[day]*AE - resp)   #What is left in stomach is used to grow or cover respiration 
                 stomach[day]    = Food_enc                 #Empty stomach, add feeding
                 stomach[day]    = min(max(0.,stomach[day]),GutLim)   #and it may fill up the stomach by feeding
                 gut[day,j,h_time] = stomach[day]/wgt_mgDW[day]
                  
             wgt_mgDW[day] = max(wgt_mgDW[day],0.018)   #Do not shrink below hatch size      
             surv[day] *= larva_surv(wgt_mgDW[day],dh)  #Survive each hour
             
             a_count += 1                           #the rest here is for plotting
             gut2[day,a_count] = stomach[day]
             L_ageH[day,a_count] = age[day]
             bodym[day,a_count] = wgt_mgDW[day]
             Psurv[day,a_count] = surv[day]

#******************************************************************************
# The rest is only plotting results
print('best day', 'survival')
print(np.argmax(surv), max(surv))

if(plotfigLightLim==True):
    fig = plt.figure(figsize=(6,5))
    left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
    ax = fig.add_axes([left, bottom, width, height]) 
    
    cp = ax.contour(LightLim)
    ax.clabel(cp, inline=True, 
              fontsize=10)
    contour_filled = plt.contourf(LightLim, Levels=0, colour = 'k')
    plt.colorbar(contour_filled)
    ax.set_title('Daylight - Light limitation on feeding')
    ax.set_xlabel('Day of year')
    ax.set_ylabel('Hour of day')
    fig.savefig('Daylight.svg', dpi = 600)

    plt.show()


fig_gut, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,
         gridspec_kw={'hspace': 0},figsize=(6,6),constrained_layout=True)
fig_gut.suptitle('Gut dynamics', fontsize=16)
bm2 = bodym.copy(); a2=L_ageH.copy(); g2= gut2.copy()
bm2[bodym == 0.] = np.nan; a2[L_ageH == 0.] = np.nan; g2[gut2 == 0.] = np.nan
ax2.plot(a2[150,1:1400],gut2[150,1:1400],label='150')
ax2.plot(a2[180,1:1400],gut2[180,1:1400],label='180')
ax2.plot(a2[200,1:1400],gut2[200,1:1400],label='200')
ax2.set_ylabel('Gut mass')
ax2.set_xlabel('Age (days)')
ax2.set_title('Stomach content')
ax2.legend()

ax1.set_title('Growth feeding stage')
ax1.set_xlabel('Age (days)')
ax1.set_ylabel('Body mass')
ax1.plot(a2[150,1:1400],bm2[150,1:1400],label='150')
#ax1.plot(a2[160,1:1400],bm2[160,1:1400],label='160')
ax1.plot(a2[180,1:1400],bm2[180,1:1400],label='180')
ax1.plot(a2[200,1:1400],bm2[200,1:1400],label='200')
#ax1.plot(a2[220,1:1200],bm2[220,1:1200],label='220')
#ax1.plot(a2[230,1:1200],bm2[230,1:1200],label='230')
ax1.legend()

ax3.set_ylabel('Age (days')
ax3.set_xlabel('Hour of day'); ax4.set_xlabel('Hour of day')
ax3.set_title('Gut, born day 150');ax4.set_title('Gut, born day 180')
ax3.contour(gut[150,150:210,:])
ax4.contour(gut[180,180:240,:])

fig_gut.savefig('fig_gut.svg', dpi = 300)
fig_gut.savefig('fig_gut.png', dpi = 300)

plt.show(fig_gut)


#Figure of survival
fig_surv, (ax1,ax2) = plt.subplots(2,1,constrained_layout=True)
fig_surv.suptitle('Survival', fontsize=16)
ps = Psurv.copy(); ps[surv == 1.] = np.nan
ax1.set_title('Survival from first feeding')
ax1.set_xlabel('Age (days)')
ax1.set_ylabel('Survival probability')
ax1.set_yscale('log')
ax1.set_ylim(0.1,1.E-9)
ax1.set_xlim(100,350)
ax1.plot(150 + a2[150,1:1200],ps[150,1:1200],label='150')
ax1.plot(160 + a2[160,1:1200],ps[160,1:1200],label='160')
ax1.plot(180 + a2[180,1:1200],ps[180,1:1200],label='180')
ax1.plot(200 + a2[200,1:1200],ps[200,1:1200],label='200')
ax1.plot(220 + a2[220,1:1200],ps[220,1:1200],label='220')
ax1.plot(240 + a2[240,1:1200],ps[240,1:1200],label='240')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
          borderaxespad=0.)
ax2.set_yscale('log')
ax2.set_ylim(0.1,1.E-9)
ax2.set_xlim(100,350)
ax2.set_title('Survival through stages')
ax2.set_xlabel('Day of the year')
ax2.set_ylabel('Survival probability')
ax2.plot(Year,Ps_yolk,label='Ps yolk')
ax2.plot(Year[132:340],surv[132:340]/Ps_yolk[132:340],label='Ps feed')
ax2.plot(Year,surv,label='Ps total')

ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
          borderaxespad=0.)
fig_surv.savefig('fig_surv.svg', dpi = 300)

plt.show()



if(plotfig1larva==True):
    sgr = np.zeros(366)
    SMR = np.zeros(366)
    tot_egg_surv= [0.]*366 #np.zeros(366) 
    food_intake = np.zeros([24,366])
    food_intake_day = np.zeros(366)
    food_intake_day24 = np.zeros(366)
    Food_intake_day_wgt010 = np.zeros(366)
    Food_intake_day_wgt050 = np.zeros(366)
    surv_Yolk = np.zeros(366)
    Yolkdays = np.zeros(366)
    wgt= 0.77 #some of the plots require a given body mass

    for i in range(0,366):
        sgr[i] = SGR_T(TempDay[i])/dh
        tot_egg_surv[i] = eggsurv[i]*pr_hatch[i]
        Yolkdays[i] = yolk_duration(TempDay[i])
        surv_Yolk[i] = larva_surv(0.018,yolk_duration(TempDay[i]))
        SMR[i] = metab_rate(wgt,TempDay[i])/(wgt*dh)
        
        for hour in range(24):
            feeding = larva_feeding(i,hour,wgt)
            food_intake[hour,i]=min(wgt*MaxGut,feeding)
            food_intake_day[i] += feeding
            food_intake_day24[i] += larva_feeding(i,12,wgt)
            Food_intake_day_wgt010[i] += larva_feeding(i,hour,0.1)
            Food_intake_day_wgt050[i] += larva_feeding(i,hour,0.5)
            
    fod=length=respW = np.zeros(100)
    for j in range(10,100):
        fod[j] = larva_feeding(180,12,j*0.01)/(wgt*dh)  #Specific daily prey enc
        length[j] = np.log(j*0.01/0.0008)/905.2
        respW[j] = metab_rate(j*0.01,26.)/(wgt*dh)      #Specific daily resp at 26 C
        
    eggdaysT= np.zeros(12);yolkdaysT= np.zeros(12);T_range = np.zeros(12)
    respT= np.zeros(12); gr_rateT= np.zeros(12);hatchprob = np.zeros(12)
    #f=open('Temperaturefunctions.txt', 'w')
    #f.write('Temperature, metabolism, eggduration, yolkduration,growthrate,hatchingsuccess \n')
    for j in range(12):
        T_range[j]=j+18
        respT[j] = metab_rate(wgt,(j+18.))/(wgt*dh)
        eggdaysT[j] = egg_dev_time(j+18.)
        yolkdaysT[j] = yolk_duration(j+18.)/dh
        gr_rateT[j] = SGR_T(j+18.)/dh
        hatchprob[j] = egg_hatching_prob(j+18.)
    
    
#    Tdata = [T_range, respT,eggdaysT,yolkdaysT,gr_rateT,hatchprob] 
#    with open("Temperaturefunctions.txt", "w") as file:
#        file.write('Temperature, metabolism, eggduration, yolkduration,growthrate,hatchingsuccess \n')
#        for x in zip(*Tdata):
#            file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(*x))
        
#    file.close()
    #food_intake_day.tofile("Food_intake_day_12hoursoflightallseason.npy")
    #surv.tofile("Surv with no food lim.npy")
    #surv_onlyT = np.fromfile("Surv with no food lim.npy")
    #foodintake12hoursoflight = np.fromfile("Food_intake_day_12hoursoflightallseason.npy")
     
    fig_larva, ax = plt.subplots(constrained_layout=True)
    fig_larva.suptitle('Fitness, feeding, SGR and egg surv', fontsize=16)

    ax1 = ax2 = ax.twinx()
    ax.set_ylabel("Diel prey encounter/body mass, wgt= %s" %  round((wgt),3))
    ax.set_xlabel("Day of year (spawning)")
    ax1.set_ylabel('Relative fitness of birthday')
    #ax.plot(Year,timetohatch,'r-',label='Time to hatch (hours)')
    ax2.plot(Year,sgr,'r-', label ='SGR(T), wgt= %s' %  round((wgt),3))
    ax.plot(Year,food_intake_day/wgt,'g-', label ='Specific prey enc')
    #ax.plot(Year,foodintake12hoursoflight,'r-', label ='Prey enc12hourslight')
    ax1.plot(Year,surv/max(surv),'b-', label ='Fitness of birthday')
    ax1.plot(Year,tot_egg_surv/max(tot_egg_surv),'y-', label ='Egg survival')
    handles, labels = ax.get_legend_handles_labels()
    fig_larva.legend(bbox_to_anchor=(0.15, 1.02, 0.7, 0.102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
    fig_larva.savefig('fig_larva.svg', dpi = 300)

    plt.show()     
    
    if(Year03==True): surv.tofile("survT2003.npy")
    if(Year04==True): surv.tofile("survT2004.npy")
    if(Year06==True): surv.tofile("survT2006.npy")


    #surv_onlyT = np.fromfile("Surv with no food lim.npy")


## Fig of feeding processes
    fig_feeding, (ax1,ax3) = plt.subplots(2,sharex='col',sharey=False,
                 gridspec_kw={'hspace': 0},figsize=(8,10),constrained_layout=False)
    #ax3=ax2.twinx()
    #ax2.set_ylabel("SGR and SMR")
    ax3.set_xlabel("Day of year")
    ax1.set_ylabel('Prey enc /W /day')
    ax3.set_ylabel('Prey encounter' '\n' 'and energy surplus')

    #ax.plot(Year,timetohatch,'r-',label='Time to hatch (hours)')
   # ax2.plot(Year,sgr, label ='SGR')
    #ax2.plot(Year,SMR, label ='SMR')
    
    ax3.set_ylim(-0.1, ymax=1.)
    ax1.set_ylim(0,1)
    ax3.plot(Year,(sgr+SMR)/AE,label ='(SGR+SMR)/AE')
    ax3.plot(Year,food_intake_day/0.77, label ='I/W')
    ax3.plot(Year,((food_intake_day/0.77)-((sgr+SMR)/AE)), label ='I/W-(SGR+SMR)/AE')
    ax3.plot(Year,((food_intake_day/0.77)-(SMR/AE)), label ='I/W-(SMR/AE)')
    import matplotlib.patches as patches
    ax3.add_patch(patches.Rectangle((160, 0), # (x,y) 
            50, # width
            1, # height
            alpha=0.3, facecolor="grey", edgecolor="black", linewidth=0, linestyle='solid'
            ))

    ax1.plot(Year,food_intake_day/wgt, label ='W=0.77',linewidth=3)
    ax1.plot(Year,food_intake_day24/wgt, label ='W=0.77 24h light')
    ax1.plot(Year,Food_intake_day_wgt050/0.5, label ='W=0.5')
    ax1.plot(Year,Food_intake_day_wgt010/0.1, label ='W=0.1')

    ax1.legend(loc="best"); ax2.legend(loc="best"); ax3.legend(loc="best") 
    fig_feeding.savefig('fig_feeding.png', dpi = 300)
    fig_feeding.savefig('fig_feeding.svg', dpi = 300)

    plt.show()       

## Figure of some temperature dependencies    
    fig_physiologyT, ax = plt.subplots(constrained_layout=True)
    ax1 = ax.twinx()
    ax.set_ylabel(r"SGR, Meatabolic rate (/day)" "\n" r"Hatching probability")
    ax.set_xlabel("Temperature)")
    ax1.set_ylabel('Stage duratation (hours)')
    ax.plot(range(18,30),gr_rateT,'r-', label ='SGR(T)')
    ax.plot(range(18,30),respT,'b-', label ='SMR(T)')
    ax.plot(range(18,30),hatchprob,'g-', label ='Hatch success')
    ax1.plot(range(18,30),eggdaysT,'y-', label ='Egg stage duration')
    ax1.plot(range(18,30),yolkdaysT,'k', label ='Yolk stage duration')
    handles, labels = ax.get_legend_handles_labels()
    fig_physiologyT.legend(bbox_to_anchor=(0.15, 1.02, 0.7, 0.102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
    fig_physiologyT.savefig('fig_physiologyT.svg')
    plt.show()         

if (plotfig_env == True):
    #Figure of data with three axes
    from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes

    fig_env = plt.figure()

    host = HostAxes(fig_env, [0.15, 0.1, 0.8, 0.8])
    par1 = ParasiteAxes(host, sharex=host)
    par2 = ParasiteAxes(host, sharex=host)
    host.parasites.append(par1)
    host.parasites.append(par2)

    host.set_ylabel("Temperature")
    host.set_xlabel("Day of year")

    host.axis["right"].set_visible(False)
    par1.axis["right"].set_visible(True)
    par1.set_ylabel("Cladocera")

    par1.axis["right"].major_ticklabels.set_visible(True)
    par1.axis["right"].label.set_visible(True)

    par2.set_ylabel("Daylength")
    offset = (50, 0)
    new_axisline = par2.get_grid_helper().new_fixed_axis
    par2.axis["right2"] = new_axisline(loc="right", axes=par2, offset=offset)

    fig_env.add_axes(host)

    host.set_xlabel("Day of the year")
    host.set_ylabel("Temperature")
    par1.set_ylabel("Cladocera")
    
    temp2003=np.fromfile('temp2003.npy'); temp2004=np.fromfile('temp2004.npy')

    p1, = host.plot(range(1,365),TempDay[1:365],'b-', label='Temp')
    p1, = host.plot(range(1,365),temp2003[1:365],'b-.', label='Temp 2003')
    p1, = host.plot(range(1,365),temp2004[1:365],'b--', label='Temp 2004')
    p4, = host.plot(range(15,350,30), TempMidMonth[0:12], 'ro', label='Temp data')
    p5, = par1.plot(range(15,350,30),CladoceraMidMonth[0:12],'go', label='Clad data')
    p2, = par1.plot(range(1,365),CladoDay[1:365],'g-', label='Cladocera')
    p3, = par2.plot(range(1,365),DayLength[1:365],'y-', label='Daylength')

    host.legend()

    host.axis["left"].label.set_color(p1.get_color())
    par1.axis["right"].label.set_color(p2.get_color())
    par2.axis["right2"].label.set_color(p3.get_color())
    plt.title('Data on temperature and Cladocerans')
    fig_env.savefig('Env.svg')
    plt.show()

if(plotfig_egg == True):
    #figure of egg hatching and survival
    fig_egg, ax = plt.subplots(constrained_layout=True)
    ax1 = ax.twinx()
    ax.set_ylabel("Stage duration (hours)")
    ax.set_xlabel("Day of year (spawning)")
    ax.set_ylim(0, ymax=130)
    ax1.set_ylabel("Survival probability")
    ax.plot(Year,timetohatch,'r-',label='Age at hatching (hours)')
    ax.plot(Year,Yolkdays/dh,'k-.',label='Yolk stage (hours)')
    ax1.plot(Year,eggsurv,'b-', label ='Egg predation survival')
    #ax1.plot(Year,eggsurvMcG,'b-.', label ='Egg predation survival McG')
    ax1.plot(Year,pr_hatch,'g-', label ='Hatching probability')
    ax1.plot(Year,surv_Yolk,'k-', label ='Survival yolk stage')
    #ax1.plot(Year,surv/max(surv),'y-', label ='Relative total survival')
    ax1.plot(Year,tot_egg_surv,'y-', label ='Total egg survival')
    #ax1.plot(Year,tot_egg_survMcG,'y-.', label ='Total egg survival McG')
    handles, labels = ax.get_legend_handles_labels()
    fig_egg.legend(bbox_to_anchor=(0.11, 1.02, 0.8, 0.102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
    plt.show()
    
#plot fitness as chance of survival -Fig 1
#**************
#from matplotlib.ticker import ScalarFormatter

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

fig_fitness, (hostA,hostB,hostC) = plt.subplots(3,
             figsize=(9,8), sharex=True, constrained_layout=False)
fig_fitness.subplots_adjust(hspace=0.1,right=0.75)

#Panel A
ax1A = hostA.twinx()
ax2A = hostA.twinx()
ax3A = hostA.twinx()
ax1A.set_ylabel(r"Body mass mg DW")
ax2A.spines["right"].set_position(("axes", 1.1))
ax2A.spines["right"].set_visible(True)
ax3A.spines["right"].set_position(("axes", 1.2))
ax3A.spines["right"].set_visible(True)
ax3A.set_yscale('log');ax3A.set_ylabel("Chance of survival with age", color='grey'),ax3A.set_ylim(1e-10,1)

hostA.set_xlim(150,260); #hostA.set_ylim(100,250)
ax2A.set_ylim(10,15); ax2A.set_ylabel("Daylenght (hours)",color='orange');hostA.set_ylabel(r"Cladocera (m" r"$^{-3}$)", color='green')
hostA.plot(range(1,365),CladoDay[1:365],color='green', label='Cladocera')
hostA.plot(range(15,350,30),CladoceraMidMonth[0:12],'go', label='Clad data')

ax2A.plot(range(1,365),TempDay[1:365],'r-', label='Temp Avg')
ax2A.plot(range(1,365),DayLength[1:365],'y--', label='Daylength')

ax3A.plot(160 + a2[160,1:1200],ps[160,1:1200],color='grey',dashes=[2, 2],label='Born 160')
ax3A.plot(180 + a2[180,1:1200],ps[180,1:1200],color='grey',dashes=[2, 2],label='Born 180')
ax3A.plot(200 + a2[200,1:1200],ps[200,1:1200],color='grey',dashes=[2, 2],label='Born 200')
ax3A.plot(220 + a2[220,1:1200],ps[220,1:1200],color='grey',dashes=[2, 2],label='Born 220')
#ax3A.plot(225 + a2[225,1:1200],ps[225,1:1200],color='grey',dashes=[2, 2],label='Born 225')

ax1A.plot(160 + a2[160,1:1200],bm2[160,1:1200],'k-',label='160')
ax1A.plot(180 + a2[180,1:1200],bm2[180,1:1200],'k-',label='180')
ax1A.plot(200 + a2[200,1:1200],bm2[200,1:1200],'k-',label='200')
ax1A.plot(220 + a2[220,1:1200],bm2[220,1:1200],'k-',label='220')
#ax1A.plot(225 + a2[225,1:1200],bm2[240,1:1200],'k-',label='225')

#Panel B***************
par1B = hostB.twinx()
par1B.set_xlim(150,260); par1B.set_ylim(0.005,1000); par1B.set_yscale('log')
#surv.tofile("surv_300_010.npy") #save results for printing
nofLim=np.fromfile('nofLim.npy'); surv_lowm01=np.fromfile('surv_lowm01.npy')
n500_C015=np.fromfile('surv_500_015.npy'); n400_C010=np.fromfile('surv_400_010.npy');n300_C010=np.fromfile('surv_300_010.npy')

# get the field data
Larvaem3 = [0.0139,0,0,0,0,0,0,5.26E-03,0,0,0.0689,9.38E-03,0.0144,0.187,0.3411,0.9099,0.2281,16.399,0.2922,17.9296,0.0445,1.6697,0.3689,0.4967,1.3426,1.4146,0.2417,1.2391,0.0362,0.149,26.7838,0.2874,2.4841,0.0408,5.07E-03,0.4819,0.8815,2.75E-03,0.1331,5.54E-03,0.0199,5.98E-03,0.1512,0.0364,0.2952,0,0.043,0,0,0.1959,0,0,0.0132,7.09E-03,6.70E-03,5.10E-03,2.91E-03,4.94E-03,0,0,0,0,0]
SampleDates = [158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,212,213,214,215,216,217,219,221,222,223]

p2B, = hostB.plot(Year,nofLim/max(nofLim),color='blue', dashes=[2,2], label ='No food lim')
p2B, = hostB.plot(Year,n500_C015/max(n500_C015),'b-',label ='n500_C015')
p2B, = hostB.plot(Year,n400_C010/max(n400_C010),color='blue', dashes=[3,3], label='n400_C010')
p2B, = hostB.plot(Year,n300_C010/max(n300_C010),color='blue', dashes=[1,1], label='n300_C010')

p4B = par1B.bar(SampleDates,Larvaem3, color='lightgrey', label='Field data')
hostB.legend(loc='lower right',prop={'size': 8})
#hostB.legend(loc='upper left'); 
par1B.set_ylabel(r"Larval density field (m" r"$^3$)", color='grey')
hostB.set_ylim(0,1.05);hostB.set_ylabel("Relative egg fitness", color='blue')

#Panel C***********************
ax1C = hostC.twinx()
hostC.set_ylabel("Egg fitness              ")
hostC.set_ylim(0,2E-5)
hostC.set_yticks([0,4E-6,8E-6,12E-6])
#hostC.ticklabel_format(style='sci')
hostC.set_xlim(150,260);hostC.set_xlabel("Day of year")
ax1C.set_yticks(ticks=[22,24,26,28]);ax1C.set_ylim(16,29);ax1C.set_ylabel("               Temperature")

temp2003=np.fromfile('temp2003.npy'); temp2004=np.fromfile('temp2004.npy'); temp2006=np.fromfile('temp2006.npy')
temp2011=np.fromfile('temp2011.npy');survT2011 = np.fromfile('survT2011.npy')
survT2004=np.fromfile('survT2004.npy');survT2003 = np.fromfile('survT2003.npy');survT2006 = np.fromfile('survT2006.npy')
ax1C.plot(range(1,365),TempDay[1:365],'green', label='Average')
ax1C.plot(range(1,365),temp2003[1:365],'red', label='2003')
ax1C.plot(range(1,365),temp2004[1:365], 'blue',label='2004')
ax1C.plot(range(1,365),temp2006[1:365], 'orange',label='2006')
#hostC.plot(range(1,365),temp2011[1:365],'r', dashes=[1, 2],label='Temp 2011')
ax1C.legend(loc = 'best')
hostC.plot(Year,n500_C015,'g--', label ='Average')
hostC.plot(Year,survT2003,'r--', label ='2003')
hostC.plot(Year,survT2004, 'b--', label ='2004')
hostC.plot(Year,survT2006, color ='orange', label ='2006')
#ax1C.plot(Year,survT2011,'k', dashes=[1, 2],label ='Surv 2011')
#ax1C.legend(loc = 'lower right'); ax1C.ticklabel_format(axis = 'y', 
#           style='sci', scilimits=(0,0),useOffset=False,useLocale=True)

#handles, labels = hostA.get_legend_handles_labels()
#fig_fitness.legend(loc='upper left',bbox_to_anchor=(0.1,0.92,1.,0.05),ncol=4)

fig_fitness.savefig('fig2ABC.png', dpi = 300)
fig_fitness.savefig('Fig2ABC.svg')  

TempFood_test = True
if(TempFood_test == True):
    inds = 6; days_i=50
    wgt_i= np.zeros((inds,days_i*24)); age_i=np.zeros((inds,days_i*24)); stomach_i=np.zeros((inds,days_i*24))
    gut_i = np.zeros((inds,days_i,24)); len_i = np.zeros((inds,days_i*24))
    Temp_i = 26.
    for ind in range(0,inds):         # 1 is temperature-max, 2 onward is food lim
        a_count = 0; dwgt=0.018; stmch=0.
        for j in range(0,days_i-1):
         #Loop over 24 hours diel cycle
         Food_enc_day = 0.; resp_day=0.
         if (ind==4): day_j = 190 + j
         if (ind==5): day_j = 230 + j                   
         for h_time in range(24): 
             a_count += 1
             if dwgt>=EndSize: continue    #If flexion, done

             if ind==0:
                 dwgt += dwgt*(0.0418*Temp_i - 0.8355)/24.
                 stmch = MaxGut*dwgt*1.
             else:
                 GutLim = MaxGut*dwgt*1.
                 growthT = dwgt*SGR_T(Temp_i)
                 #the digestion is driven by growth, metabolism and assimilation efficiency
                 resp = metab_rate(dwgt,Temp_i)
                 digestion = (resp + growthT)/AE   
                 
                 #The gut is a container which determines if the larva can grow at max rate, limited by temperature or food
                 if ind==1: Food_enc = GutLim*1.0*LightLim[h_time,180]  #larva_feeding(180,h_time,wgt_i[ind]) #food encountered, in mgDW - from a given day, or 
                 if ind==2: Food_enc = GutLim*0.3*LightLim[h_time,180]; high_xval=a_count/24.  #larva_feeding(180,h_time,wgt_i[ind]) #food encountered, in mgDW - from a given day, or 
                 if ind==3: Food_enc = GutLim*0.2*LightLim[h_time,180]  #larva_feeding(180,h_time,wgt_i[ind]) #food encountered, in mgDW - from a given day, or 
                 if ind==4: 
                      Food_enc = larva_feeding(day_j,h_time,dwgt) #food encountered, in mgDW - from a given day, or 
                      growthT = dwgt*SGR_T(TempDay[day_j])
                       #the digestion is driven by growth, metabolism and assimilation efficiency
                      resp = metab_rate(dwgt,TempDay[day_j])
                      digestion = (resp + growthT)/AE
                 if ind==5: 
                      Food_enc = larva_feeding(day_j,h_time,dwgt) #food encountered, in mgDW - from a given day, or 
                      growthT = dwgt*SGR_T(TempDay[day_j])
                       #the digestion is driven by growth, metabolism and assimilation efficiency
                      resp = metab_rate(dwgt,TempDay[day_j])
                      digestion = (resp + growthT)/AE    

                             
                 # if the growth is temperature limited; enough food in the stomach 
                 if(stmch >= digestion):
                     dwgt  += growthT
                     stmch += (Food_enc - digestion)     #Add prey remove digestion
                     stmch  = min(stmch,GutLim) #Constrain gutcontent  #digestion is driven by growth
                     gut_i[ind,j,h_time] = stmch/GutLim
                 else:
                     dwgt  += (stmch*AE - resp)   #What is left in stomach is used to grow or cover respiration 
                     stmch  = Food_enc                 #Empty stomach, add feeding
                     stmch  = min(max(0.,stmch),GutLim)   #and it may fill up the stomach by feeding
                     gut_i[ind,j,h_time] = stmch/GutLim
            
             wgt_i[ind,a_count] = dwgt   #Do not shrink below hatch size  
             len_i[ind,a_count] = np.log(dwgt/0.0008)/905.2
             age_i[ind,a_count] = dh*a_count
             stomach_i[ind,a_count]  = stmch
             
                

if(TempFood_test == True):
    wgt_i[wgt_i == 0.] = np.nan; age_i[age_i == 0.] = np.nan

    fig_ind, ((ax1,ax2),(ax7,ax8),(ax3,ax4),(ax5,ax6)) = plt.subplots(4,2,
             gridspec_kw={'hspace': 0},figsize=(10,10),constrained_layout=True)
    fig_ind.suptitle('Single individuals', fontsize=16)
    
    ax2.plot(age_i[0,:],stomach_i[0,:],label='R18')
    ax2.plot(age_i[1,:],stomach_i[1,:],label='1.g/h')
    ax2.plot(age_i[2,:],stomach_i[2,:],label='.3g/h')
    #ax2.plot(age_i[3,:],stomach_i[3,:],label='0.2 guts/h')
    ax2.plot(age_i[4,:],stomach_i[4,:],label='d180')
    ax2.plot(age_i[5,:],stomach_i[5,:],label='d200')
    ax2.set_ylabel('Gut mass')
    ax2.set_xlabel('Age (days)')
    ax2.set_title('Stomach content')
    ax2.legend()
    
    ax1.set_title('Growth feeding stage')
    ax1.set_xlabel('Age (days)')
    ax1.set_ylabel('Body mass')
    ax1.plot(age_i[0,:],wgt_i[0,:],label='R18')
    ax1.plot(age_i[1,:],wgt_i[1,:],label='1 g/h')
    ax1.plot(age_i[2,:],wgt_i[2,:],label='.3 g/h')
    #ax1.plot(age_i[3,:],wgt_i[3,:],label='0.2 guts/h')
    ax1.plot(age_i[4,:],wgt_i[4,:],label='d180')
    ax1.plot(age_i[5,:],wgt_i[5,:],label='d200')
    ax1.legend()
    
    ax7.set_title('Length feeding stage')
    ax7.set_xlabel('Age (days)')
    ax7.set_ylabel('Body length in m'); ax7.set_xlim(left=0,right=high_xval)
    age_Malca = np.array([2,3,4,5,6,7,8,9,10,11,12,13,15,17,18,19])
    len_Malca = np.array([2.284,3.153,3.537,3.506,4.237,4.552,4.902,5.287,5.532,6.194,6.855,6.79,6.863,7.697,8.293,8.059])
    ax7.plot(age_Malca[:],0.001*len_Malca[:],'ko', label='M17')
    ax7.plot(2+age_i[0,:],len_i[0,:],label='R18')
    ax7.plot(2+age_i[1,:],len_i[1,:],label='1g/h')
    ax7.plot(2+age_i[2,:],len_i[2,:],label='.3g/h')
    #ax7.plot2+(age_i[3,:],len_i[3,:],label='0.2 guts/h')
    ax7.plot(2+age_i[4,:],len_i[4,:],label='d180')
    ax7.plot(2+age_i[5,:],len_i[5,:],label='d200')
    ax7.legend()
    
    ax8.set_title('Feeding stage obs')
    y2ax8 = ax8.twinx(); y2ax8.set_ylabel('DW in mg')
    y2ax8.plot(age_i[0,:],wgt_i[0,:],label='R18')
    #y2ax8.plot(age_i[0,:],0.0126*age_i[0,:]**(1.851),label='dwG06')
    y2ax8.plot(age_i[0,:],wgt_i[4,:],label='d180')
    y2ax8.plot(age_i[0,:],wgt_i[4,:],label='d200')

    ax8.set_xlabel('Age (days)')
    ax8.set_ylabel('Body length in mm'); ax8.set_xlim(left=0,right=high_xval)
    ax8.plot(age_Malca[:],len_Malca[:],'ko', label='M17')
    ax8.plot(2+age_i[0,:],2.85+0.35*age_i[0,:],'k-',label='G06_03')
    ax8.plot(2+age_i[0,:],1.65+0.43*age_i[0,:],'k--',label='G06_04')
    ax8.plot(2+age_i[0,:],1.86+0.41*age_i[0,:],'k-.',label='G06_05')
    ax8.plot(2+age_i[0,:],len_i[0,:]*1000,label='R18')
    ax8.plot(2+age_i[4,:],len_i[4,:]*1000,label='d180')
    ax8.legend(loc='upper left',ncol=2,prop={'size': 8});y2ax8.legend(loc='lower right')
    
    ax3.set_ylabel('Age (days'); ax5.set_ylabel('Age (days')
    ax5.set_xlabel('Hour of day'); ax6.set_xlabel('Hour of day')
    ax3.set_title('Gut, ind 1');ax4.set_title('Gut, ind 2')
    ax5.set_title('Gut, day 180');ax6.set_title('Gut, day 200')

    #ax3.contour(gut_i[2,:,:])
    #ax4.contour(gut_i[3,:,:])
    ax5.contour(gut_i[4,:,:])
    ax6.contour(gut_i[5,:,:])
    
    fig_ind.savefig('fig_ind.svg', dpi = 300)
    fig_ind.savefig('fig_ind.png', dpi = 300)
    
    plt.show(fig_ind)

"""figure of ration and temperature"""
#T_range=np.linspace(20,30)
T_range = np.arange(20,30, dtype='int')
ration_range = 10; w0 = 0.05
gTr =np.zeros((30,10)); dwgt =np.zeros((30,10)); dwgt[::]=w0
for Tp in T_range:
    for ration in range(ration_range):
        Food_enc_sum = 0.
        stmch = 0.0*w0
        for i_h in range(4,28):
             h_time = i_h
             if(i_h>23): h_time = i_h-24
             GutLim = MaxGut*dwgt[Tp,ration]
             growthT = dwgt[Tp,ration]*SGR_T(Tp)
              
             #The gut is a container which determines if the larva can grow at max rate, limited by temperature or food
             Food_enc = (ration/(ration_range*15))*w0*LightLim[h_time,180]  
             if(ration==0): Food_enc = larva_feeding(180,h_time,w0) #food encountered, in mgDW - from a given day, or 
             Food_enc_sum +=Food_enc
             #the digestion is driven by growth, metabolism and assimilation efficiency
             resp = metab_rate(w0,Tp)
             digestion = (resp + growthT)/AE   
                         
             # if the growth is temperature limited; enough food in the stomach 
             if(stmch >= digestion):
                 dwgt[Tp,ration]  += growthT
                 stmch += (Food_enc - digestion)     #Add prey remove digestion
                 stmch  = min(stmch,GutLim) #Constrain gutcontent  #digestion is driven by growth
             else:
                 dwgt[Tp,ration]  += (stmch*AE - resp)   #What is left in stomach is used to grow or cover respiration 
                 stmch  = min(max(0.,Food_enc),GutLim)   #and it may fill up the stomach by feeding
             #print(h_time,Food_enc,dwgt[Tp,ration])
        #print()
        Food_enc =(ration/(ration_range))*w0
        if(ration==0): Food_enc = larva_feeding(180,12,w0)*DayLength[180] #food encountered, in mgDW - from a given day, or 
        resp = metab_rate(w0,Tp)*24.
        w_max = w0*np.exp(0.0418*Tp - 0.8355)
        dwgt[Tp,ration] = min(w0+(Food_enc*AE - resp),w0*np.exp(0.0418*Tp - 0.8355))
        gTr[Tp,ration]= min((0.0418*Tp - 0.8355), (np.log(dwgt[Tp,ration]/w0)/1.))
       
#length_m = np.log(dwgt/0.0008)/905.2
fig_gTempRation_33, ax = plt.subplots(figsize=(8,4))

#gTr[20:30,0].tofile("gTr05_n500C015.npy") #Save different food concentrations
gTr_n500C015 = np.fromfile("gTr_n500C015.npy"); gTr_n400C01= np.fromfile("gTr_n400C01.npy")
gTr_n300C01= np.fromfile("gTr_n300C01.npy"); gTr_n300C015=np.fromfile("gTr_n300C015.npy")
ax.plot(T_range, (0.0418*T_range - 0.8355),label='R&al_18')
#ax.plot(T_range,gTr[20:30,1],label='0.1')
ax.plot(T_range,gTr[20:30,3],label='0.3')
ax.plot(T_range,gTr[20:30,6],label='0.6')
ax.plot(T_range,gTr[20:30,9],label='0.9')
ax.plot(T_range,gTr_n500C015,'--',label='Prey n500C015')
#ax.plot(T_range,gTr_n300C015,'--',label='Prey n300C015')
ax.plot(T_range,gTr_n400C01,'--',label='Prey n400C01')
ax.plot(T_range,gTr_n300C01,'--',label='Prey n300C01')

ax.set_title('Larval size 0.33 mgDW')
ax.set_xlabel("Temperature"); ax.set_ylabel("Growth rate g/g/day")
ax.legend(frameon=False)
plt.show(fig_gTempRation_33)
fig_gTempRation_33.savefig('fig_gTempRation_33.svg', dpi = 300)
fig_gTempRation_33.savefig('fig_gTempRation_33.png', dpi = 300)
#length_m = np.log(dwgt/0.0008)/905.2

#Plot for a different body size (0.5):
fig_gTempRation_05, ax = plt.subplots(figsize=(8,4))

#gTr[20:30,0].tofile("gTr05_n300C01.npy") #Save different food concentrations
gTr05_n500C015 = np.fromfile("gTr05_n500C015.npy"); gTr05_n400C01= np.fromfile("gTr05_n400C01.npy")
gTr05_n300C01= np.fromfile("gTr05_n300C01.npy"); gTr05_n300C015=np.fromfile("gTr05_n300C015.npy")
ax.plot(T_range, (0.0418*T_range - 0.8355),label='R&al_18')
#ax.plot(T_range,gTr[20:30,1],label='0.1')
ax.plot(T_range,gTr[20:30,3],label='0.3')
ax.plot(T_range,gTr[20:30,6],label='0.6')
ax.plot(T_range,gTr[20:30,9],label='0.9')
ax.plot(T_range,gTr05_n500C015,'--',label='Prey n500C015')
#ax.plot(T_range,gTr05_n300C015,'--',label='Prey n300C015')
ax.plot(T_range,gTr05_n400C01,'--',label='Prey n400C01')
ax.plot(T_range,gTr05_n300C01,'--',label='Prey n300C01')

ax.set_title('Larval size = 0.05 mgDW')
ax.set_xlabel("Temperature"); ax.set_ylabel("Growth rate g/g/day")
ax.legend(frameon=False)
plt.show(fig_gTempRation_33)
fig_gTempRation_05.savefig('fig_gTempRation_05.svg', dpi = 300)
fig_gTempRation_05.savefig('fig_gTempRation_05.png', dpi = 300)


fig_gDataField, ax8 = plt.subplots(figsize=(8,4))
ax8.set_title('Field data vs alternative prey abundances')
ax8.set_xlabel('Age (days)')
ax8.set_ylabel('Body length in mm'); ax8.set_xlim(left=0,right=high_xval)
ax8.plot(age_Malca[:],len_Malca[:],'ko', label='M17')
ax8.plot(2+age_i[0,:],2.85+0.35*age_i[0,:],'-',label='G06_03')
ax8.plot(2+age_i[0,:],1.65+0.43*age_i[0,:],'-',label='G06_04')
ax8.plot(2+age_i[0,:],1.86+0.41*age_i[0,:],'-',label='G06_05')
ax8.plot(2+age_i[0,:],len_i[0,:]*1000,label='R18')
#len_i[4,:].tofile("len_i_n500C015.npy") #Save different growth trajectories
len_i_n300C01=np.fromfile("len_i_n300C01.npy");len_i_n400C01=np.fromfile("len_i_n400C01.npy");len_i_n500C015=np.fromfile("len_i_n500C015.npy")
ax8.plot(2+age_i[4,:],len_i_n500C015*1000,'--',label='n500C015')
ax8.plot(2+age_i[4,:],len_i_n400C01*1000,'--',label='n400C01')
ax8.plot(2+age_i[4,:],len_i_n300C01*1000,'--',label='n300C01')
ax8.legend(loc='upper left',ncol=2,prop={'size': 10})
plt.show(fig_gDataField)



"""
fig_fitness2, ax = plt.subplots(figsize=(8,4))
ax1 = ax.twinx()
ax1.set_ylabel("Egg fitness")
ax1.ticklabel_format(axis = 'y',style='scientific')
ax.set_xlabel("Day of year")
#ax1.set_ymin(1E-6)
ax.set_xlim(150,250)
ax.set_ylim(17,29)
ax.set_ylabel("Temperature")
#ax.plot(a2[180,1:1200],bm2[180,1:1200],label='180')
temp2003=np.fromfile('temp2003.npy'); temp2004=np.fromfile('temp2004.npy')
temp2011=np.fromfile('temp2011.npy');survT2011 = np.fromfile('survT2011.npy')
survT2004=np.fromfile('survT2004.npy');survT2003 = np.fromfile('survT2003.npy')
ax.plot(range(1,365),TempDay[1:365],'r-', label='Temp Avg')
ax.plot(range(1,365),temp2003[1:365],'r-.', label='Temp 2003')
ax.plot(range(1,365),temp2004[1:365],'r--', label='Temp 2004')
ax.plot(range(1,365),temp2011[1:365],'r', dashes=[1, 2],label='Temp 2011')
#ax1.set_yscale('log')
ax.legend(loc = 'best')
ax1.plot(Year,surv,'k-', label ='Surv Avg')
ax1.plot(Year,survT2003,'k-.', label ='Surv 2003')
ax1.plot(Year,survT2004,'k--', label ='Surv 2004')
ax1.plot(Year,survT2011,'k', dashes=[1, 2],label ='Surv 2011')
ax1.legend(loc = 'lower right')

#ax1.plot(Year,tot_egg_survMcG,'y-.', label ='Total egg survival McG')
#handles, labels = ax.get_legend_handles_labels()
#fig_fitness2.legend(bbox_to_anchor=(0.1, 0.95, 0.8, 0.102), loc='lower left',
#ncol=2, mode="expand", borderaxespad=0.)
plt.show()
"""
