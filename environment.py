# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 15:24:27 2019

@author: nfiof
"""

def TempandClad():
    import numpy as np

    TempDay = np.zeros(365, float)
    CladoDay = np.zeros(365, float)

    #Temperature and zooplankton measures, midmonth values. 
    MonthDays = [30,30,31,30,31,30,31,30,31,30,31,30]
    TempMidMonth = [14.8049,14.1126,14.1915, 15.4093,17.9667,21.6733,24.958,
                26.6576,25.8626, 23.2371,20.0049,16.7536]
    CladoceraMidMonth= [10,20,40,94,141,198,262,171,99,70,17,15]

    for i in range(0,13):
        if i==0:
            for Jday in range(1,16):
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
            for j in range(0,15):
                TempDay[Jday] = TempMidMonth[i-1] + j*DayIncr
                CladoDay[Jday] = CladoceraMidMonth[i-1] + j*DayIncrC
                Jday += 1

        return
    
    
#plot environmental variables temperature and food over the season
fig1, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax3 = ax1.twinx()
major_ticks = np.arange(0, 365, 30)
ax1.set_xlabel('Day of year')
ax1.set_xticks(major_ticks)
ax1.set_ylabel('Temperature')
ax2.set_ylabel('Cladocera')
ax1.plot(range(1,365),TempDay[1:365],'b-', label='Temp')
ax1.plot(range(15,350,30), TempMidMonth[0:12], 'ro', label='Midmonth Temp')
ax2.plot(range(15,350,30),CladoceraMidMonth[0:12],'go', label='Midmonth C')
ax2.plot(range(1,365),CladoDay[1:365],'g-', label='Cladocera')
ax3.plot(range(1,365),DayLength[1:365],'y-', label='Daylength')
ax1.legend(loc='best')
ax2.legend(loc='best')
ax3.legend(loc='lower center')
plt.title('Data on temperature and Cladocerans')
plt.show()     