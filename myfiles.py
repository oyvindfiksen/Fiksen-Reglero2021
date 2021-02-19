# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 20:34:41 2020

@author: nfiof
"""
import Tunalarvaseason

fig_1A, (ax1,ax3) = plt.subplots(2,sharex='col',sharey=False,
         gridspec_kw={'hspace': 0},figsize=(8,10),constrained_layout=False)
ax3.set_xlabel("Day of year")
ax1.set_ylabel('Prey enc /W /day')
ax3.set_ylabel('Prey encounter' '\n' 'and energy surplus')
ax3.set_ylim(-0.1, ymax=1.)
ax1.set_ylim(0,1)
ax3.plot(Year,(sgr+SMR)/AE,label ='(SGR+SMR)/AE')
ax3.plot(Year,food_intake_day/0.77, label ='I/W')
ax3.plot(Year,((food_intake_day/0.77)-((sgr+SMR)/AE)), label ='I/W-(SGR+SMR)/AE')
ax3.plot(Year,((food_intake_day/0.77)-(SMR/AE)), label ='I/W-(SMR/AE)')

ax1.plot(Year,food_intake_day/wgt, label ='W=0.77',linewidth=3)
ax1.plot(Year,food_intake_day24/wgt, label ='W=0.77 24h light')
ax1.plot(Year,Food_intake_day_wgt050/0.5, label ='W=0.5')
ax1.plot(Year,Food_intake_day_wgt010/0.1, label ='W=0.1')

ax1.legend(loc="best"); ax2.legend(loc="best"); ax3.legend(loc="best") 
fig_1A.savefig('fig_feeding.png', dpi = 300)
fig_1A.savefig('fig_feeding.svg', dpi = 300)

plt.show()    

 
