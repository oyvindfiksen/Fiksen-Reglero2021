# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 07:55:05 2019

@author: nfiof
"""

# importing pandas package
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import datetime
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

 # read field data of larvae
df = pd.read_csv("larval_data_bluefin_2001_2013.csv", sep=';')
dfT = pd.read_excel("temperature_data.xlsx", sep=';')
pd.DataFrame(dfT)
pd.DataFrame(df)
plt.plot('date','meansst', data=dfT)

dato = dfT['date'] #raw_input("Enter date: ")  ## format is 02-02-2016
#dato=replace('.','-', inplace = True)

dfT['dayofyear'] = dfT['date'].dt.dayofyear
dfT['year'] = dfT['date'].dt.year

day_of_year = dfT['dayofyear']
year = dfT['year']
temp = dfT['meansst']
temp2003 = np.array(temp[2*366:3*366])
temp2004 = np.array(temp[3*366:4*366])
temp2005 = np.array(temp[4*366:5*366])
temp2006 = np.array(temp[5*366:6*366])
temp2007 = np.array(temp[6*366:7*366])
temp2011 = np.array(temp[2554:2920])


plt.plot(range(0,366),temp2003,label='2003')
plt.plot(range(0,366),temp2004, label='2004')
plt.plot(range(0,366),temp2005,label='2005')
plt.plot(range(0,366),temp2006,label='2006')
plt.plot(range(0,366),temp2007, label='2007')
plt.plot(range(0,366),temp2004, label='2004')
plt.plot(range(0,366),temp2003,label='2003')
plt.plot(range(0,366),temp2011,'k-', label='2011')
plt.plot(range(0,366),TempDay, label='Average')
plt.scatter(range(15,350,30),TempMidMonth)
plt.legend()

plt.plot(day_of_year[1*366:2*366],temp[1*366:2*366])
plt.plot(day_of_year[3*366:4*366],temp[3*366:4*366])

temp2003.tofile('temp2003.npy')
temp2004.tofile('temp2004.npy')

str_time= "2018-06-03 08:00:00"
date_date = datetime.datetime.strptime("date", "%Y-%m-%d")
dfT['date']


"""
df.replace(-999,np.nan, inplace = True)

yolksac=df["n_thyn_yolk"]
temp = df["Temperature_mixed_layer"]

plt.scatter(temp,yolksac)
plt.ylim(1,500)
plt.yscale('log')

print(max(temp))
sns.distplot( df["n_thyn_pre"] )
sns.distplot( df["JD"] )
sns.kdeplot(df['JD'], shade=True, bw=2., color="olive")
plt.stem('JD', 'Lnm3', use_line_collection= True, data=df)
plt.yscale('log')
"""    