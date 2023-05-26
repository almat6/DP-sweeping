'''
This code should start first before the IR data collecting,
This code is used to monitor the data created from IR characterization, and copy the rawdata to S_drive, then get the peak area of monomer from IR spectrum
convert peak area into concentration and conversion of monomer.  
'''


import csv
from posixpath import split, splitext
from time import sleep
from turtle import color
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy as scipy
from scipy import optimize
import matplotlib.ticker
from matplotlib import gridspec
from matplotlib import rcParams
import matplotlib.ticker as ticker
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import os
import glob
import datetime
import time
from xdrlib import ConversionError
import watchdog
import pandas as pd
from numpy import sqrt
from pathlib import Path
import sys 
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import logging
import shutil
import matplotlib.ticker
from numpy.linalg import norm
import math


# create file to save experiment data
print("Please input the name of you Experiment :>>")
ExperimentName = input(':>>')
print("Please input the parentpath of saving your data")
ParentFolder = r"S:\Sci-Chem\PRD\IR 112\Ansila" #input(':>>')

# Obtain the year,month,date
today = datetime.datetime.now()
YearTime = today.strftime("%Y")
MonthTime = today.strftime("%m")
day = today.strftime("%d")
print(YearTime,MonthTime,day)

#Create Year folder
YearFolderPath = r'{}\{}'.format(ParentFolder,YearTime)
if not os.path.exists(YearFolderPath):
    os.makedirs(YearFolderPath)

#Create Month folder under year folder
MonthFolderPath = r'{}\{}\{}'.format(ParentFolder,YearTime, MonthTime)
if not os.path.exists(MonthFolderPath):
    os.makedirs(MonthFolderPath)

#Create day folder under Month
DayFolderPath = r'{}\{}\{}\{}'.format(ParentFolder, YearTime, MonthTime,day)
if not os.path.exists(DayFolderPath):
    os.makedirs(DayFolderPath)

# Create Experiment folder under day and create subfolders under experiment to save all data
ExperimentNameFolder = r'{}\{}\{}\{}\{}'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName)
ExperimentNameFolder_1 = r'{}\{}\{}\{}\{}_1'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName)
ExperimentNameFolder_2 = r'{}\{}\{}\{}\{}_2'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName)
ExperimentNameFolder_3 = r'{}\{}\{}\{}\{}_3'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName) 
        
if not os.path.exists(ExperimentNameFolder):
    os.makedirs(ExperimentNameFolder)
    os.makedirs(r'{}\{}\{}\{}\{}\{}'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName,'IR_RawData'))
    os.makedirs(r'{}\{}\{}\{}\{}\{}'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName,'DeconvolutionPicture'))
    os.makedirs(r'{}\{}\{}\{}\{}\{}'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName,'RealtimePicture'))
    path_t = ExperimentNameFolder
        
elif not os.path.exists(ExperimentNameFolder_1):
    os.makedirs(ExperimentNameFolder_1)
    os.makedirs(r'{}\{}\{}\{}\{}_1\{}'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName,'IR_RawData'))
    os.makedirs(r'{}\{}\{}\{}\{}_1\{}'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName,'DeconvolutionPicture'))
    os.makedirs(r'{}\{}\{}\{}\{}_1\{}'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName,'RealtimePicture'))
    path_t = ExperimentNameFolder_1

elif not os.path.exists(ExperimentNameFolder_2):
    os.makedirs(ExperimentNameFolder_2)
    os.makedirs(r'{}\{}\{}\{}\{}_2\{}'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName,'IR_RawData'))
    os.makedirs(r'{}\{}\{}\{}\{}_2\{}'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName,'DeconvolutionPicture'))
    os.makedirs(r'{}\{}\{}\{}\{}_2\{}'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName,'RealtimePicture'))
    path_t = ExperimentNameFolder_2

else:
    os.makedirs(ExperimentNameFolder_3)
    os.makedirs(r'{}\{}\{}\{}\{}_3\{}'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName,'IR_RawData'))
    os.makedirs(r'{}\{}\{}\{}\{}_3\{}'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName,'DeconvolutionPicture'))
    os.makedirs(r'{}\{}\{}\{}\{}_3\{}'.format(ParentFolder,YearTime, MonthTime,day,ExperimentName,'RealtimePicture'))
    path_t = ExperimentNameFolder_3


# concentration calculation 
def CalConcentration(peakarea)->float:
    '''
    Explanation: this function is used to convert the peak area of monomer into its concentration
    '''

    a = -0.0338
    b = 1.1486
    d = 0.2009
    c = d - float(peakarea)

    delta = b ** 2 - 4 * a * c
    x = (b - math.sqrt(delta)) / (-2 * a)

    return(x)

# conversion calculation
def CalConversion(x:float, y:float)->float:
    '''
    Explanation: This function is used to calculate the conversion of monomer:
    x: the real time concentration of Monomer  
    y: the initial concentration of Monomer
    '''
    Conversion = (1 - (x/y))*100
    return(Conversion)

# Area under curve integrate function 
def integrate(x, y):
    area = np.trapz(y=y, x=x)
    return area


# create a empty csv file to save all the calsulated data 
with open(r'{}\{}-Data.csv'.format(path_t, ExperimentName), 'w') as f:
    pass 

# cwd = os.getcwd()
# print(cwd)

# Monomer Initial Concentration calculation and data collect time calculation 
V_reactor = 6 # float(input('Please input the Volume of your reactor (ml):>> '))
V_dead = 0.1325 #float(input('\033[1;31mPlease input the dead Volume of your reactor (ml):>> \033[0m'))
V_sample = 0.05 #float(input('Please input the Volume between the sample collector and equipment (ml):>> '))
V_input = 0.1325 #float(input('Please input the Volume between mixer and the reactor (ml):>> '))
StockCon_Monomer = 6 #float(input('\033[1;31mPlease input the initial concentration of your Monomer(M):>> \033[0m'))
StockCon_RAFTAgent = 0.5 #float(input('\033[1;31mPlease input the initial concentration of your RaftAgent (M):>> \033[0m'))

# print(f'V_reactor is {V_reactor} ml,\nV_dead is {V_dead} ml,\nC_monomer is {C_monomer} M,\nC_RaftAgent is {C_RaftAgent} M, \nC_Initiator is {C_Initiator} M')
Pro_Monomer = 3 #float(input('\033[1;31mPlease input the proposed  concentration of Monomer:>> \033[0m'))
#C_Monomerend = float(input('\033[1;31mPlease input the proposed end concentration of Monomer:>> \033[0m'))

DP_Monomerstart = 50 #float(input('\033[1;31mPlease input the proposed start DP of Monomer:>> \033[0m'))
DP_Monomerend = 170 #float(input('\033[1;31mPlease input the proposed end DP of Monomer:>> \033[0m'))

# # input the equivlent of reagents
# equivlent_Monomer = float(input("Please input the equivlent of your monomer:>>  "))
# equivlent_RaftAgent = float(input("Please input the equivlent of your RaftAgent:>>  "))

# input the residence time of concentration sweep
ResidenceTime = float(input('Please input the involved residence time (seconds):>> '))
DPchangetime = 600
sleeptime = 5 #float(input('Please input the pump sleep time (seconds):>> '))
Flowrate = V_reactor* 60 / ResidenceTime #converting it to ml/min
cleartime = (V_dead+V_input+V_reactor)*60/Flowrate + 300 #extra 5 minutes to make sure the old solution in the flow reactor is cleared out

# FlowRateMonomer = []
# FlowRateRaftAgent = []
# FlowRateSolvent = []
# Dp = []
# ConcentrationMonomer = []
#
# FRMonomer = Flowrate * C_Monomerstart / StockCon_Monomer
# FlowRateMonomer.append(FRMonomer)
# FRRaftAgent_0 = Flowrate* C_Monomerstart * equivlent_RaftAgent/ (equivlent_Monomer * StockCon_RAFTAgent)
# FlowRateRaftAgent.append(FRRaftAgent_0)
# FRSolvent_0 = Flowrate - FRMonomer - FRRaftAgent_0
# FlowRateSolvent.append(FRSolvent_0)
# Dp_1 = StockCon_Monomer * FRMonomer / Flowrate / (StockCon_RAFTAgent * FRRaftAgent_0 / Flowrate)
# Dp.append(Dp_1)
#
# n = int(ResidenceTime/sleeptime)
# print(n)
# ConDecreaseStep = (C_Monomerstart - C_Monomerend)/n
# VoluDecreaseStep = ConDecreaseStep * Flowrate / StockCon_Monomer
# print(ConDecreaseStep, VoluDecreaseStep)
#
# for i in range(n):
#     FRMonomer = FRMonomer - VoluDecreaseStep
#     ConMonomer = FRMonomer * StockCon_Monomer / Flowrate
#     FRRaftAgent = ConMonomer*equivlent_RaftAgent*Flowrate/(equivlent_Monomer * StockCon_RAFTAgent)
#     FRSolvent = Flowrate - FRMonomer - FRRaftAgent
#     Dp_1 = StockCon_Monomer * FRMonomer / Flowrate / (StockCon_RAFTAgent * FRRaftAgent / Flowrate)
#     FlowRateMonomer.append(FRMonomer)
#     FlowRateRaftAgent.append(FRRaftAgent)
#     FlowRateSolvent.append(FRSolvent)
#     ConcentrationMonomer.append(ConMonomer)
#     Dp.append(Dp_1)
#
# ConcentrationMonomer.reverse()
#
# RAFT agent con is changed, Monomer con is the same
FlowRateMonomer = []
ConcentrationMonomer = []
FlowRateRAFTAgent = []
FlowRateSolvent = []
DegP = []

No_steps = int(DPchangetime / sleeptime)
DPDecreaseStep = abs((DP_Monomerstart - DP_Monomerend)) / No_steps
FRMonomer = Pro_Monomer * Flowrate / StockCon_Monomer
#ConMonomer = FRMonomer * StockCon_Monomer / Flowrate
ConcentrationMonomer.append(Pro_Monomer)

for i in range(No_steps+1):
    dp = DP_Monomerstart + i * DPDecreaseStep
    FRRaftAgent = (StockCon_Monomer * FRMonomer / dp) / StockCon_RAFTAgent
    FRSolvent = Flowrate - FRMonomer - FRRaftAgent
    FlowRateMonomer.append(FRMonomer)
    DegP.append(dp)
    FlowRateRAFTAgent.append(FRRaftAgent)
    FlowRateSolvent.append(FRSolvent)


print(f"the lenth of Flowratemonomer is {len(FlowRateMonomer)}")
# create a text file to save all the parameters of reaction
Reactionparameter = [
                    f"The name of this Experiment is {ExperimentName}",
                    f"The volume of reactor is {V_reactor} ml.",
                    f"The dead Volume of reactor is {V_dead} ml.",
                    f"The volume between the mixer and reator is {V_input} ml.",
                    f"The volume between the equipment dector and sample collector is {V_sample} ml.",
                    f"The initial concentraton of monomer is {StockCon_Monomer} M",
                    f"The initial concentration of raftgent is {StockCon_RAFTAgent} M",
                    #f"The proposed start concentration in the end of time sweep is {C_Monomerstart} M",
                    #f"The proposed end concentration in the end of time sweep is {C_Monomerend} M",
                    #f"The proposed equivalent of monomer in the end of time sweep is {equivlent_Monomer}",
                    #f"The proposed equivalent of raftagent in the end of time sweep is {equivlent_RaftAgent}",
                    f"The residence time involved in this reaction is {ResidenceTime} second",
                    f"The sleep time of pump involved in this reaction is {sleeptime} second",

]

ParameterFile = open(r'{}/ExperimentParameter.txt'.format(path_t), 'a') 
ParameterFile.writelines('\n'.join(Reactionparameter))
ParameterFile.close()
sleep(0.5)


PeakAreall = []
PeakArea = []
ScanTime = []
InitialConcentration = []
InitialDP_list = []
Concentration = []
Conversion = []
DP = []


class Handler(FileSystemEventHandler):
    def on_created(self, event):
        file_created = event.src_path
        print(file_created)
        filename = os.path.splitext(os.path.split(file_created)[1])[0]
        listfile_created = os.listdir(my_path)
        
        if len(listfile_created)*5 >= (V_dead+V_input+V_reactor)*60/Flowrate + 120:
            try:
                shutil.copy(file_created, rawdata_path)
                sleep(0.5)

                with open(file_created) as csvfile:
                    rawdata = list(csv.reader(csvfile, delimiter = ","))
                
                # set scan time 
                listfile = os.listdir(rawdata_path)
                j = len(listfile)
                scantime = j*5/60 
                ScanTime.append(scantime)
                print(f"{j} scan data collected ")


                # print(rawdata)
                exampledata = np.array(rawdata[1:],dtype = np.float64) # rawdata[1:], first line is the head, so we start from the second line
                # print(exampledata)
                data_x = exampledata[:,0]
                x_data = list(data_x)
                xdata = []

                for i in data_x:
                    if 1588 <= i <= 1652:
                        xdata.append(i)
                # print(xdata)
                indexlow = x_data.index(xdata[0])
                indexhigh = x_data.index(xdata[-1])
                # print(indexhigh, indexlow)
                ydata_1 = exampledata[indexlow:indexhigh+1,1]
                
                x_base = [1652, 1588]
                y_base = [ydata_1[0], ydata_1[-1] ]

                Allarea = abs(integrate(xdata, ydata_1))
                Area2 = abs(integrate(x_base,y_base))
                peakarea = Allarea - Area2

                # print(Allarea, Area2, peakarea)
                PeakArea.append(peakarea)
                # PeakAreall.append(Allarea)

                # Calculate concentration and conversion
                concentration = CalConcentration(peakarea)
                Concentration.append(concentration)
                InitialConcentration.append(Pro_Monomer)

                T1 = 180 #to stablise the curve
                cleartime_filenum = int((T1/5))
                DPindex = int(5/sleeptime) #5 seconds
                print (f"1. {cleartime_filenum}, 2. {DPindex}")
                # change this part to calculate DP
                if j*5 <= T1:
                    InitialDP = DegP[0]
                    DP.append(InitialDP)

                    # InitialCon = ConcentrationMonomer[0]
                    # InitialConcentration.append(InitialCon)
                    # DP_1 = Dp[0]
                    # DP.append(DP_1)



                elif T1 < j*5 < T1 + DPchangetime:
                    index = int( DPindex * (j-cleartime_filenum))
                    print(f"index is {index}")
                    InitialDP = DegP[index]
                    DP.append(InitialDP)

                    # InitialCon = ConcentrationMonomer[index]
                    # InitialConcentration.append(InitialCon)
                    # DP_1 = Dp[index]
                    # DP.append(DP_1)

                else:
                    InitialDP = DegP[-1]
                    DP.append(InitialDP)

                    # InitialCon = ConcentrationMonomer[-1]
                    # InitialConcentration.append(InitialCon)
                    # DP_1 = DegP[-1]
                    # DP.append(DP_1)

                conversion = CalConversion(concentration,Pro_Monomer) #intialcon will be constant
                Conversion.append(conversion) 


                # save the calculated data to a new CSV file under ExperimentName folder 
                data = {
                    'scantime/minute':ScanTime,
                    'Peakarea':PeakArea,
                    'InitialConcentration':InitialConcentration,
                    'DP':DP,
                    'Concentration/M': Concentration,
                    'Conversion/%':Conversion,
                }
                column_names = ['scantime/minute','Peakarea','InitialConcentration','DP','Concentration/M', 'Conversion/%']
                df = pd.DataFrame(data, columns = column_names) #columns = column_names
                df.to_csv(r'{}\{}-Data.csv'.format(path_t,ExperimentName), columns = column_names)

                fig = plt.figure()
                ax = fig.add_subplot(111)
                for axis in ['top', 'bottom', 'left', 'right']:
                    ax.spines[axis].set_linewidth(2) 
                ax.tick_params(axis = 'both', which = 'major', labelsize = 12)
                ax.tick_params(axis = 'both', which = 'minor', labelsize = 12)
                plt.plot(ScanTime, PeakArea, color = 'red', linewidth = 4, linestyle = ':')
                plt.xlabel('Scantime/minute',fontsize = 14)
                plt.ylabel('Peakarea',fontsize = 14)
                plt
                plt.savefig(f'{path_t}/RealtimePicture/Scantime_Peakarea')
                plt.clf()
                plt.close()
                
                fig = plt.figure()
                ax = fig.add_subplot(111)
                for axis in ['top', 'bottom', 'left', 'right']:
                    ax.spines[axis].set_linewidth(2) 
                ax.tick_params(axis = 'both', which = 'major', labelsize = 12)
                ax.tick_params(axis = 'both', which = 'minor', labelsize = 12)
                plt.plot(ScanTime, Concentration, color = 'blue', linewidth = 3, linestyle = 'dotted')
                plt.plot(ScanTime, InitialConcentration, color = 'red', linewidth = 3, linestyle = 'dotted')
                plt.legend(["real time", "Initial"])
                plt.xlabel('Scantime/minute', fontsize = 14)
                plt.ylabel('Concentration',fontsize = 14)
                plt.savefig(f'{path_t}/RealtimePicture/Scantime_Concentration')
                plt.clf()
                plt.close()
                
                fig = plt.figure()
                ax = fig.add_subplot(111)
                for axis in ['top', 'bottom', 'left', 'right']:
                    ax.spines[axis].set_linewidth(2) 
                ax.tick_params(axis = 'both', which = 'major', labelsize = 12)
                ax.tick_params(axis = 'both', which = 'minor', labelsize = 12)
                plt.plot(ScanTime, Conversion, color = 'black', linewidth = 3, linestyle = '-')
                plt.xlabel('Scantime/minute', fontsize = 14)
                plt.ylabel('Conversion', fontsize = 14)
                plt.savefig(f'{path_t}/RealtimePicture/Scantime_Conversion')
                plt.clf()
                plt.close()

            except Exception:
                pass
        else:
            pass

my_path = r"C:\Users\IR112\Documents\iC IR Experiments\{}".format(ExperimentName) #folder we monitor
rawdata_path = r'{}/IR_RawData'.format(path_t) #folder for rawdata


observer = Observer()
event_handler = Handler()
observer.schedule(event_handler, my_path, recursive = True)
print('Observer start')
observer.start()

try:
    while True:
        sleep(0.05)
except KeyboardInterrupt:
    observer.stop()
observer.join()
