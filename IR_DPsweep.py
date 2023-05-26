from datetime import datetime
from turtle import clear
import serial
import serial.tools.list_ports
from time import sleep
from time import *
from syringepump import SyringePump


# Find the COMs that syringe pumps connected to computer
PortData = serial.tools.list_ports.comports()
print(PortData)

for port in PortData:
    print(f"\033[1;31m{port}\033[0m")


# Pumps Assigned
# PumpMonomer = SyringePump('COM3', 'PumpMonomer')
# PumpRaftAgent = SyringePump('COM6', 'PumpRaftAgent')
# PumpSolvent = SyringePump('COM7','PumpSolvent')

# Experiment condition input 
V_reactor = 0.9 #float(input('Please input the Volume of your reactor (ml):>> '))
V_dead = 0 #float(input('\033[1;31mPlease input the dead Volume of your reactor (ml):>> \033[0m'))
V_sample = 0 #float(input('\033[1;31mPlease input the Volume between the sample collector and the equipment (ml):>> \033[0m'))
V_input = 0 #float(input('\033[1;31mPlease input the Volume between the mixer and the reactor (ml):>> \033[0m'))
StockCon_Monomer = float(input('\033[1;31mPlease input the initial concentration of your Monomer(M):>> \033[0m'))
StockCon_RAFTAgent = float(input('\033[1;31mPlease input the initial concentration of your RaftAgent (M):>> \033[0m'))

# print(f'V_reactor is {V_reactor} ml,\nV_dead is {V_dead} ml,\nC_monomer is {C_monomer} M,\nC_RaftAgent is {C_RaftAgent} M, \nC_Initiator is {C_Initiator} M')

Pro_Monomer = float(input('\033[1;31mPlease input the propsed concentration of Monomer:>> \033[0m'))
# C_Monomerend = float(input('\033[1;31mPlease input the propsed end concentration of Monomer:>> \033[0m'))

# # input the equivlent of reagents
# equivlent_Monomer = float(input("Please input the initial equivlent of your monomer:>>  "))
# equivlent_RaftAgent = float(input("Please input the initial equivlent of your RaftAgent:>>  "))

DP_Monomerstart = float(input('\033[1;31mPlease input the proposed start DP of Monomer:>> \033[0m'))
DP_Monomerend = float(input('\033[1;31mPlease input the proposed end DP of Monomer:>> \033[0m'))

# input the residence time of timesweeps
ResidenceTime = float(input('Please input the involved residence time (seconds):>> '))
ConcentrationChangetime = 600
sleeptime = 5 #float(input('Please input the pump sleep time (seconds):>> '))

DPchangetime = 600
Flowrate = V_reactor* 60 / ResidenceTime #converting it to ml/min
cleartime = (V_dead+V_input+V_reactor)*60/Flowrate + 300 #extra 5 minutes to make sure the old solution in the flow reactor is cleared out

Flowrate = V_reactor* 60 / ResidenceTime



# FlowRateMonomer = []
# FlowRateRaftAgent = []
# FlowRateSolvent = []
# DP = []
#
# FRMonomer =  Flowrate*C_Monomerstart / C_Monomer
# FlowRateMonomer.append(FRMonomer)
# FRRaftAgent_0 = Flowrate* C_Monomerstart * equivlent_RaftAgent/ (equivlent_Monomer * C_RaftAgent)
# FlowRateRaftAgent.append(FRRaftAgent_0)
# FRSolvent_0 = Flowrate - FRMonomer - FRRaftAgent_0
# FlowRateSolvent.append(FRSolvent_0)
# DP_1 = C_Monomer*FRMonomer/Flowrate/(C_RaftAgent*FRRaftAgent_0/Flowrate)
# DP.append(DP_1)

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
 
n = int(DPchangetime/sleeptime)

#work out the calculated DP to make sure the code is correct
# CalDP = []
# for i in range(No_steps):
#     DP = FRMonomer * StockCon_Monomer / (FlowRateRAFTAgent[i] * StockCon_RAFTAgent)
#     CalDP.append(DP)
#
# print(CalDP)

# ConDecreaseStep = (C_Monomerend - C_Monomerstart)/n
# VoluDecreaseStep = ConDecreaseStep*Flowrate/C_Monomer
# print(ConDecreaseStep, VoluDecreaseStep)
# for i in range(n+1):
#     FRMonomer = FRMonomer + VoluDecreaseStep
#     ConMonomer = FRMonomer * C_Monomer/Flowrate
#     FRRaftAgent = ConMonomer*equivlent_RaftAgent*Flowrate/(equivlent_Monomer * C_RaftAgent)
#     FRSolvent = Flowrate - FRMonomer - FRRaftAgent
#     DP_1 = C_Monomer*FRMonomer/Flowrate/(C_RaftAgent*FRRaftAgent/Flowrate)
#     FlowRateMonomer.append(FRMonomer)
#     FlowRateRaftAgent.append(FRRaftAgent)
#     FlowRateSolvent.append(FRSolvent)
#     DP.append(DP_1)


# FlowRateMonomer.reverse()
print("1\n")
print(FlowRateMonomer)
print("2\n")
print(FlowRateRAFTAgent)
print("3\n")
print( FlowRateSolvent)
# FlowRateRaftAgent.reverse()
# FlowRateSolvent.reverse()

# print(
#     len(FlowRateMonomer),
#     f"\nthe flow rate of \033[1;31mMonomer\033[0m involved is {FlowRateMonomer}\n", 
#     len(FlowRateRaftAgent),
#     f"\nthe flow rate of \033[1;31mRaftagent\033[0m involved is {FlowRateRaftAgent}\n",
#     len(FlowRateSolvent),
#     f"\nthe flow rate of \033[1;31mSolvent\033[0m involved is {FlowRateSolvent}\n",
#     len(DP),
#     f"\nthe \033[1;31mdegree of polymerization\033[0m involved is {DP}\n",    
#     )



# calculate sleep time 

ClearTime = (V_dead+V_input+V_reactor)*60/Flowrate +300
SleepTime = []

SleepTime.append(180)

for i in range(n-1):
    Sleep = sleeptime
    SleepTime.append(Sleep)

SleepTime.append(ClearTime)
print(len(SleepTime))
print(f"The sleep time of pumps involved is {SleepTime}")

Totalflowmonomer = 0 
TotalflowRaft= 0 
TotalflowSolvent = 0 

for i in range(len(SleepTime)):
    Totalflowmonomer += SleepTime[i]/60*FlowRateMonomer[i] 
    TotalflowRaft += SleepTime[i]/60*FlowRateRAFTAgent[i]
    TotalflowSolvent += SleepTime[i]/60*FlowRateSolvent[i]  

print(
    f"the volume of Monomer needed in this experiment is {Totalflowmonomer}", 
    f"\nthe volume of Raftagent needed in this experiment is {TotalflowRaft}", 
    f"\nthe volume of Solvent needed in this experiment is {TotalflowSolvent}", 
)




# PumpSolvent.start(), PumpMonomer.start(),PumpRaftAgent.start()
# sleep(0.5)
# PumpSolvent.changeFlowrate(0),  PumpMonomer.changeFlowrate(0), PumpRaftAgent.changeFlowrate(0)
# sleep(0.5)

# # Start experiment
# for i in range(n+1):
#     # get current time
#     today = datetime.now()
#     CurrentTime = today.strftime("%H:%M:%S")
#     print(f'pumps start at {CurrentTime}')
#     PumpMonomer.changeFlowrate(FlowRateMonomer[i])
#     PumpSolvent.changeFlowrate(FlowRateSolvent[i])
#     PumpRaftAgent.changeFlowrate(FlowRateRAFTAgent[i])
#     sleep(SleepTime[i])
#
#
# PumpSolvent.stop(), PumpMonomer.stop(), PumpRaftAgent.stop()

