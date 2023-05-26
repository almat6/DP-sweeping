from datetime import datetime
from turtle import clear
import serial
import serial.tools.list_ports
from time import sleep
from time import *
from syringepump import SyringePump


# # Find the COMs that syringe pumps connected to computer
# PortData = serial.tools.list_ports.comports()
# print(PortData)

# for port in PortData:
#     print(f"\033[1;31m{port}\033[0m")


#Pumps Assigned
PumpMonomer = SyringePump('COM7', 'PumpMonomer')
PumpRaftAgent = SyringePump('COM6', 'PumpRaftAgent')
PumpSolvent = SyringePump('COM5','PumpSolvent')
PumpInitiator = SyringePump('COM8', 'PumpInitiator')

# Experiment condition input 
V_reactor = 6 #float(input('Please input the Volume of your reactor (ml):>> '))
V_dead = 0.1325 #float(input('\033[1;31mPlease input the dead Volume of your reactor (ml):>> \033[0m'))
V_sample = 0.05 #float(input('\033[1;31mPlease input the Volume between the sample collector and the equipment (ml):>> \033[0m'))
V_input = 0.1705 #float(input('\033[1;31mPlease input the Volume between the mixer and the reactor (ml):>> \033[0m'))
StockCon_Monomer = 6 #float(input('\033[1;31mPlease input the initial concentration of your Monomer(M):>> \033[0m'))
StockCon_RAFTAgent = 0.372 #float(input('\033[1;31mPlease input the initial concentration of your RaftAgent (M):>> \033[0m'))
StockCon_Initiator = 0.05
StockCon_Monomer_inRA = 0.28 

# print(f'V_reactor is {V_reactor} ml,\nV_dead is {V_dead} ml,\nC_monomer is {C_monomer} M,\nC_RaftAgent is {C_RaftAgent} M, \nC_Initiator is {C_Initiator} M')

Pro_Monomer = 3 #float(input('\033[1;31mPlease input the propsed concentration of Monomer:>> \033[0m'))
# C_Monomerend = float(input('\033[1;31mPlease input the propsed end concentration of Monomer:>> \033[0m'))

# # input the equivlent of reagents
# equivlent_Monomer = float(input("Please input the initial equivlent of your monomer:>>  "))
# equivlent_RaftAgent = float(input("Please input the initial equivlent of your RaftAgent:>>  "))

DP_Monomerstart = 50 #float(input('\033[1;31mPlease input the proposed start DP of Monomer:>> \033[0m'))
DP_Monomerend = 170 #float(input('\033[1;31mPlease input the proposed end DP of Monomer:>> \033[0m'))

# input the residence time of timesweeps
ResidenceTime = float(input('Please input the involved residence time (seconds):>> '))
ConcentrationChangetime = 600
sleeptime = 5 #float(input('Please input the pump sleep time (seconds):>> '))

MonCon_inRA = 0.28

DPchangetime = 600
Flowrate = V_reactor* 60 / ResidenceTime #converting it to ml/min
print(Flowrate)
cleartime = (V_dead+V_input+V_reactor)*60/Flowrate + 300 #extra 5 minutes to make sure the old solution in the flow reactor is cleared out

Flowrate = V_reactor* 60 / ResidenceTime

FlowRateMonomer = []
ConcentrationMonomer = []
FlowRateRAFTAgent = []
FlowRateSolvent = []
FlowRateInitiator = []
DegP = []

No_steps = int(DPchangetime / sleeptime)
DPDecreaseStep = abs((DP_Monomerstart - DP_Monomerend)) / No_steps
RAFT_con = 3/DP_Monomerstart
FRRaftAgent = RAFT_con*Flowrate/StockCon_RAFTAgent
MonCon_inRA = FRRaftAgent *StockCon_Monomer_inRA/Flowrate
FRMonomer = (Pro_Monomer-MonCon_inRA) * Flowrate / StockCon_Monomer
FRInitiator = (Pro_Monomer*Flowrate / 1000) / StockCon_Initiator
FRSolvent = Flowrate - FRMonomer - FRRaftAgent - FRInitiator

FlowRateMonomer.append(FRMonomer); FlowRateRAFTAgent.append(FRRaftAgent);FlowRateInitiator.append(FRInitiator); FlowRateSolvent.append(FRSolvent)

#ConMonomer = FRMonomer * StockCon_Monomer / Flowrate
ConcentrationMonomer.append(Pro_Monomer)


for i in range(No_steps+1):
    dp = DP_Monomerstart + i * DPDecreaseStep
    FRRaftAgent = (Pro_Monomer*Flowrate / dp) / StockCon_RAFTAgent
    MonCon_inRA = FRRaftAgent *  StockCon_Monomer_inRA / Flowrate
    FRMonomer = (Pro_Monomer - MonCon_inRA) * Flowrate / StockCon_Monomer
    FRInitiator = ( Pro_Monomer*Flowrate/ 1000) / StockCon_Initiator
    FRSolvent = Flowrate - FRMonomer - FRRaftAgent - FRInitiator
    FlowRateMonomer.append(FRMonomer)
    DegP.append(dp)
    FlowRateRAFTAgent.append(FRRaftAgent)
    FlowRateSolvent.append(FRSolvent)
    FlowRateInitiator.append(FRInitiator)
 
n = int(DPchangetime/sleeptime)


# FlowRateMonomer.reverse()
print("1\n")
print(FlowRateMonomer)
print("2\n")
print(FlowRateRAFTAgent)
print("3\n")
print( FlowRateSolvent)
print("3\n")
print( FlowRateInitiator)
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

SleepTime.append(420)

for i in range(n-1):
    Sleep = sleeptime
    SleepTime.append(Sleep)

SleepTime.append(ClearTime)
print(len(SleepTime))
print(f"The sleep time of pumps involved is {SleepTime}")

Totalflowmonomer = 0 
TotalflowRaft= 0 
TotalflowSolvent = 0 
TotalflowInitiator = 0

for i in range(len(SleepTime)):
    Totalflowmonomer += SleepTime[i]/60*FlowRateMonomer[i] 
    TotalflowRaft += SleepTime[i]/60*FlowRateRAFTAgent[i]
    TotalflowSolvent += SleepTime[i]/60*FlowRateSolvent[i]
    TotalflowInitiator += SleepTime[i]/60*FlowRateInitiator[i]  

print(
    f"the volume of Monomer needed in this experiment is {Totalflowmonomer}", 
    f"\nthe volume of Raftagent needed in this experiment is {TotalflowRaft}", 
    f"\nthe volume of Solvent needed in this experiment is {TotalflowSolvent}",
    f"\nthe volume of Initiator needed in this experiment is {TotalflowInitiator}", 
)




# PumpSolvent.start(), PumpMonomer.start(),PumpRaftAgent.start(),PumpInitiator.start()
# sleep(0.5)
# PumpSolvent.changeFlowrate(0),  PumpMonomer.changeFlowrate(0), PumpRaftAgent.changeFlowrate(0), PumpInitiator.changeFlowrate(0)
# sleep(0.5)

# Start experiment
for i in range(n+1):
    # get current time
    today = datetime.now()
    CurrentTime = today.strftime("%H:%M:%S")
    print(f'pumps start at {CurrentTime}')
    PumpMonomer.changeFlowrate(FlowRateMonomer[i])
    PumpSolvent.changeFlowrate(FlowRateSolvent[i])
    PumpRaftAgent.changeFlowrate(FlowRateRAFTAgent[i])
    PumpInitiator.changeFlowrate(FlowRateInitiator[i])
    sleep(SleepTime[i])


PumpSolvent.stop(), PumpMonomer.stop(), PumpRaftAgent.stop(), PumpInitiator.stop()

