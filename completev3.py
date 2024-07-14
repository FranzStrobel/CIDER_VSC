# ==================================================================
# Initialising the packages
# ==================================================================
import py_dss_interface
import numpy as np
import csv
from fmpy import *
from fmpy.util import plot_result  # import the plot function
from fmpy import simulate_fmu
import matplotlib.pyplot as plt
import os
import pandas as pd
from datetime import datetime
# ==================================================================

ini_datetime = datetime.now()

#C:\Program Files\OpenDSS

fmu = 'VSC_FZ.Testes.SAEB1_S2_2.fmu'
#model_description = read_model_description(fmu)
#dump(fmu)  # get information+

# Set simulation options
start = 0.0 # Start time of the simulation (seconds) #fmu start
output_interval = 0.00001 # Step size for the fmu simulation (seconds) 
OpenDSS_step = 0.008 # OpenDSS step
duration = 1.2 # total simulation time

 # Simulation's duration (seconds)



# Results
tempo_f=[]
pca_f=[]
ipca_f=[]
step=[]
Vth=[]
Rr=[]
Lr=[]
comprimento=[]
indices=[]
index0 = 0

powerP_f=[]
powerQ_f=[]

bus_candidate = "BMT163605568"


 # ==================================================================
 # .dss first simulation to generate the initial grid state 
 # ==================================================================

# Solve OpenDSS interface
# Initialize OpenDSS interface
dss = py_dss_interface.DSSDLL(r"C:\Program Files\OpenDSS")
dss_file_path = (r"D:\Documents\Documents\2 - Mestrado\CEFET\0 - Regular\Hibrido\Franz\BHMR27_FZ.dss")
dss.text(f"compile [{dss_file_path}]")
dss.text("Buscoords BHMR27CoordMT.csv")
dss.text("Solve mode=daily")
dss.text("Set Number=12")
dss.solution_solve()
dss.text("export voltages")
dss.text("Solve mode=faultstudy")
dss.solution_solve()
dss.text("export seqz")
# ==================================================================


# ==================================================================
# Complete simulation loop (.dss + FMU)  
# ==================================================================

simulation_time = OpenDSS_step
while simulation_time <= duration+output_interval:  

    # ==================================================================
    # Reading results from .dss file
    # ==================================================================
    
    # Voltage reading function from .dss
    def read_csv_file_V(file_path, columns):
        voltage_data = []
        with open(file_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                extracted_row = {col: row[col] for col in columns}
                voltage_data.append(extracted_row)
        return voltage_data

    file_path = r'D:\Documents\Documents\2 - Mestrado\CEFET\0 - Regular\Hibrido\Franz\alimBHMR27_EXP_VOLTAGES.csv'

    # ==================================================================
    va = read_csv_file_V(file_path, [' Magnitude1'])
    va = [item[' Magnitude1'] for item in va]
    vb = read_csv_file_V(file_path, [' Magnitude2'])
    vb = [item[' Magnitude2'] for item in vb]
    vc = read_csv_file_V(file_path, [' Magnitude3'])
    vc = [item[' Magnitude3'] for item in vc]
    
    teta_a = read_csv_file_V(file_path, [' Angle1'])
    bus_numbers = [item[' Angle1'] for item in teta_a]
    teta_b = read_csv_file_V(file_path, [' Angle2'])
    teta_b = [item[' Angle2'] for item in teta_b]
    teta_c = read_csv_file_V(file_path, [' Angle3'])
    teta_c = [item[' Angle3'] for item in teta_c]

    v_base = read_csv_file_V(file_path, [' BasekV'])
    v_base = [item[' BasekV'] for item in v_base]
    Bus_name = read_csv_file_V(file_path, ['Bus'])
    Bus_name = [item['Bus'] for item in Bus_name]

    try:
        # Attempt to remove the file
        os.remove(file_path)
    except OSError as e:
        # Handle the case where the file couldn't be removed
        print(f"Error: {file_path} - {e}")
    # ==================================================================


    # Impedance reading fuction from .dss
    def read_csv_file_Z(file_path, columns):
        impedance_data = []
        with open(file_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                extracted_row = {col: row[col] for col in columns}
                impedance_data.append(extracted_row)
        return impedance_data

    file_path = r'D:\Documents\Documents\2 - Mestrado\CEFET\0 - Regular\Hibrido\Franz\alimBHMR27_EXP_SEQZ.csv'
    #impedance_data = read_csv_file_Z(file_path)

    # ==================================================================
    R1 = read_csv_file_Z(file_path, [' R1'])
    R1 = [item[' R1'] for item in R1]
    X1 = read_csv_file_Z(file_path, [' X1'])
    X1  = [item[' X1'] for item in X1 ]  

    try:
        # Attempt to remove the file
        os.remove(file_path)
    except OSError as e:
        # Handle the case where the file couldn't be removed
        print(f"Error: {file_path} - {e}")
    # ==================================================================

    
    # Selecting the bus candidate´s data to load the variables: V, corresponding_R1, corresponding_L1

    bus_position = Bus_name.index(bus_candidate)
    # Access the corresponding values from R1 and X1
    corresponding_R1 = R1[bus_position]
    corresponding_X1 = X1[bus_position]    
    corresponding_L1 = float(corresponding_X1) / (2 * 3.14159 * 60)
    V = float(va[bus_position])*1.414214
    
    Rr.append(corresponding_R1)
    Lr.append(corresponding_L1)
    Vth.append(V)
        
       
    # ==================================================================
    # FMU simulation
    # ==================================================================

    # Inputs for FMU simulation
    parameters = {'Vth': V, 'Rr': corresponding_R1, 'Lr': corresponding_L1} #Loading FMU variables (left) from .dss variables (right) , 'rampP':750000,'rampQ':0
   
    # Specify the quantities to export in the result
    outputs = ['time', 'powerCalc_pcc.Q', 'powerCalc_pcc.P', 'ipca.i', 'pca.v']   # Values requested from FMU 
        
    # Simulate the FMU with modified parameters
    print("simulation_time =", simulation_time)
       
    result = simulate_fmu(fmu, start_time=start,
                          stop_time=simulation_time,
                          relative_tolerance = 0.000001,
                          output_interval=output_interval,
                          start_values=parameters,
                          output=outputs)

       
    # Retrieve numerical data from the simulation result
    time = result['time']  # Time vector
    powerQ = result['powerCalc_pcc.Q']  
    powerP = result['powerCalc_pcc.P']
    ipca = result['ipca.i']
    pca = result['pca.v']
        
    index1 = len(time)- index0
    index0 = len(time)
    round(OpenDSS_step/(time[2]-time[1])) #calculation of elements that will be stored from the last fmu simulation
   
    # stored values for final results
    comprimento.append(len(time))
    indices.append(index1)
    step.append((time[2]-time[1]))
    powerP_f.extend(powerP[-index1:]) 
    powerQ_f.extend(powerQ[-index1:])
    ipca_f.extend(ipca[-index1:])
    pca_f.extend(pca[-index1:])
    tempo_f.extend(time[-index1:])
      
   
    # ==================================================================
    # .dss simulation considering the updated outputs from FMU simulation
    # ==================================================================

    # Updating OpenDSS File 
    # Initialize OpenDSS interface
    dss = py_dss_interface.DSSDLL(r"C:\Program Files\OpenDSS")
    dss_file_path = (r"D:\Documents\Documents\2 - Mestrado\CEFET\0 - Regular\Hibrido\Franz\BHMR27_FZ.dss")
    dss.text(f"compile [{dss_file_path}]")
    #dss.text(f"Edit Load.BMT167029211a Bus1=BMT167029211.1    Phases=1 Conn=Wye  Model=1 kV=2.4  kW={100*current[-1]}    kvar={100*current[-1]}")
    #dss.text(f"Edit Load.BMT167029211b Bus1=BMT167029211.2    Phases=1 Conn=Wye  Model=1 kV=2.4  kW={100*current[-1]}    kvar={100*current[-1]}")
    #dss.text(f"Edit Load.BMT167029211c Bus1=BMT167029211.3    Phases=1 Conn=Wye  Model=1 kV=2.4  kW={100*current[-1]}  kvar={100*current[-1]}")
    print("P =", powerP_f[-1])   
    dss.text(f"Edit load.BESS1 bus1=BMT163605568.1.2.3.0,Phases=3,kv=13.800000000,kW={-powerP_f[-1]/1000},pf=1,Vminpu=0.93,Vmaxpu=1.5,model=1,status=variable")
    #dss.text(f"Edit Storage.BESS1 bus1=BMT9697784B,kVA=750,kW={powerP[-1]/1000},kvar={0}")
    #dss.text(f"Edit Storage.BESS2 bus1=BMT181448789,kVA=450,kW={current[-1]/1000},kvar={current[-1]/1000}")

    dss.text("Buscoords BHMR27CoordMT.csv")
    dss.text("Solve mode=daily")
    dss.text("Set Number=12")
    dss.solution_solve()
    dss.text("export voltages")
    dss.text("Solve mode=faultstudy")
    dss.solution_solve()
    dss.text("export seqz")
    # ==================================================================

   
    simulation_time += OpenDSS_step
  
# ==================================================================
# Complete simulation loop close (.dss + FMU)  
# ==================================================================


# Save results
#df_general = pd.DataFrame({'startfmu': start,'intervalo fmu':output_interval, 'passo OpenDSS':OpenDSS_step, 'duração total':duration})
df_resultados_OpenDSS = pd.DataFrame({'Vth': Vth, 'Rr': Rr, 'Lr':Lr})
df_resultados_fmu = pd.DataFrame({'Tempo':tempo_f, 'ipca': ipca_f, 'vpca':pca_f,'PowerP':powerP_f, 'PowerQ':powerQ_f})

df_resultados_OpenDSS.to_csv('valida1_OpenDSS_0_004.csv', index=False)
df_resultados_fmu.to_csv('valida1_fmu_0_004.csv', index=False)
#df_general.to_csv('simulacao1_geral.csv', index=False)

# Plot results
fin_datetime = datetime.now()

print("Início =",ini_datetime)
print("Fim =",fin_datetime)
print("Passos =",step)
print("Comprimentos =",comprimento)
print("indices =",indices)
print("Tempofinal =", len(tempo_f))
print("Pfinal =",len(powerP_f))
print("Ipcafinal =",len(ipca_f))
print("Vpcafinal =",len(pca_f))
print("Vth =",Vth)
print("Rr =",Rr)
print("Lr =",Lr)

plt.plot(tempo_f, ipca_f)
plt.xlabel('Time')
plt.ylabel('Current')
plt.title('ipca')
plt.grid(True)
plt.show()

plt.plot(tempo_f, pca_f)
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.title('pca')
plt.grid(True)
plt.show()