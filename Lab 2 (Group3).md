Lynn Li (ll643) and Allison Tran (ant42)
Group 3
Lab 3

#Question 1
```python
from aguaclara.core.units import unit_registry as u
import aguaclara.research.environmental_processes_analysis as epa
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
data_file_path="https://raw.githubusercontent.com/ll643/Lab-2/master/Lab%202%20Data%20(Group3).tsv"
df = pd.read_csv(data_file_path,delimiter='\t')
array = np.array(df)
print(df)
V = 3.785*(u.L)
Q =(0.00537)*(u.L/u.s)
theta = (V/Q)
time = (array[:,0]*(u.day)-0.60090697*(u.day))#0.60090697 is time 0
time = time.to(u.s)
hydraulic_residence_time = (time/theta).to(u.dimensionless)
pH = array[:,1]
plt.figure()
plt.plot(hydraulic_residence_time, pH,'-')
plt.xlabel('Hydraulic Residence Time (dimensionless)')
plt.ylabel('pH of Lake')
plt.title('pH of Lake vs. Hydraulic Residence Time')
plt.show()
```
![figure1](https://github.com/ll643/Lab-2/blob/master/Lab%203%20(Group3)%20Pic1.png)
Figure 1

#Question 2
##Givens:
$$Volume_{lake}=4quarts=3.785L$$
$$pH_{lake(init)}=7$$
$$Q=161ml/30s=0.161L/30s=0.00537L/s$$
##Work:
$${ANC\; }={\; }\left[{HCO}_{{3}}^{{-}} \right]+{\; 2}\left[{CO}_{{3}}^{{-2}} \right]+\left[{OH}^{{-}} \right]{\; -}\left[{H}^{+} \right]$$
$${ANC_{0} }={\; }\left[{HCO}_{{3}}^{{-}} \right]-[{H}^{+}]$$
$$[{HCO}_{{3}}^{{-}}]=mol_{HCO_3^{-}}/Volume_{lake}=(0.623g/\frac{84.007g}{mol})/3.785L=1.96*10^{-3}mol/L$$
$$pOH=14-pH=14-7=7$$
$$[OH^{-1}]=10^{-pH}=10^{-7}$$
$$\Rightarrow-[H^+]=-10^{-14}/10^{-7}=-1*10^{-7}mol/L$$
$$\Rightarrow ANC_{0}=(1.96*10^{-3}mol/L)-1*10^{-7}mol/L=1.96*10^{-3}mol/L$$

```python
data_file_path="https://raw.githubusercontent.com/ll643/Lab-2/master/Lab%202%20Data%20(Group3).tsv"
df = pd.read_csv(data_file_path,delimiter='\t')
array = np.array(df)
V = 3.785*(u.L)
Q =(0.00537)*(u.L/u.s)
theta = (V/Q)
time = (array[:,0]*(u.day)-0.60090697*(u.day))#0.60090697 is time 0
time = time.to(u.s)
hydraulic_residence_time = (time/theta).to(u.dimensionless)
ANC_0 = 0.00196
ANC_in = -0.001
ANC_con = (ANC_in*(1-np.exp(-hydraulic_residence_time)))+(ANC_0*np.exp(-hydraulic_residence_time))
plt.figure()
plt.plot(hydraulic_residence_time, ANC_con, label='Conservative ANC')
plt.xlabel('Hydraulic Residence Time (dimensionless)')
plt.ylabel('Conservative ANC (eq/L)')
plt.title('Conservative ANC vs. Hydraulic Residence Time')
plt.legend()
plt.show()
```
#Question 3
```python
K1 = 10**-6.3
K2 = 10**-10.3
Kw = 10**(-14)
P_CO2 = 10**(-3.5)
Kh = 10**(-1.5)
V = 3.785*(u.L)
Q = (0.00537)*(u.L/u.s)
theta = V/Q
data_file_path="https://raw.githubusercontent.com/ll643/Lab-2/master/Lab%202%20Data%20(Group3).tsv"
df = pd.read_csv(data_file_path,delimiter='\t')
array = np.array(df)
time = (array[:,0]*(u.day)-0.60090697*(u.day))#0.60090697 is time 0
time = time.to(u.s)
hydraulic_residence_time = (time/theta).to(u.dimensionless)
Conc_NaCO3 = 0.00196
pH= array[:,1]
Conc_H=10**(-1*pH)
a0 = 1/(1+(K1/Conc_H)+(K1*K2/(Conc_H**2)))
a1 = 1/((Conc_H/K1)+1+(K2/Conc_H))
a2 = 1/(((Conc_H**2)/(K1*K2))+(Conc_H/K2)+1)
Closed_ANC = (Conc_NaCO3*(a1+(2*a2)))+(Kw/Conc_H)-Conc_H
plt.plot(hydraulic_residence_time, Closed_ANC, label = 'Closed ANC')
plt.xlabel('Hydraulic Residence Time (dimensionless)')
plt.ylabel('Closed ANC (eq/L)')
plt.title('Closed ANC vs. Hydraulic Residence Time')
plt.legend()
plt.show()
```
#Question 4
```python
K1 = 10**-6.3
K2 = 10**-10.3
Kw = 10**(-14)
P_CO2 = 10**(-3.5)
Kh = 10**(-1.5)
V = 3.785*(u.L)
Q = (0.00537)*(u.L/u.s)
theta = V/Q
data_file_path="https://raw.githubusercontent.com/ll643/Lab-2/master/Lab%202%20Data%20(Group3).tsv"
df = pd.read_csv(data_file_path,delimiter='\t')
array = np.array(df)
time = (array[:,0]*(u.day)-0.60090697*(u.day))#0.60090697 is time 0
time = time.to(u.s)
hydraulic_residence_time = (time/theta).to(u.dimensionless)
Conc_NaCO3=0.00196
pH= array[:,1]
Conc_H=10**(-1*pH)
a0 = 1/(1+(K1/Conc_H)+(K1*K2/(Conc_H**2)))
a1 = 1/((Conc_H/K1)+1+(K2/Conc_H))
a2 = 1/(((Conc_H**2)/(K1*K2))+(Conc_H/K2)+1)
Open_ANC = (((P_CO2*Kh)/a0)*(a1+(2*a2)))+(Kw/Conc_H)-Conc_H
plt.plot(hydraulic_residence_time, Open_ANC, label = 'Open ANC')
plt.xlabel('Hydraulic Residence Time (dimensionless)')
plt.ylabel('Open ANC (eq/L)')
plt.title('Open ANC vs. Hydraulic Residence Time')
plt.legend()
plt.show()
```
#Final ANC Graph
```python
K1 = 10**-6.3
K2 = 10**-10.3
Kw = 10**(-14)
P_CO2 = 10**(-3.5)
Kh = 10**(-1.5)
V = 3.785*(u.L)
Q = (0.00537)*(u.L/u.s)
theta = V/Q
data_file_path="https://raw.githubusercontent.com/ll643/Lab-2/master/Lab%202%20Data%20(Group3).tsv"
df = pd.read_csv(data_file_path,delimiter='\t')
array = np.array(df)
time = (array[:,0]*(u.day)-0.60090697*(u.day))#0.60090697 is time 0
time = time.to(u.s)
hydraulic_residence_time = (time/theta).to(u.dimensionless)
Conc_NaCO3 = 0.00196
pH= array[:,1]
Conc_H=10**(-1*pH)
a0 = 1/(1+(K1/Conc_H)+(K1*K2/(Conc_H**2)))
a1 = 1/((Conc_H/K1)+1+(K2/Conc_H))
a2 = 1/(((Conc_H**2)/(K1*K2))+(Conc_H/K2)+1)
ANC_0 = 0.00196
ANC_in = -0.001
ANC_con = (ANC_in*(1-np.exp(-hydraulic_residence_time)))+(ANC_0*np.exp(-hydraulic_residence_time))
Closed_ANC = (Conc_NaCO3*(a1+(2*a2)))+(Kw/Conc_H)-Conc_H
Open_ANC = (((P_CO2*Kh)/a0)*(a1+(2*a2)))+(Kw/Conc_H)-Conc_H
plt.plot(hydraulic_residence_time, ANC_con, label='Conservative ANC')
plt.plot(hydraulic_residence_time, Closed_ANC, label = 'Closed ANC')
plt.plot(hydraulic_residence_time, Open_ANC, label = 'Open ANC')
plt.xlabel('Hydraulic Residence Time (dimensionless)')
plt.ylabel('ANC (eq/L)')
plt.title('ANC vs. Hydraulic Residence Time')
plt.legend()
plt.show()
```
![figure2](https://github.com/ll643/Lab-2/blob/master/Lab%203%20(Group3)%20Pic2.png)

#Question 5
##pH
```python
from aguaclara.core.units import unit_registry as u
import aguaclara.research.environmental_processes_analysis as epa
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
data_file_path="https://raw.githubusercontent.com/ll643/Lab-2/master/Lab%202%20Round%202%20Data%20(group3).tsv"
df = pd.read_csv(data_file_path,delimiter='\t')
array = np.array(df)
time = (array[:,0]*(u.day)-0.628707843*(u.day))#0.628707843 is time 0
time = time.to(u.s)
hydraulic_residence_time = (time/theta).to(u.dimensionless)
V = 3.785*(u.L)
Q =(0.00537)*(u.L/u.s)
theta = (V/Q)
pH = array[:,1]
plt.figure()
plt.plot(hydraulic_residence_time, pH,'-')
plt.xlabel('Hydraulic Residence Time (dimensionless)')
plt.ylabel('pH of Lake')
plt.title('pH of Lake vs. Hydraulic Residence Time')
plt.show()
```
![figure3](https://github.com/ll643/Lab-2/blob/master/Lab%203%20(Group3)%20Pic3.png)

##ANC
```python
K1 = 10**-6.3
K2 = 10**-10.3
Kw = 10**(-14)
P_CO2 = 10**(-3.5)
Kh = 10**(-1.5)
V = 3.785*(u.L)
Q = (0.00537)*(u.L/u.s)
theta = V/Q
data_file_path="https://raw.githubusercontent.com/ll643/Lab-2/master/Lab%202%20Round%202%20Data%20(group3).tsv"
df = pd.read_csv(data_file_path,delimiter='\t')
array = np.array(df)
time = (array[:,0]*(u.day)-0.628707843*(u.day))#0.628707843 is time 0
time = time.to(u.s)
hydraulic_residence_time = (time/theta).to(u.dimensionless)
Conc_NaCO3 = 0.00196
pH= array[:,1]
Conc_H=10**(-1*pH)
a0 = 1/(1+(K1/Conc_H)+(K1*K2/(Conc_H**2)))
a1 = 1/((Conc_H/K1)+1+(K2/Conc_H))
a2 = 1/(((Conc_H**2)/(K1*K2))+(Conc_H/K2)+1)
ANC_0 = 0.00196
ANC_in = -0.001
ANC_con = (ANC_in*(1-np.exp(-hydraulic_residence_time)))+(ANC_0*np.exp(-hydraulic_residence_time))
Closed_ANC = (Conc_NaCO3*(a1+(2*a2)))+(Kw/Conc_H)-Conc_H
Open_ANC = (((P_CO2*Kh)/a0)*(a1+(2*a2)))+(Kw/Conc_H)-Conc_H
plt.plot(hydraulic_residence_time, ANC_con, label='Conservative ANC')
plt.plot(hydraulic_residence_time, Closed_ANC, label = 'Closed ANC')
plt.plot(hydraulic_residence_time, Open_ANC, label = 'Open ANC')
plt.xlabel('Hydraulic Residence Time (dimensionless)')
plt.ylabel('ANC (eq/L)')
plt.title('ANC vs. Hydraulic Residence Time')
plt.legend()
plt.show()
```
![figure4](https://github.com/ll643/Lab-2/blob/master/Lab%203%20(Group3)%20Pic4.png)

What we learned in the second trial is that by doubling the initial ANC amount (from 0.623g of NaCO3 to 1.246g of NaCO3), the lake was better at buffering the acid rain. Based off of the figures produced, we can see that the final pH of the lake after 20 mins with the extra ANC (6) is much higher than the final pH of the lake after 20 mins with the normal dosage of ANC (4).

#Question 1
##What do you think would happen if enough NaHCO3 were added to the lake to maintain an ANC greater than 50μeq/L for 3 residence times with the stirrer turned off? How much NaHCO3 would need to be added?

If enough NaHCO3 was added to the lake to maintain an ANC greater than 50μeq/L for 3 residence times with the stirrer turned off, the lake would be able to withstand the the acid rain for a longer time before it reached its buffering threshold and surrended to acidification. However, because the lake is not being stirred, certain areas of the lake will not become as effective at buffering the incoming acid rain. Visually, this change would cause the lake to stay blue longer before it started to become a clear yellow/green color.

```python
from aguaclara.core import utility as ut
ANC_out=0.00005 #u.eq/u.L
ANC_in=-0.001 #u.eq/u.L
hydraulic_residence_time = 3
volume_lake = 4
molar_mass_NaHCO3 = 84.007
ANC_0=(ANC_out-(ANC_in*(1-np.exp(-hydraulic_residence_time))))/(np.exp(-hydraulic_residence_time))
print(ANC_0)
Conc_H_Plus=1*(10**(-7))
Conc_HCO3 = ANC_0+Conc_H_Plus
Mass_NaHCO3=(Conc_HCO3*volume_lake)*(molar_mass_NaHCO3)
print('We would need', ut.round_sf(Mass_NaHCO3,4),'g of NaHCO3.')
```

#Question 2
##What are some of the complicating factors you might find in attempting to remediate a lake using CaCO3? Below is a list of issues to consider.
*extent of mixing
*solubility of CaCO3 (find the solubility and compare with NaHCO3)
*density of CaCO3 slurry (find the density of CaCO3)

One of the complicating factors that we might find in attempting to remediate a lake using CaCO3 is the extremely low solubility of CaCO3 in pure water (15mg/L) which would make it difficult to dissolve therefore making it not as effective of a remediating agent. In addition its solubility is variant with the temperature of the water (it is inversly proportional to temperature). In contrast, NaHCO3 is extremely water-soluble no matter the temperature of the water. However, it is important to note that due to this descrepency, CaCO3 could potentially be more soluble than NaHCO3 if we are using cold water. Another complicating factor we must take into consideration is that fact that CaCO3 as a slurry is denser than NaHCO3 (2.74 g/cm3 vs 2.2g/cm3) and therefore would be more difficult to mix thoroughly throughout the lake. These factors should be kept in mind when deciding on a remediation process with CaCO3, as calcium carbonate is more selectively soluble and denser than NaHCO3.
