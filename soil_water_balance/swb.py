"""
A rainfall - runoff-  infiltration model. Like SWB.
"""
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
from math import pow

df = pd.read_excel("panda2.xlsx")
elev = np.array(df['elev'])
hsg = np.array(df['hsg'])
ksat = []
hsg_reclass = []
for i in hsg:
    if i == 0:
        hsg_reclass.append('A')
        ksat.append(1000)
    elif i == 1:
        hsg_reclass.append('B')
        ksat.append(500)
    elif i == 2:
        hsg_reclass.append('C')
        ksat.append(1500)
    else:
        hsg_reclass.append('D')
        ksat.append(1200)


lc = np.array(df['lc'])
elevation = np.reshape(elev,(30,30))

class Cell (object):
    def __init__(self, id, downstream_cell, KSAT, INTC, LC, HSG):
        self.id = id # id indexing starts from top left corner end ends on bottom right corner of the grid
        self.KSAT = KSAT # Saturated vertical hydraulic conductivity of the topsoil layer, m/day
        self.INTC = INTC # Interception capacity, m
        self.LC = LC # land cover type
        self.HSG = HSG # Hydrological soil group
        self.five_day_pp = [] # Latest 5-day precipitation list
        self.curve_number = curve_number(self.LC, self.HSG) # The curve_number function which requires land cover type index and hydrological soil group is called
        self.downstream_cell = int(downstream_cell)
    def set_initial(self, pp, evap, int_stor, s_stor):
        self.PP = float(pp) #Precipitation, m/day
        self.init_INT_STOR = int_stor - evap * 0.85 if int_stor - evap * 0.85 > 0 else 0. # initial intercepted water storage, m
        self.init_S_STOR = s_stor - evap * 0.15 if s_stor - evap * 0.15 > 0 else 0.# initial surface water storage, m
        self.runoff_in = 0. # Inflow runoff prior the time-step begin, m
        self.runoff_out = 0. # Outflow runof prior the time-step begin, m
        self.ACT_EVAP = evap if evap < int_stor + s_stor else int_stor + s_stor
        # 5-day precipitation ammount is required to apply curve number index correction
        self.five_day_pp.append(self.PP)
        if self.five_day_pp > 5: # Limiting list to 5 elements
            self.five_day_pp = self.five_day_pp[1:]

    def calculate (self):

        self.res_INT_STOR = self.init_INT_STOR
        self.res_S_STOR = self.init_S_STOR
        self.res_INFIL = 0. # Infiltration from the land surface to soil, m
        self.delta_INT_STOR = 0.
        self.delta_S_STOR = 0.

        S = (1000./self.curve_number) - 10
        Ia = 0.05 * S
        Q = (((self.runoff_in + self.PP)/25) * ((self.runoff_in + self.PP)/25)) / ((self.runoff_in + self.PP)/25 - Ia + S)
        self.runoff_out = Q*25

        if self.runoff_in + self.PP - self.runoff_out < self.INTC:
            self.delta_INT_STOR = self.runoff_in + self.PP - self.runoff_out
            self.res_INT_STOR = self.init_INT_STOR + self.delta_INT_STOR if self.init_INT_STOR + self.delta_INT_STOR > 0 else 0.
        else:
            self.delta_INT_STOR = self.INTC - self.init_INT_STOR
            self.res_INT_STOR = self.INTC

        if self.runoff_in + self.PP + self.init_S_STOR - self.runoff_out - self.delta_INT_STOR < self.KSAT:
            self.res_INFIL = self.runoff_in + self.PP + self.init_S_STOR - self.runoff_out - self.delta_INT_STOR if self.runoff_in + self.PP + self.init_S_STOR - self.runoff_out - self.delta_INT_STOR > 0 else 0.
        else:
            self.res_INFIL = self.KSAT

        self.delta_S_STOR = self.runoff_in + self.PP - self.runoff_out - self.delta_INT_STOR - self.res_INFIL
        self.res_S_STOR = self.init_S_STOR + self.delta_S_STOR if self.init_S_STOR + self.delta_S_STOR > 0 else 0.

        self.runoff_in = 0.

    def give_runoff(self, model_object):
        model_object.cells[self.downstream_cell].take_runoff(self.runoff_out)
        self.runoff_out = 0.

    def take_runoff(self, runoff_surplus):
        self.runoff_in += runoff_surplus



class Model(object):

    def __init__(self, ncol, nrow, downstream_list, KSAT, LC, HSG):
        self.ids = range(ncol * nrow)
        self.ncol = ncol
        self.nrow = nrow
        self.cells = []
        for cell_id in self.ids:
            self.cells.append(Cell(id = cell_id, downstream_cell = downstream_list[cell_id], KSAT = KSAT, INTC = 5., LC = LC[cell_id], HSG = HSG[cell_id]))



def model_run (model_object, pp_list, et_list):
    # checking for valid input time-series
    if len(pp_list) != len(et_list) or len(pp_list) == 0:
        print 'invalid data'
        return
    runtime = range(len(pp_list)) # List of time step indexes
    infiltration ={} # Dictionary that will collect computed infiltration data
    runoff = {}
    act_evap ={}
    INT_STOR = []
    S_STOR = []
    for step in runtime:
        for index, cell_x in enumerate(model_object.cells):
            cell_x.set_initial(pp_list[step], et_list[step], int_stor = INT_STOR[index] if len(INT_STOR) > 0 else 0. , s_stor = S_STOR[index] if len(S_STOR) > 0 else 0.)

        runoff2distribute = True
        iteration = 0

        while runoff2distribute:
            infiltration_ts = []
            act_evap_ts = []
            INT_STOR = []
            S_STOR = []
            OUTFLOW = 0.
            sum_runoff_in_prior = 0.
            sum_runoff_in_posterior = 0.
##########################################################################################################
            for cell in model_object.cells:
                sum_runoff_in_prior += cell.runoff_in
            for cell in model_object.cells:
                cell.calculate()
                INT_STOR.append(cell.res_INT_STOR)
                S_STOR.append(cell.res_S_STOR)
                act_evap_ts.append(cell.ACT_EVAP)
                infiltration_ts.append(cell.res_INFIL)
                if cell.id == cell.downstream_cell:
                    OUTFLOW += cell.runoff_out
                    cell.runoff_out = 0.

                if cell.runoff_out != 0 : # avoiding zero division
                    cell.give_runoff(model_object)

            for cell in model_object.cells:
                sum_runoff_in_posterior += cell.runoff_in
#########################################################################################################
            if abs(sum_runoff_in_prior - sum_runoff_in_posterior) <=0.1:
                runoff2distribute = False
            print step, iteration, sum_runoff_in_prior - sum_runoff_in_posterior
            iteration += 1
        # Calculated infiltration dictionary in format {time-step n:[cell 1 id...cell i id]}
        runoff[step] = OUTFLOW
        infiltration[step] = np.reshape(infiltration_ts,(model_object.ncol,model_object.nrow))
        act_evap[step] = np.reshape(act_evap_ts,(model_object.ncol,model_object.nrow))

    return (runoff, infiltration, act_evap)

def curve_number(lc, hsg):
    """
    Reads the excel cunve numbers catalog and returnes a CN value corresponding to given landcover type and hydrological soil group. Used to compute generated runoff at the model cells - Cell.calculate()
    """
    df = pd.read_excel("cn.xlsx") # excel file containing table with curve numbers is read and saved into a pandas dataframe object
    return int(df[hsg][df['id'] == lc]) # reading curve number corresponding to given lc and hsg from the dataframe


def flow_direction (elevation):
    """
    Function for computing flow direction raster from the given elevation raster. Returnes a list of downstream cells for every cell index.
    """
    downstream_list = []
    for i, row in enumerate(elevation):
        for j, cell in enumerate(row):

            minimal = cell
            index = len(row) * i + j

            hohe = elevation[i-1][j-1] if i != 0 and j != 0 else None
            if hohe <= minimal and hohe != None:
                index = len(row) * (i - 1) + (j - 1)
                minimal = hohe

            hohe = elevation[i-1][j] if i != 0 else None
            if hohe <= minimal and hohe != None:
                index = len(row) * (i - 1) + (j)
                minimal = hohe

            hohe = elevation[i-1][j+1] if i != 0 and j != len(row) - 1 else None
            if hohe <= minimal and hohe != None:
                index = len(row) * (i - 1) + (j + 1)
                minimal = hohe

            hohe = elevation[i][j-1] if j != 0 else None
            if hohe <= minimal and hohe != None:
                index = len(row) * (i) + (j - 1)
                minimal = hohe

            hohe = elevation[i][j+1] if j != len(row) - 1 else None
            if hohe <= minimal and hohe != None:
                index = len(row) * (i) + (j + 1)
                minimal = hohe

            hohe = elevation[i+1][j-1] if i != len(elevation) - 1 and j != 0 else None
            if hohe <= minimal and hohe != None:
                index = len(row) * (i + 1) + (j - 1)
                minimal = hohe

            hohe = elevation[i+1][j] if i != len(elevation) - 1 else None
            if hohe <= minimal and hohe != None:
                index = len(row) * (i + 1) + (j)
                minimal = hohe

            hohe = elevation[i+1][j+1] if i != len(elevation) - 1 and j != len(row) - 1 else None
            if hohe <= minimal and hohe != None:
                index = len(row) * (i + 1) + (j + 1)
                minimal = hohe

            downstream_list.append(index)
    print 'FLOW DIRECTION RASTER IS CREATED'
    return downstream_list



######################### Test #################################################
#elevation = [[100.0,101.1,101.2],[99.9,102.0,103.0],[99.0,101.9,105.0]]
pp_list = []
et_list = []
with open('meteo_m.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        pp_list.append(float(row['pp']))
        et_list.append(float(row['et']))

fd = flow_direction(elevation)
model = Model(len(elevation),len(elevation[0]),fd,ksat,lc,hsg_reclass)
runoff, infiltration, act_evap = model_run(model, pp_list, et_list)
total_infiltration = []
total_runoff = []
total_act_evap = []

for step in infiltration:
     total_act_evap.append(np.mean(act_evap[step]))
for step in infiltration:
     total_infiltration.append(np.mean(infiltration[step]))
for step in runoff:
     total_runoff.append(runoff[step]/(len(elevation)*len(elevation[0])))
d = {'precip':pp_list, 'evap':total_act_evap, 'runoff':total_runoff, 'infiltr': total_infiltration}
summary = pd.DataFrame(d)


ksat = np.array(ksat)
hsg_reclass = np.array(hsg_reclass)
ksat = np.reshape(ksat,(30,30))
hsg_reclass = np.reshape(hsg_reclass,(30,30))
lc = np.reshape(lc,(30,30))
print 'ELEVATION MAP'
plt.imshow(elevation, interpolation='nearest')
plt.colorbar()
plt.show()
print 'Kf MAP'
plt.imshow(ksat, interpolation='nearest')
plt.colorbar()
plt.show()
print 'LANDCOVER MAP'
plt.imshow(lc, interpolation='nearest')
plt.colorbar()
plt.show()
print 'SOIL MAP'
plt.imshow(np.reshape(hsg,(30,30)), interpolation='nearest')
plt.colorbar()
plt.show()
print 'INFILTRATION DAY 290'
plt.imshow(infiltration[290], interpolation='nearest')
plt.colorbar()
plt.show()
print 'INFILTRATION DAY 380'
plt.imshow(infiltration[380], interpolation='nearest')
plt.colorbar()
plt.show()
print summary.describe()
print 'balance error, %: ' + str((sum(d['precip'])-sum(d['infiltr'])-sum(d['runoff'])-sum(d['evap']))/sum(d['precip']) * 100)
