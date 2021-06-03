# Reaction Mechanism Analysis
# William Pearson
# With the help of Dr. Ramanathan Srinivasan and Dr. Fathima Fasmin 

###Packages###
import csv  # Read and Write Files
import math # Use Math Functions (sqrt, e)
import numpy as np # Use Matrices
from scipy.optimize import minimize # Optimise Function
import matplotlib.pyplot as plt # Plot Graphs

###Functions###

def OpenPotentialCurrentDataFile(FileName):
    ##Opens a file using the 'openwith' function,
    ##  Turns the text into a float and adds the
    ##  Potenial and Current values to individual lists
    
    xValues=[]  # Create empty lists ro return
    yValues=[]

    with open(FileName) as DataSet: # Opens the file
        DataSetList = csv.reader(DataSet, delimiter='\t')   #Reads and seperates the file by 'tab'
        for Row in DataSetList:
            xValues.append(float(Row[0]))   # Adds the first number in row(i) to the xValues lists (Potential)
            yValues.append(float(Row[1]))   # Adds the second number in row(i) to the yValues lists (Current)
    return [xValues, yValues]


def OpenEISDataFiles(FileNames):
    ## Opens up each EIS file using the 'openwith' function,
    ##  Turns the text into a float and adds the
    ##  Frequency, ZReal, and ZImaginary values to individual lists

    AllEISData = [] # Create empty list to return

    for IndividualEISDataFile in FileNames:
        

        xValues=[]  # Create empty lists
        yValues=[]
        zValues=[]
        IndividualEISData=[]
        
        with open(IndividualEISDataFile) as DataSet:    # Opens the file
            DataSetList = csv.reader(DataSet, delimiter='\t')   # Reads and seperates the file by 'tab'
            for Row in DataSetList:
                xValues.append(float(Row[0]))   # Adds the first number in row(i) to the xValues lists (Frequency)
                yValues.append(float(Row[1]))   # Adds the second number in row(i) to the yValues lists (Zreal)
                zValues.append(float(Row[2]))   # Adds the second number in row(i) to the yValues lists (Zimaginary)

            IndividualEISData = [xValues, yValues, zValues] # Puts all the values from above together
            AllEISData.append(IndividualEISData)    # Adds the individual data to one list
            
    return AllEISData # Returns list of all EIS Data


def ExtractDataFromEISDataSet(ExperimentalEISData):
    ## 
    ExtractedData = []  # Create empty list to return 
    
    for DataSet in ExperimentalEISData:
        ImpedanceDataSet = []   # Create empty list
        for Row in range(len(DataSet[0])): # Loop through every Row in each EIS Data Set.  [0] to loop through the length of the set
            newImpedance = complex(DataSet[1][Row], DataSet[2][Row])    # Combine the Zr and Zim into a complex number
            ImpedanceDataSet.append(newImpedance)   # Add the new complex number to a list for each EIS Data Set
            
        ExtractedData.append(ImpedanceDataSet)  # Add the new whole list to the return List

    ExtractedData = np.array(ExtractedData, dtype=object) # Turn into a numpy array to make matrix calculations easier later
    
    return ExtractedData    # Returns list of extracted Zr and Zim as a complex number


def CalculateImpedanceWeightFactor(PotentialSet, FrequencySet, IWF):
    ##  Each potential has its own EIS[freq] data
    ##  Calculate the weight factor depending on frequency for every frequency at every potential
    
    ImpedanceWeightFactor = []  # Create empty list to return

    for i in range(len(PotentialSet)):  # Go through each freqeuncy data
        NewWeightFactorList = []     # Create empty list
        
        for j in FrequencySet[i]:   # Go through every frequency entry
            IndividualWeightFactor = 1 / (j**IWF[i])    # Calculate the weight factor
            NewWeightFactorList.append(IndividualWeightFactor)  # Add every value to a temporary list 
            
        ImpedanceWeightFactor.append(NewWeightFactorList)   # Add the temporary list to the return list

    ImpedanceWeightFactor = np.array(ImpedanceWeightFactor, dtype=object)

    return ImpedanceWeightFactor    # Return List of all Weight Factors


def CalculateError(Parameters):
    ## Error is caluclated as (Zimulated Data - Actual Data) multiplied by weighing factors
    
    ErrorValues = np.array([]) # Create empty list
    
    for i in range(len(PotentialSet)):
        ZSimulated = SingleVoltage_EISSimulation(   # Returns a numpy array of simulated EIS data of exectly the same size as ZDataSet
                                PotentialSet[i], 
                                Parameters,RSolSet[i],QSet[i],AlphaSet[i],FrequencySet[i],FirstSpeciesImpVector[i],SecondSpeciesImpVector[i])


        NewErrorValues = (ZSimulated - ZDataSet[i])*ImpedanceWeightFactor[i]*DataWeighingFactor[i] # Calculates the Error Vector
        ErrorValues = np.append(ErrorValues, np.array(NewErrorValues)) # Adds the error vector to all the list of error vectors

    NormalisedErrorValue = np.linalg.norm(ErrorValues)  # Calculates the magnitude of the Vector

    return NormalisedErrorValue     # This is the value we want to be minimised


def NonLinearConstraint(Parameters):
    ## Calculates the Non-linear constraints
    ##  In this case the absolute difference between the
    ##  Experimental Current and simulated current is calculated
    ##  It is then subtracted from the FreqConstrainst number
    ##  If the value is less than 0 then these parameters are not allowed
    SimulatedCurrent = CurrentSimulation(ExperimentalPolarisationData[0], Parameters)   # Simulate the Current
    Constraint = CurrentConstraint - (abs(np.array(SimulatedCurrent) - np.array(ExperimentalPolarisationData[1])))  # Calculate the difference
    return Constraint 


def SingleVoltage_EISSimulation(Potential,
                                Parameters, RSol, Q, Alpha, FrequencyList, FirstSpeciesImp, SecondSpeciesImp):
    ## Caluclates EIS based on all the input values and return the total impedance in the form Zr + jZi
    
    ZtValues = []
    for frequency in FrequencyList:
        newZtValue = Ta_HFMechEIS(Potential,
                                    Parameters, RSol, Q, Alpha, frequency, FirstSpeciesImp, SecondSpeciesImp)
        ZtValues.append(newZtValue)

    Zt = np.array(ZtValues)
        
    return Zt


def CurrentSimulation(PotentialList, Parameters):
    ## Simulates all the current values for the given set of potential values
    CurrentList = []

    for Potential in PotentialList:
        
        NewSimulatedCurrent = Ta_HFMechDC(Potential, Parameters)
        CurrentList.append(NewSimulatedCurrent)

    return CurrentList


def Ta_HFMechDC(Potential,
                Parameters):

    k10         = 10**Parameters[0]
    b1          = Parameters[1]
    k20         = 10**Parameters[2]
    b2          = Parameters[3]
    k30prime    = 10**Parameters[4]
    k40prime    = 10**Parameters[5]
    b4          = Parameters[6]
    km10        = 10**Parameters[7]
    km20        = 10**Parameters[8]
    k400prime   = 10**Parameters[9]
    alpha       = Parameters[10]
    beta        = Parameters[11]
    delta       = Parameters[12]

    k30 = k30prime * ((SecondSpeciesPola)**alpha)
    k40 = k40prime * ((FirstSpeciesPola)**beta) + k400prime*((SecondSpeciesPola)**delta)

    k1 = k10*math.exp(b1*Potential)
    k2 = k20*math.exp(b2*Potential)
    k3 = k30
    k4 = k40*math.exp(b4*Potential)
    km1 = km10*math.exp(-(2*38-b1)*Potential)
    km2 = km20*math.exp(-(3*38-b2)*Potential)

    Theta1SS = (k1*(k3+km2))/((k1+km1+k2+k4)*k3+(k1-km2)*k2)
    Theta2SS = Theta1SS * k2 / (k3+km2)
    ThetaVSS = 1-Theta1SS-Theta2SS

    CurrentValue = 5*F*(((k2 * Theta1SS) + (k4 * Theta1SS) - (km2*Theta2SS)))

    return CurrentValue
    
    

def Ta_HFMechEIS(Potential,
                 Parameters, RSol, Q, Alpha, Freq, FirstSpeciesImp, SecondSpeciesImp):

    k10         = 10**Parameters[0]
    b1          = Parameters[1]
    k20         = 10**Parameters[2]
    b2          = Parameters[3]
    k30prime    = 10**Parameters[4]
    k40prime    = 10**Parameters[5]
    b4          = Parameters[6]
    km10        = 10**Parameters[7]
    km20        = 10**Parameters[8]
    k400prime   = 10**Parameters[9]
    alpha       = Parameters[10]
    beta        = Parameters[11]
    delta       = Parameters[12]
    Tau         = 10**Parameters[13]
    
    omega = 2 * math.pi * Freq
    k30 = k30prime*((SecondSpeciesImp)**alpha) 
    k40 = k40prime*((FirstSpeciesImp)**beta) + k400prime*((SecondSpeciesImp)**delta)

    k1 = k10 *math.exp(b1*Potential)
    k2 = k20 *math.exp(b2*Potential)
    k3 = k30
    k4 = k40 *math.exp(b4*Potential)
    km1 = km10 *math.exp(-(2*38-b1)*Potential)
    km2 = km20 *math.exp(-(3*38-b2)*Potential)

    Theta1SS = (k1*(k3+km2))/((k1+km1+k2+k4)*k3+(k1-km2)*k2)
    Theta2SS = Theta1SS * k2 / (k3+km2)

    bm1=-(2*38-b1)
    bm2=-(3*38-b2)

    A1 = Tau*complex(0, omega)+k1+km1+k2+k4-km2
    B1 = k1-km2
    C1 = (b1*k1*(1-Theta1SS-Theta2SS)-bm1*km1*Theta1SS-b2*k2*Theta1SS-b4*k4*Theta1SS+bm2*km2*Theta2SS)
    A2 = k2
    B2 = -(k3+Tau*complex(0, omega)+km2)
    C2 = -(b2*k2*Theta1SS)+km2*bm2*Theta2SS

    dtheta1 = (C2*B1-B2*C1)/(B1*A2-B2*A1)
    dtheta2 = (C2-(A2*dtheta1))/B2

    RtInv = F * (2*(b1*k1*(1-Theta1SS-Theta2SS)-bm1*km1*Theta1SS) + 3*(b2*k2*Theta1SS -km2*bm2*Theta2SS+ b4 * k4 *Theta1SS))

    Rt = 1 / RtInv

    ZfInv = RtInv - F *((2* (k1+km1)- 3 * (k2+k4)) * dtheta1+ (2 * k1 +3*km2)* dtheta2)

    
    Zt = (RSol + (1 /(ZfInv + Q*complex(0, omega)**Alpha)))
    

    return Zt
        
   
def OptimizationProgram(NumberOfTrials):
    ###Optimisation Loop###
    ParameterVariability = 10
    Parameters = np.array(InitialParameters)
    MinError = InitialError

    ## Loop through number of trials, set new parameter values based on jiggled values then run the optimisation again.
    for i in range(NumberOfTrials):
        print(f'Trial {i+1} of {NumberOfTrials}')    
    
        NewParameters = Parameters * (1 + (np.random.random_sample(len(InitialParameters))-0.5) * ParameterVariability/100) # Add Vairability

        Solution = minimize(CalculateError, NewParameters, method='SLSQP', bounds = AllBounds, constraints=Constraints, options={'maxiter':1000})
        print(f'New Error = {Solution.fun}')
        print(Solution)
        
        if Solution.fun < MinError:
            print('New Error is Less')
            MinError = Solution.fun
            Parameters = Solution.x
        else:
            print('New Error is More')
            
        print(f'Current Error = {MinError}')
    
    print(f'MinError = {MinError}')
    FinalParameters = Parameters
    print(f'Final Parameters = {FinalParameters}')

    return FinalParameters


def SaveParameters(Parameters, FileName):
    with open(FileName, mode='w') as File:
        CSVWriter = csv.writer(File, delimiter=',')

        CSVWriter.writerow(Parameters)


def LoadParameters(FileName):
    with open(FileName, mode='r') as File:
        CSVReader = csv.reader(File)
        LineCount = 0

        ParameterList = []
        for row in CSVReader:
            if LineCount == 0:
                for Parameter in row:
                    ParameterList.append(float(Parameter))
                    LineCount =+ 1
            else:
                return ParameterList

            
def PotentialPlot(FinalParameters):
    plt.plot(ExperimentalPolarisationData[0], np.array(ExperimentalPolarisationData[1])*1000)   # Turn into a np.array and multiply by 1000 to get into mA
    plt.plot(ExperimentalPolarisationData[0], np.array(CurrentSimulation(ExperimentalPolarisationData[0], FinalParameters))*1000, color="orange")   # Turn into a np.array and multiply by 1000 to get into mA
    plt.xlabel('Potential (V vs OCP)')
    plt.ylabel('Current Density (mA/cm^2)')
    plt.show()

    
def EISPlot(FinalParameters):
    for i in range(3):

        Zt = SingleVoltage_EISSimulation(PotentialSet[i],FinalParameters, RSolSet[i], QSet[i], AlphaSet[i], FrequencySet[i], FirstSpeciesImpVector[i], SecondSpeciesImpVector[i])

        plt.plot(ExperimentalEISData[i][1], ExperimentalEISData[i][2])
        plt.plot(np.array(Zt.real), np.array(Zt.imag))
        plt.xlabel('Zre (Ohm/cm2)')
        plt.ylabel('Zimag (Ohm/cm2)')
        plt.gca().invert_yaxis()
        plt.show()
        plt.cla()

        
###Constants###
F = 96485 # Faradays Constant


###Initialise Files
PotentialCurrentFile = 'Ta_0.5MHF_PDP.txt'

EISFiles = [
            'Ta_500mMHF_50mV_EIS.txt',
            'Ta_500mMHF_200mV_EIS.txt',
            'Ta_500mMHF_300mV_EIS.txt',
            ]


###Data From Other Programs###
PotentialSet = [0.05, 0.2, 0.3]   # List of Potentials the EIS was captued at
RSolSet = [5.6, 5.6, 5.6]   # Solution Resistance Gained from CPP and EEC
QSet = [9.87e-5, 5.95e-5, 5.58e-5]   # Effective Capacitance Gained from CPP and EEC
AlphaSet = [0.83, 0.91, 0.89]   # Capacitance Alpha Value Gained from CPP and EEC
RtVectorSet = [206.6, 296.6, 429]   # Charge Transfer Resistance Gained from CPP and EEC



###Initial Values and Upper/lower Boundaries for Parameters###
TInitP = [  #TemporaryInitialParameters                
            4.1295985e-09,   # k10
            1.2014046e+01,   # b1
            2.4842893e-11,   # k20
            1.5717619e+01,   # b2
            9.5701493e-11,   # k30Prime
            1.2789314e-08,   # k40Prime
            1.1370702e+01,   # b4
            4.1301441e-08,   # km10
            2.4842893e-11,   # km20
            1.4463068e-09,   # k400
            5.8208725e-01,   # alpha
            1.4773024e+00,   # beta
            1.7900219e+00,   # delta
            2.6710180e-09    # Gamma
        ]

InitialParameters = [ # Multiply by logs so that when the number is changed it does it on a better scale.
                    math.log10(TInitP[0]),  # k10
                    TInitP[1],              # b1
                    math.log10(TInitP[2]),  # k20
                    TInitP[3],              # b2
                    math.log10(TInitP[4]),  # k30prime
                    math.log10(TInitP[5]),  # k40prime
                    TInitP[6],              # b4
                    math.log10(TInitP[7]),  # km10
                    math.log10(TInitP[8]),  # km20
                    math.log10(TInitP[9]),  # k400
                    TInitP[10],             # alpha
                    TInitP[11],             # beta
                    TInitP[12],             # delta
                    math.log10(TInitP[13])  # Gamma
                    ]

UpperBoundParameters = [
                        math.log10(1e-2),   #k10
                        2*38,               #b1
                        math.log10(1e-2),   #k20
                        3*38,               #b2
                        math.log10(1e-2),   #k30prime
                        math.log10(1e-2),   #k40prime
                        3*38,               #b4
                        math.log10(1e-2),   #km10
                        math.log10(1e-2),   #km10
                        math.log10(1e-2),   #k400
                        3,                  #alpha
                        3,                  #beta
                        3,                  #delta
                        math.log10(1e-2),   #Gamma
                        ]

LowerBoundParameters = [
                        math.log10(1e-20),   #k10
                        0,               #b1
                        math.log10(1e-20),   #k20
                        0,               #b2
                        math.log10(1e-20),   #k30prime
                        math.log10(1e-20),   #k40prime
                        0,               #b4
                        math.log10(1e-20),   #km10
                        math.log10(1e-25),   #km10
                        math.log10(1e-25),   #k400
                        0,                  #alpha
                        0,                  #beta
                        0,                  #delta
                        math.log10(1e-30),   #Gamma
                        ]

# Put all bounds together so that the minimise function can read it
AllBounds = [[LowerBoundParameters[0], UpperBoundParameters[0]],
             [LowerBoundParameters[1], UpperBoundParameters[1]],
             [LowerBoundParameters[2], UpperBoundParameters[2]],
             [LowerBoundParameters[3], UpperBoundParameters[3]],
             [LowerBoundParameters[4], UpperBoundParameters[4]],
             [LowerBoundParameters[5], UpperBoundParameters[5]],
             [LowerBoundParameters[6], UpperBoundParameters[6]],
             [LowerBoundParameters[7], UpperBoundParameters[7]],
             [LowerBoundParameters[8], UpperBoundParameters[8]],
             [LowerBoundParameters[9], UpperBoundParameters[9]],
             [LowerBoundParameters[10], UpperBoundParameters[10]],
             [LowerBoundParameters[11], UpperBoundParameters[11]],
             [LowerBoundParameters[12], UpperBoundParameters[12]],
             [LowerBoundParameters[13], UpperBoundParameters[13]],
            ]


###Extract Data###

#Extract Polarisation Data
ExperimentalPolarisationData = OpenPotentialCurrentDataFile(PotentialCurrentFile)   #[[Potential, Potential, ...], [Current, Current, ...]]

# Extract EIS Data
ExperimentalEISData = OpenEISDataFiles(EISFiles)    #[Frequency, Zreal, Zimag]
ZDataSet = ExtractDataFromEISDataSet(ExperimentalEISData)   #Gets all Impedances together and converts into complex number
FrequencySet = [ExperimentalEISData[0][0], ExperimentalEISData[1][0], ExperimentalEISData[2][0]] #Gets all frequencies together


###Species Concentration###
# First Species concentration for polarisation analysis.  This is "X" in the equation 15 of the case study
FirstSpeciesPola=0.019541                    
# First Species concentration for Impedance analysis corrsponding to each Vdc
FirstSpeciesImpVector=[0.019541, 0.019541, 0.019541]    

# Second Species concentration for polarisation analysis. % This is "Y" in equation 16 of the case study
SecondSpeciesPola=0.43711
# Second Species concentration for polarisation analysis. 
SecondSpeciesImpVector=[0.43711, 0.43711, 0.43711] 


###Constraints###
CurrentConstraint = 0.000015

# No Equality Constraints used
# Only In-Equality Constraints used
Constraints = {'type':'ineq', 'fun':NonLinearConstraint}


###Weighing Factors###
# Impedance weighing Factors, frequency is in the denominator, and is raised to the power given here
IWF = [0, 0, 0]

# Weighing Factors for Edc data
DataWeighingFactor = np.array([0.5, 0.1, 0.1])
ImpedanceWeightFactor = CalculateImpedanceWeightFactor(PotentialSet, FrequencySet, IWF)


###Calculate Initial Error###
InitialError = CalculateError(InitialParameters)
print(f'Initial Error: {InitialError}')


###Optimisation Program###
Run = False # If True - Optimisation runs, if False - Optimisation does not run
NumberOfTrials = 10

if Run:
    FinalParameters = OptimizationProgram(NumberOfTrials)

    SaveParameters(FinalParameters, 'FinalParameters10It.csv') # WARNING - this will override any previously saved parameters unless the filename is changed.



###Plotting Graphs###
FinalParameters = LoadParameters('FinalParameters10It.csv')

PotentialPlot(FinalParameters)
EISPlot(FinalParameters)





























