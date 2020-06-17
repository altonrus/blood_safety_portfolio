# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 17:50:45 2020

@author: alton
"""

#from scipy.stats import norm  # normal distribution for CIs
import pandas as pd
import numpy as np
#import time
# import math
# import sys
# sys.path.append('G:\My Drive\Blood Transfusion\Optimal portfolio\blood portfolio r project\NMC_sim')

"""
 Parameterize
"""

    
SCANDATdata = pd.read_csv('G:/My Drive/Blood Transfusion/Optimal portfolio/blood portfolio r project/NMC_sim/SCANDATA.csv')
# Define parameters as global variables
p_death_healthy = [0.00586390262469649, 0.000396284798625857, 0.00026205790345557, 
                        0.000196736858924851, 0.000158240465680137, 0.000150831911014393, 
                        0.00013492931611836, 0.000121125871373806, 0.000107538355223369, 
                        9.52993868850172E-05, 8.87355563463643E-05, 9.49519744608551E-05, 
                        0.000122189288958907, 0.000175166132976301, 0.000248543423367664, 
                        0.00032804548391141, 0.000409930275054649, 0.000501902832183987, 
                        0.000602643529418856, 0.000706460501533002, 0.000813615217339247, 
                        0.000914035015739501, 0.00099424587097019, 0.00104815745726228, 
                        0.00108298438135535, 0.00111166411079466, 0.00114323885645717, 
                        0.00117668299935758, 0.00121571309864521, 0.00125990901142359, 
                        0.00130612566135824, 0.00135253288317472, 0.00140105350874364, 
                        0.00145133980549872, 0.0015040592988953, 0.00156606524251401, 
                        0.00163467996753752, 0.00170090352185071, 0.00176273647230119, 
                        0.00182732986286283, 0.00190737284719944, 0.00201099156402051, 
                        0.0021355883218348, 0.0022804201580584, 0.00244479230605066, 
                        0.00262052961625159, 0.00282053602859378, 0.00306605803780258, 
                        0.00336857605725527, 0.00371979991905391, 0.00409019365906715, 
                        0.00447356421500444, 0.00489137275144458, 0.00534395314753056, 
                        0.00582200475037098, 0.00631902134045959, 0.00682523706927896, 
                        0.0073414696380496, 0.00787698198109865, 0.00844705477356911, 
                        0.00906533654779196, 0.00973428599536419, 0.0104441214352846, 
                        0.0111765656620264, 0.0119275134056807, 0.0127070974558592, 
                        0.0135565148666501, 0.0145145915448666, 0.0156336780637503, 
                        0.0169610474258661, 0.018596213310957, 0.0205059554427862, 
                        0.02256422303617, 0.0246546063572168, 0.0268777962774038, 
                        0.0293662250041962, 0.0322047881782055, 0.0353653393685818, 
                        0.0388779826462269, 0.0431636124849319, 0.047922533005476, 
                        0.0529773235321045, 0.0588391572237015, 0.0655613914132118, 
                        0.0729490742087364, 0.0811935216188431, 0.0908139199018478, 
                        0.101376831531525, 0.112927429378033, 0.12550188601017, 
                        0.139124721288681, 0.153806000947952, 0.169538646936417, 
                        0.186296001076698, 0.204029887914658, 0.222669392824173, 
                        0.242120549082756, 0.262267202138901, 0.28297284245491, 
                        0.304084002971649, 1]

# multiplier derived from SCANDAT by bin
surv_mult = pd.Series([-0.817572418110032, -0.813276839365234,
                            -0.797255332490873, -0.796849824185898, -0.79389313639905, -0.788511714568704,
                            -0.785958671736801, -0.778819817626416, -0.776996855799508, -0.77468002220649,
                            -0.773511707057567, -0.770957286490019, -0.768091531084629, -0.767861689953426,
                            -0.764936741259688, -0.761923124147346, -0.760928033300991, -0.757204208586802,
                            -0.752917080591849, -0.752746349586526, -0.75140021756903, -0.750034076821807,
                            -0.747413823684727, -0.743547158650365, -0.737859557401385, -0.737720943842025,
                            -0.737275764221944, -0.735683527739676, -0.735319502027239, -0.734885190074792,
                            -0.733364801161838, -0.732552162235254, -0.731571583485507, -0.7299628395368,
                            -0.727997197668908, -0.727705820084743, -0.725136499691371, -0.724432636786181,
                            -0.723955979443069, -0.721608410206325, -0.720420333352485, -0.719413966075051,
                            -0.714501341533066, -0.711670127108074, -0.710399485143846, -0.709618768639444,
                            -0.709577753842831, -0.709515845538552, -0.708520965144798, -0.707221322118684,
                            -0.705628462903537, -0.704279072938498, -0.699977726376569, -0.693087484256031,
                            -0.691340408861523, -0.681419540343925, -0.679104793883406, -0.670283295124512,
                            -0.661069074272555, -0.657565790561358, -0.646130514189894, -0.642647996239707,
                            -0.633325515629018, -0.631280742765057, -0.623842152698512, -0.596598212512847,
                            -0.594285903987529, -0.582534174204996, -0.581110787220343, -0.514117833707892,
                            -0.48637077676596, -0.449511020731904, -0.368944671501742, -0.366154771453697,
                            -0.292639044770105, -0.280746675672264, -0.153440011414076, -0.032664173435345,
                            0.518510271429626, 1.4651623712012, 1.89076717961945],
                           [3331, 3231, 3233, 3133, 3333, 3321, 3123, 3332, 3132, 3131,
                            3323, 3232, 2133, 3223, 2123, 2332, 2333, 2331, 3213, 3313,
                            2321, 2233, 2323, 2232, 2231, 1233, 2223, 2132, 2313, 1333,
                            2213, 3221, 1332, 3113, 1331, 1123, 1321, 2131, 3322, 1223,
                            1323, 2322, 1132, 1133, 1232, 2111, 1313, 1322, 1131, 1231,
                            1213, 3111, 2113, 2221, 1122, 1222, 1221, 1113, 3222, 2222,
                            3121, 3122, 1312, 2312, 2122, 3312, 2121, 1121, 1111, 2212,
                            1212, 3212, 1311, 1112, 2112, 3112, 2311, 3311, 3211, 2211,
                            1211])
surv_mult_max = max(surv_mult)

#  Kleinman unadjusted
surv_unadj = [[0.895, 0.877, 0.862, 0.862, 0.848],
                   [0.729, 0.667, 0.630, 0.601, 0.572],
                   [0.618, 0.507, 0.435, 0.400, 0.312]]
                   



"""
 Define utility sub-functions
"""

    
def serializeParams(params, indices):
    for key in indices.keys():
        params[key] = pd.Series(data = params[key], index = indices[key])
    return params

def discFac(t1, t2, params):
    return (1 - np.exp(-1 * params["rdisc"] * (t2 - t1))) / (params["rdisc"] * (t2 - t1))

#@profile
def getSurvival(age, unitsRBC, unitsPLT, unitsFFP, params):
    # Convert age, units RBC, FFP, PLT into bins to pull in survival multiplier
    if age < 45:
        agebin = 1
    elif age < 65:
        agebin = 2
    else:
        agebin = 3

    if unitsRBC == 0:
        RBCbin = 1
    elif unitsRBC < 5:
        RBCbin = 2
    else:
        RBCbin = 3

    if unitsPLT == 0:
        PLTbin = 1
    elif unitsPLT < 5:
        PLTbin = 2
    else:
        PLTbin = 3

    if unitsFFP == 0:
        FFPbin = 1
    elif unitsFFP < 5:
        FFPbin = 2
    else:
        FFPbin = 3
        

    lookup = agebin * 1000 + RBCbin * 100 + PLTbin * 10 + FFPbin
    mult = surv_mult[lookup]
    # print("lookup, mult", lookup, mult)
    # Generate survival cumprob
    n = 100 - age
    cumSurv = np.zeros(int(n + 1))
    healthySurv = 1 - p_death_healthy[age]
    recip_surv_unadj = surv_unadj[agebin - 1][0]
    #print(recip_surv_unadj)
    cumSurv[0] = recip_surv_unadj + (mult / surv_mult_max) * (healthySurv - recip_surv_unadj)

    # print("healthy, unadj, cumSurv[0]", healthySurv, surv_unadj, cumSurv[0])

    # Will first fill cumSurv with the probability of surviving to a year
    for i in range(1, 5):
        healthySurv = healthySurv * (1 - p_death_healthy[age + i])
        recip_surv_unadj = surv_unadj[agebin - 1][i]
        cumSurv[i] = recip_surv_unadj + (mult / surv_mult_max) * (healthySurv - recip_surv_unadj)
    for i in range(5, n):
        cumSurv[i] = cumSurv[i - 1] * (1 - p_death_healthy[age + i])

    # Then do '1-' in order to give the probability of dying in a given year
    cumSurv = 1 - cumSurv


    survival = next(x[0] for x in enumerate(cumSurv)
                         if x[1] > np.random.uniform()) + np.random.uniform()
    
    
    bl_QALY = survival * params["uBaseline"] * discFac(0, survival, params)
    
    return bl_QALY, survival
    # print("cumSurv", cumSurv)
    # print("survival", self.survival)
    
#@profile
def addProdCost(duration, ageStart, mort, params):
    #print("duration", duration)
    #If mortality, first do the consumption
    cLoss = 0 #tallies consumption loss
    if mort==1:
        yearsAdded = 0   #tracks years of productivity that have been added
        cDuration = duration #duration for adding consumption
        while cDuration > 0:
            currentBracket = params["consumpByAge"].index[params["consumpByAge"].index<=ageStart+yearsAdded].max()
            nextBracket = params["consumpByAge"].index[params["consumpByAge"].index>ageStart+yearsAdded].min()
            if ageStart + cDuration > nextBracket:
                timeLost = nextBracket - currentBracket
            else:
                timeLost = cDuration
            yearsAdded += timeLost
            cLoss += timeLost*params["consumpByAge"][currentBracket]*discFac(max(currentBracket, ageStart) - ageStart, yearsAdded, params)
            cDuration = cDuration - timeLost
        #print("Consumption lost", cLoss)
        
    
    yearsAdded = 0   #tracks years of productivity that ahve been added
    prodLoss = 0 #tracks productivity loss
    while duration > 0:
        currentBracket = params["prodByAge"].index[params["prodByAge"].index<=ageStart+yearsAdded].max()
        nextBracket = params["prodByAge"].index[params["prodByAge"].index>ageStart+yearsAdded].min()
        if ageStart + duration > nextBracket:
            timeLost = nextBracket - currentBracket
        else:
            timeLost = duration
        yearsAdded += timeLost
        prodLoss += timeLost*params["prodByAge"][currentBracket]*discFac(max(currentBracket, ageStart) - ageStart, yearsAdded, params)
        duration = duration - timeLost

        #print("current", currentBracket, "next", nextBracket,"timeLost", timeLost, "prodLoss", prodLoss)
    #print("ageStart", ageStart, "Mort", mort, "ProdLoss", prodLoss, "cLoss", cLoss)
    return prodLoss - cLoss

def getRecip(params):
    recipIndex = int(np.floor(np.random.uniform() * 790824))
    #print(self.SCANDATdata.iloc[[recipIndex]])
    age = int(SCANDATdata.iat[recipIndex, 4] + np.floor(np.random.uniform() * 5))
    # print("age", self.age)
    #sex = SCANDATdata.get_value(recipIndex, 0]
    unitsFFP = SCANDATdata.iat[recipIndex, 3]*0.408765293
    unitsRBC = SCANDATdata.iat[recipIndex, 1]*0.408765293
    unitsPLT = SCANDATdata.iat[recipIndex, 2]*0.408765293
    
    #print("Age, RBC, PLT, FFP", self.age, self.unitsRBC, self.unitsPLT, self.unitsFFP)
    bl_QALY, survival = getSurvival(age, unitsRBC, unitsPLT, unitsFFP, params)
    #print("Survival", self.survival)
    return age, unitsFFP, unitsRBC, unitsPLT, bl_QALY, survival

def recipOutcomes(age, unitsFFP, unitsRBC, unitsPLT, bl_QALY, survival, params):
    randRecip = np.random.uniform()
    randDeath = np.random.uniform();
    randDisability = np.random.uniform()
    costRecip = 0; costProd = 0;
    QALYLrecip = 0
    encephalitis = 0
    meningitis = 0
    AFP = 0
    fever = 0
    
    immunocomp_mult = 1
    if unitsPLT > 0:
        if np.random.uniform() < params["p_immunocomp_PLT"]:
            immunocomp_mult = params["RR_immunocomp"]
    else:
        if np.random.uniform() < params["p_immunocomp_noPLT"]:
            immunocomp_mult = params["RR_immunocomp"]
    
    if np.random.uniform() < params["p_symptoms"] * immunocomp_mult:
        # clinical symptoms  
        p_fev = params["p_fever"].asof(age)
        p_enceph = params["p_encephalitis"].asof(age)
        p_menin = params["p_meningitis"].asof(age)
        #params["p_AFP"] = params["p_AFP"].get_value(params["p_AFP"].index.asof(age))
        

        # WNV fever
        if randRecip < p_fev:
            fever = 1
            costRecip = params["cost_fever_init"] * discFac(0, min(survival, params["dFever"]), params)
            if randDeath < params["p_death_fever"].asof(age):
                QALYLrecip = bl_QALY
                costProd = addProdCost(survival, age, 1, params)
            else:
                QALYLrecip = min(survival, params["dFever"]) * params["uBaseline"] * (1 - params["uFever"]) * discFac(0, min(survival, params["dFever"]), params) 
                costProd = addProdCost(min(survival, params["dFever"]), age, 0, params)
                if survival > 0.5:
                    costRecip += params["cost_fever_5yr"] * (min(survival, 5.0)/5.0) * discFac(0, min(survival, 5.0), params)
                if age < 50:
                    if randDisability < params["p_disability_fever_u50"]:
                        QALYLrecip += min(survival, params["dDisability"]) * params["uBaseline"] * (1 - params["uFeverDisability"]) * discFac(0, min(survival, params["dDisability"]), params)
                else:
                    if randDisability < params["p_disability_fever_50up"]:
                         QALYLrecip += min(survival, params["dDisability"]) * params["uBaseline"] * (1 - params["uFeverDisability"]) * discFac(0, min(survival, params["dDisability"]), params)
        
        # Meningitis    
        elif randRecip < p_fev + p_menin:
             
            meningitis = 1
            costRecip = params["cost_mening_init"] * discFac(0, min(survival, params["dMeningitis"]), params)
            if randDeath < params["p_death_meningitis"].asof(age):
                QALYLrecip = bl_QALY
                costProd = addProdCost(survival, age, 1, params)
            else:
                QALYLrecip = min(survival, params["dMeningitis"]) * params["uBaseline"] * (1 - params["uMeningitis"]) * discFac(0, min(survival, params["dMeningitis"]), params) 
                costProd = addProdCost(min(survival, params["dMeningitis"]), age, 0, params)
                if survival > 0.5:
                    costRecip += params["cost_menin_5yr"] * (min(survival, 5.0)/5.0) * discFac(0, min(survival, 5.0), params)
                if randDisability < params["p_immunocomp_noPLT"]:
                    QALYLrecip += min(survival, params["dDisability"]) * params["uBaseline"] * (1 - params["uMeningitisDisability"]) * discFac(0, min(survival, params["dDisability"]), params)
    
        # Encephalitis    
        elif randRecip < p_fev + p_menin + p_enceph:
            encephalitis = 1
            costRecip = params["cost_enceph_init"] * discFac(0, min(survival, params["dEncephalitis"]), params)
            if randDeath < params["p_death_encephalitis"].asof(age):
                QALYLrecip = bl_QALY
                costProd = addProdCost(survival, age, 1, params)
            else:
                QALYLrecip = min(survival, params["dEncephalitis"]) * params["uBaseline"] * (1 - params["uEncephalitis"]) * discFac(0, min(survival, params["dEncephalitis"]), params) 
                costProd = addProdCost(min(survival, params["dEncephalitis"]), age, 0, params)
                if survival > 0.5:
                    costRecip += params["p_disability_encephalitis"] * (min(survival, 5.0)/5.0) * discFac(0, min(survival, 5.0), params)
                if randDisability < params["p_immunocomp_noPLT"]:
                    QALYLrecip += min(survival, params["dDisability"]) * params["uBaseline"] * (1 - params["uEncephalitisDisability"]) * discFac(0, min(survival, params["dDisability"]), params)

        # Acute flaccid paralysis    
        else:
            AFP = 1
            costRecip = params["cost_AFP_init"] * discFac(0, min(survival, params["dAFP"]), params)
            if randDeath < params["p_death_AFP"].asof(age):
                QALYLrecip = bl_QALY
                costProd = addProdCost(survival, age, 1, params)
            else:
                QALYLrecip = min(survival, params["dAFP"]) * params["uBaseline"] * (1 - params["uAFP"]) * discFac(0, min(survival, params["dAFP"]), params) 
                costProd = addProdCost(min(survival, params["dEncephalitis"]), age, 0, params)
                if survival > 0.5:
                    costRecip += params["cost_AFP_5yr"] * (min(survival, 5.0)/5.0) * discFac(0, min(survival, 5.0), params)
                if randDisability < params["p_disability_AFP"]:
                    QALYLrecip += min(survival, params["dDisability"]) * params["uBaseline"] * (1 - params["uAFPDisability"]) * discFac(0, min(survival, params["dDisability"]), params)

    return costRecip, costProd, QALYLrecip, fever, meningitis, encephalitis, AFP



def doRep(params):# Runs simulation for one recipient

        #Get next recipient
        age, unitsFFP, unitsRBC, unitsPLT, bl_QALY, survival = getRecip(params)

        # Get prob pregnant
        # Get outcomes assuming no pregnancy
        costRecip, costProd, QALYLrecip, fever, meningitis, encephalitis, AFP = recipOutcomes(age, unitsFFP, unitsRBC, unitsPLT, bl_QALY, survival, params)
        
        return (costRecip + costProd, QALYLrecip, fever, 
                meningitis, encephalitis, AFP, 
                unitsFFP, unitsRBC, unitsPLT)



"""

Run model 

"""
# def test(obj):
#     print(obj)
#     print(type(obj))


def runmodelWNV_separate_indices(N_iter, N_rep, params, indices, fname = None, to_file = 1):
  params = serializeParams(params, indices)
  runModelWNV(N_iter, N_rep, params, fname = None, to_file = 1)



def runModelWNV(N_iter, N_rep, params, fname = None, to_file = 1):
    
    output_cols = ["Tot_RBC",
                   "Tot_PLT",
                   "Tot_FFP",
                   "cost_times_RBC",
                   "cost_times_PLT",
                   "cost_times_FFP",
                   "QALY_times_RBC",
                   "QALY_times_PLT",
                   "QALY_times_FFP",
                   "fever_times_RBC",
                   "fever_times_PLT",
                   "fever_times_FFP",
                   "mening_times_RBC",
                   "mening_times_PLT",
                   "mening_times_FFP",
                   "enceph_times_RBC",
                   "enceph_times_PLT",
                   "enceph_times_FFP",
                   "AFP_times_RBC",
                   "AFP_times_PLT",
                   "AFP_times_FFP",
                   "n_recip",
                   "iter"
                   ]

        
    rows = np.zeros(shape=(N_iter, 23))
    
    for itn in range(0 , N_iter):
        cost, QALYs, fever, meningitis, encephalitis, AFP, unitsFFP, unitsRBC, unitsPLT = doRep(params)
    
        row = np.array([unitsRBC, unitsPLT, unitsFFP, 
                                    cost*unitsRBC, cost*unitsPLT, cost*unitsFFP,
                                    QALYs*unitsRBC, QALYs*unitsPLT, QALYs*unitsFFP,
                                    fever*unitsRBC, fever*unitsPLT, fever*unitsFFP,
                                    meningitis*unitsRBC, meningitis*unitsPLT, meningitis*unitsFFP,
                                    encephalitis*unitsRBC, encephalitis*unitsPLT, encephalitis*unitsFFP,
                                    AFP*unitsRBC, AFP*unitsPLT, AFP*unitsFFP])
        
        for reps in range(1, N_rep):
            cost, QALYs, fever, meningitis, encephalitis, AFP, unitsFFP, unitsRBC, unitsPLT = doRep(params)
            row  += np.array([unitsRBC, unitsPLT, unitsFFP,
                                    cost*unitsRBC, cost*unitsPLT, cost*unitsFFP,
                                    QALYs*unitsRBC, QALYs*unitsPLT, QALYs*unitsFFP,
                                    fever*unitsRBC, fever*unitsPLT, fever*unitsFFP,
                                    meningitis*unitsRBC, meningitis*unitsPLT, meningitis*unitsFFP,
                                    encephalitis*unitsRBC, encephalitis*unitsPLT, encephalitis*unitsFFP,
                                    AFP*unitsRBC, AFP*unitsPLT, AFP*unitsFFP])
        rows[itn, ] = np.append(row, [N_rep,  itn])
        #print("iter ", itn, " complete")
        
        
      

    output = pd.DataFrame(data=rows, columns=output_cols)
      
    if to_file == 1:
        output.to_csv(fname, sep=',')
    else:
        return(output)







def PSA_on_Sherlock(iter_start, n_iters, N_rep, fname):
  # Read in PSA input
  inputs = pd.read_csv("NMC_sim/output/WNV_PSA_inputs.csv", index_col = "key")
  
  # Loop over indivated sets of inout
  for iter in range(iter_start, iter_start + n_iters):
    # Construct the parameter dictionary
    param_dict = dict()
    for param_key in  param_keys:
      if inputs.loc[param_key,"V"+ str(iter)].size == 1:
        param_dict[param_key] = inputs.loc[param_key,"V"+ str(iter)]
      else:
        param_dict[param_key] = pd.Series(data = inputs.loc[param_key,"V"+ str(iter)].array, index =  inputs.loc[param_key,"Index"].array)
        
    #Run the model
    this_output = runModelWNV(1, N_rep, params, to_file = o)
    
    #Concatinate the iterations
    if iter == iter_start:
      agg_output = this_output
    else:
      agg_output.append(this_output)
  #Write to file
  agg_output.to_csv(fname, sep=',')
  
  
  

  
  







