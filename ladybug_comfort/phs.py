def predictedHeatStrain(Ta, mrt, Tdp, rh, ws, SR, N, Tground, Rprim, vapourPressure, Epot, age, sex, heightCM, heightM, weight, bodyPosition, Icl, ac, acclimated, Met, activityDuration, HRrates, dehydrationRiskRates, climate):
    # inputs: (Ta, mrt, ws, vapourPressure, heightM, weight, bodyPosition, Icl, acclimated, Met, activityDuration):
    # based on: Dr. Jacques Malchaire Quick Basic code from:
    # "Ergonomics of the thermal environment - Analytical determination and interpretation of heat stress using calculation of predicted heat strain", ISO 7933, 2004
    
    drink = 1  # if drink = 1: water replacement is sufficient, and workers can drink freely; otherwise: drink = 0
    if acclimated == "acclimated":  # 100 if acclimatised person, 0 if unacclimatised person
        accl = 100
    elif acclimated == "unacclimated":
        accl = 0
    
    # Effective radiating area of the body (1 = sitting, 2 = standing, 3 = crouching)
    if bodyPosition == "sitting":
        Ardu = 0.7  # dimensionless
    elif bodyPosition == "standing":
        Ardu = 0.77  # dimensionless
    elif bodyPosition == "crouching":
        Ardu = 0.67  # dimensionless
    
    Pa = vapourPressure * 0.1  # partial water vapour pressure, converted from hectopascals to kilopascals
    Work = 0  # effective mechanical power, in W/m2
    imst = 0.38  # static moisture permeability index, dimensionless
    Ap = 0.54  # fraction of the body surface covered by the reflective clothing, dimensionless
    Fr = 0.97  # emissivity of the reflective clothing, dimensionless
    # walking
    defspeed = 1  # 1 if walking speed entered, 0 otherwise
    Walksp = 0  # walking speed, in m/s
    defdir = 1  # 1 if walking direction entered, 0 otherwise
    THETA = 0  # angle between walking direction and wind direction degrees
    
    Adu = 0.202 * weight ** 0.425 * heightM ** 0.725  # body surface area in m2
    spHeat = 57.83 * weight / Adu
    
    SWp = 0
    SWtot = 0; Tre = 36.8; Tcr = 36.8; Tsk = 34.1; Tcreq = 36.8; TskTcrwg = 0.3
    Dlimtre = 0; Dlimloss50 = 0; Dlimloss95 = 0
    Dmax50 = 0.075 * weight * 1000
    Dmax95 = 0.05 * weight * 1000
    
    # EXPONENTIAL AVERAGING CONSTANTS
    # Core temperature as a function of the metabolic rate: time constant: 10 minutes
    ConstTeq = math.exp(-1 / 10)
    # Skin Temperature: time constant: 3 minutes
    ConstTsk = math.exp(-1 / 3)
    # Sweat rate: time constant: 10 minutes
    ConstSW = math.exp(-1 / 10)
    
    for time in range(1, activityDuration+1):  # activityDuration - the duration of the work sequence in minutes
        # INITIALISATION MIN PER MIN
        Tsk0 = Tsk; Tre0 = Tre; Tcr0 = Tcr; Tcreq0 = Tcreq; TskTcrwg0 = TskTcrwg
        
        # EVALUATION OF THE MAXIMUM SWEAT RATE AS A FUNCTION OF THE METABOLIC RATE
        SWmax = (Met - 32) * Adu
        if SWmax > 400: SWmax = 400
        if SWmax < 250: SWmax = 250
        # For acclimatised subjects (accl=100), the maximum Sweat Rate is greater by 25%
        if accl >= 50: SWmax = SWmax * 1.25
        if accl < 50: Wmax = 0.85
        else: Wmax = 1
        
        # EQUILIBRIUM CORE TEMPERATURE ASSOCIATED TO THE METABOLIC RATE
        Tcreqm = 0.0036 * Met + 36.6
        # Core temperature at this minute, by exponential averaging
        Tcreq = Tcreq0 * ConstTeq + Tcreqm * (1 - ConstTeq)
        # Heat storage associated with this core temperature increase during the last minute
        dStoreq = spHeat * (Tcreq - Tcreq0) * (1 - TskTcrwg0)
        
        # SKIN TEMPERATURE PREDICTION
        # Skin Temperature in equilibrium
        # Clothed model
        Tskeqcl = 12.165 + 0.02017 * Ta + 0.04361 * mrt + 0.19354 * Pa - 0.25315 * ws
        Tskeqcl = Tskeqcl + 0.005346 * Met + 0.51274 * Tre
        # Nude model
        Tskeqnu = 7.191 + 0.064 * Ta + 0.061 * mrt + 0.198 * Pa - 0.348 * ws
        Tskeqnu = Tskeqnu + 0.616 * Tre
        # Value at this minute, as a function of the clothing insulation
        if Icl >= 0.6: Tskeq = Tskeqcl
        else: Tskeq = Tskeqnu + 2.5 * (Tskeqcl - Tskeqnu) * (Icl - 0.2)  # Interpolation between the values for clothed and nude subjects, if  0.2 < clo < 0.6
        
        if Icl <= 0.2: Tskeq = Tskeqnu
        else: Tskeq = Tskeqnu + 2.5 * (Tskeqcl - Tskeqnu) * (Icl - 0.2)  # Interpolation between the values for clothed and nude subjects, if  0.2 < clo < 0.6
        
        # Skin Temperature at this minute, by exponential averaging
        Tsk = Tsk0 * ConstTsk + Tskeq * (1 - ConstTsk)
        # Saturated water vapour pressure at the surface of the skin
        Psk = 0.6105 * math.exp(17.27 * Tsk / (Tsk + 237.3))
        
        # CLOTHING INFLUENCE ON EXCHANGE COEFFICIENTS
        # Static clothing insulation
        Iclst = Icl * 0.155
        # Clothing area factor
        fcl = 1 + 0.3 * Icl
        # Static boundary layer thermal insulation in quiet air
        Iast = 0.111
        # Total static insulation
        Itotst = Iclst + Iast / fcl
        
        # Relative velocities due to air velocity and movements
        if defspeed > 0:
            if defdir == 1:
                # Unidirectional walking
                Var = abs(ws - Walksp * math.cos(3.14159 * THETA / 180))
            else:
                # Omni-directional walking
                if ws < Walksp:
                    Var = Walksp
                else: Var = ws
        else:
            # Stationary or undefined speed
            Walksp = 0.0052 * (Met - 58)
            if Walksp > 0.7:
                Walksp = 0.7
            Var = ws
        
        # Dynamic clothing insulation
        # Clothing insulation correction for wind (Var) and walking (Walksp) 
        Vaux = Var
        if Var > 3:
            Vaux = 3
        Waux = Walksp
        if Walksp > 1.5:
            Waux = 1.5
        CORcl = 1.044 * math.exp((.066 * Vaux - 0.398) * Vaux + (.094 * Waux - 0.378) * Waux)
        if CORcl > 1:
            CORcl = 1
        CORia = math.exp((.047 * Var - 0.472) * Var + (.117 * Waux - 0.342) * Waux)
        if CORia > 1:
            CORia = 1
         
        CORtot = CORcl
        if Icl <= 0.6:
            CORtot = ((.6 - Icl) * CORia + Icl * CORcl) / .6
        
        Itotdyn = Itotst * CORtot
        IAdyn = CORia * Iast
        Icldyn = Itotdyn - IAdyn / fcl
        
        # Permeability index
        # Correction for wind and walking
        CORe = (2.6 * CORtot - 6.5) * CORtot + 4.9
        imdyn = imst * CORe
        if imdyn > 0.9:
            imdyn = 0.9
        # Dynamic evaporative resistance
        Rtdyn = Itotdyn / imdyn / 16.7
        
        # HEAT EXCHANGES
        # Heat exchanges through respiratory convection and evaporation
        # temperature of the expired air
        Texp = 28.56 + 0.115 * Ta + 0.641 * Pa
        Cres = 0.001516 * Met * (Texp - Ta)
        Eres = 0.00127 * Met * (59.34 + 0.53 * Ta - 11.63 * Pa)
        
        # Mean temperature of the clothing: Tcl
        # Dynamic convection coefficient
        Z = 3.5 + 5.2 * Var
        if Var > 1:
            Z = 8.7 * Var ** 0.6
        Hcdyn = 2.38 * abs(Tsk - Ta) ** 0.25
        if Z > Hcdyn:
            Hcdyn = Z
        
        auxR = 5.67E-08 * Ardu
        FclR = (1 - Ap) * 0.97 + Ap * Fr
        Tcl = mrt + 0.1
        
        for k in range(100):
            # Radiation coefficient
            Hr = FclR * auxR * ((Tcl + 273) ** 4 - (mrt + 273) ** 4) / (Tcl - mrt)
            Tcl1 = ((fcl * (Hcdyn * Ta + Hr * mrt) + Tsk / Icldyn)) / (fcl * (Hcdyn + Hr) + 1 / Icldyn)
            
            if abs(Tcl - Tcl1) > 0.001:
                Tcl = (Tcl + Tcl1) / 2
                continue
            else:
                break
        
        # Convection and Radiation heat exchanges
        Conv = fcl * Hcdyn * (Tcl - Ta)
        Rad = fcl * Hr * (Tcl - mrt)
        # Maximum Evaporation Rate
        Emax = (Psk - Pa) / Rtdyn
        # Required Evaporation Rate
        Ereq = Met - dStoreq - Work - Cres - Eres - Conv - Rad
         
        # INTERPRETATION
        # Required wettedness
        wreq = Ereq / Emax
        
        # Required Sweat Rate
        #    If no evaporation required: no sweat rate
        if Ereq <= 0:
            Ereq = 0; SWreq = 0;
        else:
            #    If evaporation is not possible, sweat rate is maximum
            if Emax <= 0:
                Emax = 0; SWreq = SWmax;
            else:
                #    If required wettedness greater than 1.7: sweat rate is maximum
                if wreq >= 1.7:
                    wreq = 1.7; SWreq = SWmax;
                else:
                    #    Required evaporation efficiency
                    Eveff = (1 - wreq ** 2 / 2)
                    if wreq > 1:
                        Eveff = (2 - wreq) ** 2 / 2
                    SWreq = Ereq / Eveff
                    if SWreq > SWmax:
                        SWreq = SWmax
        
        # Predicted Sweat Rate, by exponential averaging
        SWp = SWp * ConstSW + SWreq * (1 - ConstSW)
        if SWp <= 0:
            Ep = 0; SWp = 0;
        else:
            # Predicted Evaporation Rate
            k = Emax / SWp
            wp = 1
            if k >= 0.5:
                wp = -k + math.sqrt(k * k + 2)
            if wp > Wmax:
                wp = Wmax
            Ep = wp * Emax
        
        # Heat Storage
        dStorage = Ereq - Ep + dStoreq
        
        # PREDICTION OF THE CORE TEMPERATURE
        Tcr1 = Tcr0
        for g in range(50):
            # Skin - Core weighting
            TskTcrwg = 0.3 - 0.09 * (Tcr1 - 36.8)
            if TskTcrwg > 0.3:
                TskTcrwg = 0.3
            if TskTcrwg < 0.1:
                TskTcrwg = 0.1
            
            Tcr = dStorage / spHeat + Tsk0 * TskTcrwg0 / 2 - Tsk * TskTcrwg / 2
            Tcr = (Tcr + Tcr0 * (1 - TskTcrwg0 / 2)) / (1 - TskTcrwg / 2)
            if abs(Tcr - Tcr1) > 0.001:
                Tcr1 = (Tcr1 + Tcr) / 2;
                continue
            else:
                break
        
        # PREDICTION OF THE CENTRAL (RECTAL) TEMPERATURE
        Tre = Tre0 + (2 * Tcr - 1.962 * Tre0 - 1.31) / 9  # in Celsius degrees
        
        if Dlimtre == 0 and Tre >= 38: Dlimtre = time
        # Total water loss rate during the minute (in W/m2)
        SWtot = SWtot + SWp + Eres
        SWtotg = SWtot * 2.67 * Adu / 1.8 / 60
        
        if Dlimloss50 == 0 and SWtotg >= Dmax50: Dlimloss50 = time
        if Dlimloss95 == 0 and SWtotg >= Dmax95: Dlimloss95 = time
        if drink == 0:
            Dlimloss95 = Dlimloss95 * 0.6;
            Dlimloss50 = Dlimloss95
        continue
    
    if Dlimloss50 == 0:
        Dlimloss50 = activityDuration
    if Dlimloss95 == 0:
        Dlimloss95 = activityDuration
    if Dlimtre == 0:
        Dlimtre = activityDuration
    
    PHSresults = [SWtotg, Dlimloss95, Dlimtre]
    
    
    if (Dlimloss95 >= activityDuration) and (Dlimtre >= activityDuration):
        effectPHS = 0
        comfortable = 1
    elif (Dlimtre <= 30) or (Dlimloss95 <= 30):
        effectPHS = 4
        comfortable = 0
    elif ((Dlimtre > 30) and (Dlimtre <= 120)) or ((Dlimloss95 > 30) and (Dlimloss95 <= 120)):
        effectPHS = 3
        comfortable = 0
    elif ((Dlimtre > 120) and (Dlimtre < (activityDuration - activityDuration*0.015))) or ((Dlimloss95 > 120) and (Dlimloss95 < (activityDuration - activityDuration*0.015))):
        effectPHS = 2
        comfortable = 0
    elif ((Dlimtre >= (activityDuration - activityDuration*0.015)) and (Dlimtre < activityDuration)) or ((Dlimloss95 >= (activityDuration - activityDuration*0.015)) and (Dlimloss95 < activityDuration)):
        effectPHS = 1
        comfortable = 0
    
    return Tre, effectPHS, comfortable, []
