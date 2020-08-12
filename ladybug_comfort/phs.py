import enum
import math


class BodyPosition(enum.Enum):
    """Effective radiating area of the body (1 = sitting, 2 = standing, 3 = crouching)"""

    sitting = 1
    standing = 2
    crouching = 3


def predictedHeatStrain(
    Ta: float,  # air temperature
    mrt: float,  # mean radiant temperature
    Tdp: float,  # dewpoint temperature
    wind_speed: float,  # in m/s
    SR: float,  # solar radiation, W/m2
    vapour_pressure_hPa: float,
    heightM: float,  # in m
    weight: float,  # in kg
    body_position: BodyPosition,  # 1 = sitting, 2 = standing, 3 = crouching
    insulation: float,  # clothing insulation factor 0-1
    metabolic_rate: float,  # in Watts?
    activityDuration: int,  # activity duration in minutes
    can_drink: bool = True,
    acclimatization: float = 100,  # how acclimatized is the subject? 0-100 where 100 is fully acclimatized
    walk_speed: float = 0,  # in m/s
    use_walk_speed: bool = True,  # if walking speed entered, 0 otherwise
    walk_angle: float = 0,  # angle between walking direction and wind direction in degrees
    use_walk_angle: bool = False,
    reflective_clothing: float = 0.54,  # fraction of the body surface covered by the reflective clothing, dimensionless
    reflective_clothing_emissivity: float = 0.97,  # emissivity of the reflective clothing, dimensionless
    work: float = 0,  # effective mechanical power, in W/m2
    imst: float = 0.38,  # static moisture permeability index, dimensionless
):
    """
    Calculate predicted heat strain from inputs, as per ISO 7933.

    Ta: float,  # air temperature
    mrt: float,  # mean radiant temperature
    Tdp: float,  # dewpoint temperature
    wind_speed: float, # in m/s
    SR: float,  # solar radiation, W/m2
    vapour_pressure_hPa: float,
    heightM: float, # in m
    weight: float, # in kg
    body_position: BodyPosition, # 1 = sitting, 2 = standing, 3 = crouching
    insulation: float,  # clothing insulation factor 0-1
    metabolic_rate: float,  # in Watts?
    activityDuration: float,  # activity duration in minutes
    can_drink: bool = True,
    acclimatization: float = 100,  # how acclimatized is the subject? 0-100 where 100 is fully acclimatized
    walk_speed: float = 0,  # in m/s
    use_walk_speed: bool = True,  # if walking speed entered, 0 otherwise
    walk_angle: float = 0,  # angle between walking direction and wind direction in degrees
    use_walk_angle: bool = False,
    reflective_clothing: float = 0.54,  # fraction of the body surface covered by the reflective clothing, dimensionless
    reflective_clothing_emissivity: float = 0.97,  # emissivity of the reflective clothing, dimensionless
    work: float = 0,  # effective mechanical power, in W/m2
    imst: float = 0.38,  # static moisture permeability index, dimensionless

    adapted from ladybug-legacy https://github.com/ladybug-tools/ladybug-legacy/blob/master/src/Ladybug_Thermal%20Comfort%20Indices.py#L1323-L1609
    based on: Dr. Jacques Malchaire Quick Basic code from:
    "Ergonomics of the thermal environment - Analytical determination and interpretation of heat stress using calculation of predicted heat strain", ISO 7933, 2004
    """

    if acclimatization > 100 or acclimatization < 0:
        raise ValueError("acclimatization out of range")

    # Effective radiating area of the body (1 = sitting, 2 = standing, 3 = crouching)
    ardu: float = 0
    if body_position is BodyPosition.sitting:
        ardu = 0.7  # dimensionless
    elif body_position is BodyPosition.standing:
        ardu = 0.77  # dimensionless
    elif body_position is BodyPosition.crouching:
        ardu = 0.67  # dimensionless
    else:
        raise ValueError("body_position invalid")

    if vapour_pressure_hPa < 0:
        raise ValueError("vapour_pressure_hPa cannot be negative")
    vapour_pressure_Pa = (
        vapour_pressure_hPa * 0.1
    )  # partial water vapour pressure, converted from hectopascals to kilopascals

    if reflective_clothing > 1 or reflective_clothing < 0:
        raise ValueError("reflective_clothing invalid")
    if reflective_clothing_emissivity > 1 or reflective_clothing_emissivity < 0:
        raise ValueError("reflective_clothing_emissivity invalid")

    # walking
    if work < 0:
        raise ValueError("work cannot be < 0")
    if imst < 0 or imst > 1:
        raise ValueError("imst invalid")

    if weight < 0:
        raise ValueError("weight cannot be < 0")
    if heightM < 0:
        raise ValueError("height cannot be < 0")
    body_surface_area = (
        0.202 * weight ** 0.425 * heightM ** 0.725
    )  # body surface area in m2
    spHeat = 57.83 * weight / body_surface_area

    SWp: float = 0
    SWtot: float = 0
    temperature_rectal = 36.8
    Tcr = 36.8
    Tsk = 34.1
    Tcreq = 36.8
    TskTcrwg = 0.3
    Dlimtre: float = 0
    Dlimloss50: float = 0
    Dlimloss95: float = 0
    Dmax50 = 0.075 * weight * 1000
    Dmax95 = 0.05 * weight * 1000

    # EXPONENTIAL AVERAGING CONSTANTS
    # Core temperature as a function of the metabolic rate: time constant: 10 minutes
    ConstTeq = math.exp(-1 / 10)
    # Skin Temperature: time constant: 3 minutes
    ConstTsk = math.exp(-1 / 3)
    # Sweat rate: time constant: 10 minutes
    ConstSW = math.exp(-1 / 10)

    for time in range(
        1, activityDuration + 1
    ):  # activityDuration - the duration of the work sequence in minutes
        # INITIALISATION MIN PER MIN
        Tsk0 = Tsk
        temperature_rectal0 = temperature_rectal
        Tcr0 = Tcr
        Tcreq0 = Tcreq
        TskTcrwg0 = TskTcrwg

        # EVALUATION OF THE MAXIMUM SWEAT RATE AS A FUNCTION OF THE METABOLIC RATE
        SWmax: float = (metabolic_rate - 32) * body_surface_area
        if SWmax > 400:
            SWmax = 400
        if SWmax < 250:
            SWmax = 250
        # For acclimatizationimatised subjects (acclimatization=100), the maximum Sweat Rate is greater by 25%
        # is there any reason for it to be a float not a bool?
        if acclimatization >= 50:
            SWmax = SWmax * 1.25
        if acclimatization < 50:
            Wmax = 0.85
        else:
            Wmax = 1

        # EQUILIBRIUM CORE TEMPERATURE ASSOCIATED TO THE METABOLIC RATE
        Tcreqm = 0.0036 * metabolic_rate + 36.6
        # Core temperature at this minute, by exponential averaging
        Tcreq = Tcreq0 * ConstTeq + Tcreqm * (1 - ConstTeq)
        # Heat storage associated with this core temperature increase during the last minute
        dStoreq = spHeat * (Tcreq - Tcreq0) * (1 - TskTcrwg0)

        # SKIN TEMPERATURE PREDICTION
        # Skin Temperature in equilibrium
        # Clothed model
        Tskeqcl = (
            12.165
            + 0.02017 * Ta
            + 0.04361 * mrt
            + 0.19354 * vapour_pressure_Pa
            - 0.25315 * wind_speed
        )
        Tskeqcl = Tskeqcl + 0.005346 * metabolic_rate + 0.51274 * temperature_rectal
        # Nude model
        Tskeqnu = (
            7.191
            + 0.064 * Ta
            + 0.061 * mrt
            + 0.198 * vapour_pressure_Pa
            - 0.348 * wind_speed
        )
        Tskeqnu = Tskeqnu + 0.616 * temperature_rectal
        # Value at this minute, as a function of the clothing insulation
        if insulation > 1 or insulation < 0:
            raise ValueError("insulation invalid")
        if insulation >= 0.6:
            Tskeq = Tskeqcl
        else:
            Tskeq = Tskeqnu + 2.5 * (Tskeqcl - Tskeqnu) * (
                insulation - 0.2
            )  # Interpolation between the values for clothed and nude subjects, if  0.2 < clo < 0.6

        if insulation <= 0.2:
            Tskeq = Tskeqnu
        else:
            Tskeq = Tskeqnu + 2.5 * (Tskeqcl - Tskeqnu) * (
                insulation - 0.2
            )  # Interpolation between the values for clothed and nude subjects, if  0.2 < clo < 0.6

        # Skin Temperature at this minute, by exponential averaging
        Tsk = Tsk0 * ConstTsk + Tskeq * (1 - ConstTsk)
        # Saturated water vapour pressure at the surface of the skin
        Psk = 0.6105 * math.exp(17.27 * Tsk / (Tsk + 237.3))

        # CLOTHING INFLUENCE ON EXCHANGE COEFFICIENTS
        # Static clothing insulation
        insulation_static = insulation * 0.155
        # Clothing area factor
        fcl = 1 + 0.3 * insulation
        # Static boundary layer thermal insulation in quiet air
        Iast = 0.111
        # Total static insulation
        Itotst = insulation_static + Iast / fcl

        # Relative velocities due to air velocity and movements
        if use_walk_speed:
            if use_walk_angle:
                # Unidirectional walking
                air_flow = abs(
                    wind_speed - walk_speed * math.cos(3.14159 * walk_angle / 180)
                )
            else:
                # Omni-directional walking
                if wind_speed < walk_speed:
                    air_flow = walk_speed
                else:
                    air_flow = wind_speed
        else:
            # Stationary or undefined speed
            walk_speed = 0.0052 * (metabolic_rate - 58)
            if walk_speed > 0.7:  # how acclimatized is the subject? 0-100
                walk_speed = 0.7
            air_flow = wind_speed

        # Dynamic clothing insulation
        # Clothing insulation correction for wind (air_flow) and walking (walk_speed)
        Vaux = air_flow
        if air_flow > 3:
            Vaux = 3
        Waux = walk_speed
        if walk_speed > 1.5:
            Waux = 1.5
        CORcl = 1.044 * math.exp(
            (0.066 * Vaux - 0.398) * Vaux + (0.094 * Waux - 0.378) * Waux
        )
        if CORcl > 1:
            CORcl = 1
        CORia = math.exp(
            (0.047 * air_flow - 0.472) * air_flow + (0.117 * Waux - 0.342) * Waux
        )
        if CORia > 1:
            CORia = 1

        CORtot = CORcl
        if insulation <= 0.6:
            CORtot = ((0.6 - insulation) * CORia + insulation * CORcl) / 0.6

        Itotdyn = Itotst * CORtot
        IAdyn = CORia * Iast
        insulation_dyn = Itotdyn - IAdyn / fcl

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
        Texp = 28.56 + 0.115 * Ta + 0.641 * vapour_pressure_Pa
        Cres = 0.001516 * metabolic_rate * (Texp - Ta)
        Eres = (
            0.00127 * metabolic_rate * (59.34 + 0.53 * Ta - 11.63 * vapour_pressure_Pa)
        )

        # Mean temperature of the clothing: Tcl
        # Dynamic convection coefficient
        Z = 3.5 + 5.2 * air_flow
        if air_flow > 1:
            Z = 8.7 * air_flow ** 0.6
        Hcdyn = 2.38 * abs(Tsk - Ta) ** 0.25
        if Z > Hcdyn:
            Hcdyn = Z

        auxR = 5.67e-08 * ardu
        FclR = (
            1 - reflective_clothing
        ) * 0.97 + reflective_clothing * reflective_clothing_emissivity
        Tcl = mrt + 0.1

        for iteration in range(100):
            # Radiation coefficient
            Hr = FclR * auxR * ((Tcl + 273) ** 4 - (mrt + 273) ** 4) / (Tcl - mrt)
            Tcl1 = ((fcl * (Hcdyn * Ta + Hr * mrt) + Tsk / insulation_dyn)) / (
                fcl * (Hcdyn + Hr) + 1 / insulation_dyn
            )

            if abs(Tcl - Tcl1) > 0.001:
                Tcl = (Tcl + Tcl1) / 2
                continue
            else:
                break

        # Convection and Radiation heat exchanges
        Conv = fcl * Hcdyn * (Tcl - Ta)
        Rad = fcl * Hr * (Tcl - mrt)
        # Maximum Evaporation Rate
        Emax = (Psk - vapour_pressure_Pa) / Rtdyn
        # Required Evaporation Rate
        Ereq = metabolic_rate - dStoreq - work - Cres - Eres - Conv - Rad

        # INTERPRETATION
        # Required wettedness
        wreq = Ereq / Emax

        # Required Sweat Rate
        #    If no evaporation required: no sweat rate
        SWreq: float = 0
        if Ereq <= 0:
            Ereq = 0
            SWreq = 0
        else:
            #    If evaporation is not possible, sweat rate is maximum
            if Emax <= 0:
                Emax = 0
                SWreq = SWmax
            else:
                #    If required wettedness greater than 1.7: sweat rate is maximum
                if wreq >= 1.7:
                    wreq = 1.7
                    SWreq = SWmax
                else:
                    #    Required evaporation efficiency
                    Eveff = 1 - wreq ** 2 / 2
                    if wreq > 1:
                        Eveff = (2 - wreq) ** 2 / 2
                    SWreq = Ereq / Eveff
                    if SWreq > SWmax:
                        SWreq = SWmax

        # Predicted Sweat Rate, by exponential averaging
        SWp = SWp * ConstSW + SWreq * (1 - ConstSW)
        if SWp <= 0:
            Ep: float = 0
            SWp = 0
        else:
            # Predicted Evaporation Rate
            k: float = Emax / SWp
            wp: float = 1
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
                Tcr1 = (Tcr1 + Tcr) / 2
                continue
            else:
                break

        # PREDICTION OF THE CENTRAL (RECTAL) TEMPERATURE
        temperature_rectal = (
            temperature_rectal0 + (2 * Tcr - 1.962 * temperature_rectal0 - 1.31) / 9
        )  # in Celsius degrees

        if Dlimtre == 0 and temperature_rectal >= 38:
            Dlimtre = time
        # Total water loss rate during the minute (in W/m2)
        SWtot = SWtot + SWp + Eres
        SWtotg = SWtot * 2.67 * body_surface_area / 1.8 / 60

        if Dlimloss50 == 0 and SWtotg >= Dmax50:
            Dlimloss50 = time
        if Dlimloss95 == 0 and SWtotg >= Dmax95:
            Dlimloss95 = time
        if can_drink == 0:
            Dlimloss95 = Dlimloss95 * 0.6
            Dlimloss50 = Dlimloss95
        continue

    if Dlimloss50 == 0:
        Dlimloss50 = activityDuration
    if Dlimloss95 == 0:
        Dlimloss95 = activityDuration
    if Dlimtre == 0:
        Dlimtre = activityDuration

    is_comfortable: bool = True
    effectPHS: float = 0

    if (Dlimloss95 >= activityDuration) and (Dlimtre >= activityDuration):
        effectPHS = 0
        is_comfortable = True
    elif (Dlimtre <= 30) or (Dlimloss95 <= 30):
        effectPHS = 4
        is_comfortable = False
    elif ((Dlimtre > 30) and (Dlimtre <= 120)) or (
        (Dlimloss95 > 30) and (Dlimloss95 <= 120)
    ):
        effectPHS = 3
        is_comfortable = False
    elif (
        (Dlimtre > 120) and (Dlimtre < (activityDuration - activityDuration * 0.015))
    ) or (
        (Dlimloss95 > 120)
        and (Dlimloss95 < (activityDuration - activityDuration * 0.015))
    ):
        effectPHS = 2
        is_comfortable = False
    elif (
        (Dlimtre >= (activityDuration - activityDuration * 0.015))
        and (Dlimtre < activityDuration)
    ) or (
        (Dlimloss95 >= (activityDuration - activityDuration * 0.015))
        and (Dlimloss95 < activityDuration)
    ):
        effectPHS = 1
        is_comfortable = False

    return temperature_rectal, effectPHS, is_comfortable
