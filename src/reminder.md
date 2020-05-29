# REMINDER 

## Acronyms
| Acronym | Full Name |
| :-----: |:---------:|
| RRE     | Regular Refraction with reflected Expansion |
| RRR     | Regular Refraction with Reflected shock |
| BPR     | Bound Precursor Refraction |
| FNR     | Free precursor von Neumann Refraction |
| TNR     | Twin von Neumann Refraction |
| LSR     | Lambda Shock Refraction |
| FPR     | Free Precursor Refraction |

## Matlab programs

1. `rrrtorreBoundary.m` file

Inputs : xi, w, gamma1, gamma2, mu1, mu2

Analytically computes the value of the boundary function between a RRE system and a RRR system. When this value is zero, then the transition occurs

2. `isBeforeSfRRRToBPR.m` file

Inputs : xi, w, gamma1, gamma2, mu1, mu2 (+ temperature ratio if needed)

Determines whether the refraction system is regular or not by comparing the polars of incident and transmitted shock waves. When the trasmitted shock wave polar has points inside the incident one, then the refraction is irregular
Global RRR -> BPR transition line is obtained by a dichotomy method in the program computing those boundaries.

3. `bPRFNROmega.m` file

Inputs : xi, gamma1, gamma2, mu1, mu2 (+ temperature ratio if needed)

Analytically computes the w value for which BPR -> FNR transition occurs, at a given xi

Uses `xiToSqMach.m` to get value of the Mach number of the shock

4. `postShockMachSq.m` file

Inputs : xi, M1, gamma

Analytically computes the value of Mach number after the shock of xi pressure ratio. When the Mach number after the reflected wave equals 1, then TNR -> LSR transition occurs.

Is used in `isSfStrongRRE.m`

5. `isBeforeFPRToTNR.m` file

Inputs : xi, w, gamma1, gamma2, mu1, mu2 (+ temperature ratio if needed)
Determines whether the refraction system is FPR or TNR by using the method developped by Abd-el-Fattah and Henderson in thein 1978 article.
Global FPR -> TNR transition line is obtained by a dichotomy method in the program computing those boundaries.

6. `deltaMax.m` file

Inputs : M1 pre-shock Mach, gamma of the current phase

Calculates the maximum deflection angle after a shock

Uses `xiDeltaMax.m` to get
Uses `tanDefSq.m` to get

Is used in 

7. `getExpPolar.m` file

Inputs : Me, pre-expansion Mach, gamma of the current phase, pressure jump and deviation before expansion

Calculates the polar of the expansion fan by computing each deviation for a range of xi

Uses `prandtlMeyer.m` to get the deviation angle knowing the Mach number before and after the expansion and gamma of the current phase

Is used in `isSfStrongRRE.m`

8. `getPolar.m` file

Inputs : M1, pre-shock Mach, gamma of the current phase, mode of calculation (+pressure jump and deviation before shock if needed)

Calculates the polar of the any shock by computing each deviation for a range of xi. Different modes to get different shocks.

Uses `tanDefSq.m` to get the deviation angle knowing the Mach number before the shock, the pressure jump and gamma of the current phase

Is used in `isSfStrongRRE.m`

9. `getPolarPoint.m` file

Inputs : xis and deltas of a certain polar, mode of computation, coordinates of the point (+ boolean to know whether it is an expansion or not, if needed)

Returns the two closest values of a delta if mode = 0 (xi is given) or the two closest values of a xi if mode = 1 (delta is given)

Is used in `isSfStrongRRE.m`

10. `isSfStrongRRE.m` file

Inputs : xi, w, gamma1, gamma2, mu1, mu2 (+ temp_ratio if needed)

Determines whether the refraction structure is RRE or not by computing the graphical condition exhibited by Abd-el-Fattah and Henderson on polar shock intersections.

Uses `xiToSqMach.m` to calculate several Mach numbers
Uses `tanDefSq.m` to calculates deviations at several shocks
Uses `postShockMachSq.m` to calculate Mach number after the incident shock
Uses `getExpPolar.m`, `getPolar.m` and `getPolarPoint.m` to detrmine the intersection of the expansion polar and the transmitted polar

11. `machToXi.m` file

Inputs : Mach of the shock, gamma of the current phase, phi shock - pre-shock flow angle

Returns pressure jump through the shock

12. `plotPolar.m` file

Inputs : M1, pre-shock Mach, gamma of the current phase, mode of calculation (+pressure jump and deviation before shock if needed)

Computes the xis and deltas of a polar so that it can be directly plotted

13. `prandtlMeyer.m` file

Inputs : Mach number, gamma of the current phase

Computes the Prandtl-Meyer function for expansion waves

14. `tanDefSq.m` file

Inputs : xi pressure jump, mach number before shock, gamma of the current phase

Calculates the square of the tangent of deviation after shock

Is used in `isSfStrongRRE.m`

15. `xiDeltaMax.m` file

Inputs : M1 pre shock Mach, gamma of the current phase

Calculates the pressure jump through shock when the deflection angle is max

Is used in `getPolar.m`

16. `xiLim.m` file

Inputs : M1 pre shock Mach number, gamma of the current phase

Calculates the maximum pressure jump through shock

Is used in `getPolar.m`

17. `xiToSqMach.m` file

Inputs : pressure jump, gamma of the current phase, phi angle between shock and pre shock flow

Calculates sqaure of pre shock Mach number

Is used in `isSfStrongRRE.m`

18. `xiToTempj.m` file

Inputs : pressure jump, gamma of the current phase

Calculates the temperature jump through the shock
