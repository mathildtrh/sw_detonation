1. `rrrtorreBoundary.m` file

Inputs : xi, w, gamma1, gamma2, mu1, mu2

Analytically computes the value of the boundary function between a RRE system and a RRR system. When this value is zero, then the transition occurs

2. `isBeforeSfRRRToBPR.m` file

Inputs : xi, w, gamma1, gamma2, mu1, mu2 (+ temperature ratio if needed)

Determines whether the refraction system is regular or not by comparing the polars of incident and transmitted shock waves. When the trasmitted shock wave polar has points inside the incident one, then the refraction is irregular
Global transition line is obtained by a dichotomy method in the program computing those boundaries.

3. `bPRFNROmega.m` file

Inputs : xi, gamma1, gamma2, mu1, mu2 (+ temperature ratio if needed)

Analytically computes the w value  for which transition occurs, at a given xi

4. `postShockMachSq.m` file

Inputs : xi, M1, gamma

Analytically computes the value of Mach number after the shock of xi pressure ratio. When the Mach number after the reflected wave equals 1, then transition occurs.

5. `isBeforeFPRToTNR.m` file

Inputs : xi, w, gamma1, gamma2, mu1, mu2 (+ temperature ratio if needed)
Determines whether the refraction system is FPR or TNR by using the method developped by Abd-el-Fattah and Henderson in thein 1978 article.
Global transition line is obtained by a dichotomy method in the program computing those boundaries.

| Acronym | Full Name |
| :-----: |:---------:|
| RRE     | Regular Refraction with reflected Expansion |
| RRR     | Regular Refraction with Reflected shock |
| BPR     | Bound Precursor Refraction |
| FNR     | Free precursor von Neumann Refraction |
| TNR     | Twin von Neumann Refraction |
| LSR     | Lambda Shock Refraction |
| FPR     | Free Precursor Refraction |
