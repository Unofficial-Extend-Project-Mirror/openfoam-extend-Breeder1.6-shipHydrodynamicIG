case to demonstrate wave boundary condition.
Mind the outflow conditions.
inletOutlet on alpha1. If using zeroGradient, the tank slowly runs empty.
zeroGradient or buoyantPressure on p.

A velocity of 2 m/s is superimposed on the orbital velocity of the water (e.g. to simulate a velocity of a boat against the waves). In this case the outflow is always directed outwards.
If using a "base" velocity of e.g. 0, the outflow will be sometimes outward and sometimes inward. This will cause wavereflections.
A dedicated outflow BC is yet to be found/made for this.

