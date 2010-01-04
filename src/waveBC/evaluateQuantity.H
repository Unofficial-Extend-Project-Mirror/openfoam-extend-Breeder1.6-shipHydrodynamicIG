/* This file sets the boundary condition for the velocity field
   based on area information for each and every boundary face on the wave
   maker patch.
   DRY: The face is above the generated wave.
   WET: The face is below the free surface.
   INTERFACE: The face has points on both sides of the free surface.
*/

    if (faceType == "DRY")
    {
        /* For the dry part, a wall condition is chosen */

        this->refValue()[cI] = refU_;
        this->valueFraction()[cI] = 1; 
    }
    else 
	{
        /* For the wet part the velocity from the adopted wave theory is imposed - 1st order Stokes */
        UX = sum(
//		mathematicalConstant::pi * ampl_ / wavePeriod_ 
//		* Foam::cosh(k_ * (c[cI].component(zComp) - seaLevel_ + depth_)) / Foam::sinh(k_ * depth_)  
		ampl_*omega_*Foam::exp(k_* (c[cI].component(zComp) - seaLevel_))
                * Foam::cos(k_*distance - omega_ * db().time().value() + phase_) 
		* softStartFactor);		//Deep water

        UZ = UX;//Deep water
//sum(
//		mathematicalConstant::pi * ampl_ / wavePeriod_
//                * Foam::sinh(k_ * (c[cI].component(zComp) - seaLevel_ + depth_)) / Foam::sinh(k_ * depth_)  
//		ampl_*Foam::exp(k_* (c[cI].component(zComp) - seaLevel_))
//                * Foam::sin(k_*distance - omega_ * db().time().value() + phase_) 
//		* softStartFactor);
	}

if (faceType == "WET")
    {
        this->valueFraction()[cI] = 1;
	this->refValue()[cI] = xDir * UX + zDir * UZ + refU_;
    }
    else if (faceType == "INTERFACE")
    {
        this->valueFraction()[cI] = 1;
	this->refValue()[cI] = interfaceRatio * (xDir * UX + zDir * UZ) + refU_;
    }



