
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "surfaceWaveVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "dimensionedVector.H"
#include "volMesh.H"
#include "IFstream.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceWaveVelocityFvPatchVectorField::surfaceWaveVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
	mixedFvPatchField<vector>(p, iF),
	g_(vector::zero),
	startPoint_(vector::zero),
	refU_(vector::zero),
	waveType_("regular"),
	ampl_(0),
	k_(0),
	wavePeriod_(0),
	phase_(0),
	xDir(vector::zero),
	zDir(vector::zero),
	xComp(-1),
	zComp(-1),
	seaLevel_(0),
	distance(0),
	omega_(0)
{
    this->refValue() = pTraits<vector>::zero;
//    this->value() = pTraits<vector>::zero;
    this->valueFraction() = 0.0;
}


surfaceWaveVelocityFvPatchVectorField::surfaceWaveVelocityFvPatchVectorField
(
    const surfaceWaveVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
	mixedFvPatchField<vector>(ptf, p, iF, mapper),
	g_(ptf.g_),
	startPoint_(ptf.startPoint_),
	refU_(ptf.refU_),
	waveType_(ptf.waveType_),
	ampl_(ptf.ampl_),
	k_(ptf.k_),
	wavePeriod_(ptf.wavePeriod_),
	phase_(ptf.phase_),
	xDir(ptf.xDir),
	zDir(ptf.zDir),
	xComp(ptf.xComp),
	zComp(ptf.zComp),
	seaLevel_(ptf.seaLevel_),
	distance(ptf.distance),
	omega_(ptf.omega_)
{}


surfaceWaveVelocityFvPatchVectorField::surfaceWaveVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
	mixedFvPatchField<vector>(p, iF),
	g_(vector::zero),
	startPoint_(vector::zero),
	refU_(vector::zero),
	waveType_("regular"),
//	ampl_(0),
	k_(0),
	wavePeriod_(0),
	phase_(0),
	xDir(vector::zero),
	zDir(vector::zero),
	xComp(-1),
	zComp(-1),
	seaLevel_(0),
	distance(0),
	omega_(0)
{
//    const uniformDimensionedVectorField& g =
//        db().lookupObject<uniformDimensionedVectorField>("g");
	uniformDimensionedVectorField g
	(
	    IOobject
	    (
		  "g",
		  this->db().time().constant(),
		  this->db(),
		  IOobject::MUST_READ,
		  IOobject::NO_WRITE
	    )
	);    

	IOdictionary waveDict_
	(
	    IOobject
	    (
		  "waveDict",
		  this->db().time().constant(),
		  this->db(),
		  IOobject::MUST_READ,
		  IOobject::NO_WRITE
	    )
	);    
	g_ = g.value();
	startPoint_=vector(waveDict_.lookup("startPoint"));
	refU_=vector(waveDict_.lookup("refU"));
	waveType_=word(waveDict_.lookup("waveType"));

	if (waveType_ 	== "regular")
	{
		const dictionary& waveDict = waveDict_.subDict("regularCoeffs");
		ampl_ 		= 0.5 * scalarField(1,dimensionedScalar(waveDict.lookup("waveHeight")).value());
		wavePeriod_     = scalarField(1,dimensionedScalar(waveDict.lookup("wavePeriod")).value());
		phase_		= scalarField(1,0);
		xDir 		= waveDict.lookup("waveTravelDirection");
	}
	else if (waveType_ == "irregular")
	{
		const dictionary& waveDict = waveDict_.subDict("irregularCoeffs");
		fileName waveGroupFileName(waveDict.lookup("waveGroupFileName"));
		scalar maxOmega(2.5);
		scalar n(100);
		scalar dOmega = maxOmega / n;
		IFstream waveStream(waveGroupFileName);
		if (waveStream.good())
		{
			List<vector> waveGroup
			(
				waveStream
			);
			ampl_.setSize(waveGroup.size());
			wavePeriod_.setSize(waveGroup.size());
			phase_.setSize(waveGroup.size());
			forAll(waveGroup, i)
			{
				ampl_[i] 	= Foam::sqrt(2*waveGroup[i][0]*dOmega);
				wavePeriod_[i] 	= 2*mathematicalConstant::pi/waveGroup[i][1];
				phase_[i] 	= waveGroup[i][2];
			}
			xDir 			= waveDict.lookup("waveTravelDirection");
		}	

	}
	solveDispersionRelation();
	determineDirections();
	seaLevel_=startPoint_[zComp];
	evaluate();
}


surfaceWaveVelocityFvPatchVectorField::surfaceWaveVelocityFvPatchVectorField
(
    const surfaceWaveVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
	mixedFvPatchField<vector>(ptf, iF),
	g_(ptf.g_),
	startPoint_(ptf.startPoint_),
	refU_(ptf.refU_),
	waveType_(ptf.waveType_),
	ampl_(ptf.ampl_),
	k_(ptf.k_),
	wavePeriod_(ptf.wavePeriod_),
	phase_(ptf.phase_),
	xDir(ptf.xDir),
	zDir(ptf.zDir),
	xComp(ptf.xComp),
	zComp(ptf.zComp),
	seaLevel_(ptf.seaLevel_),
	distance(ptf.distance),
	omega_(ptf.omega_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void surfaceWaveVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    scalar UX(0), UZ(0);
#   include "interfaceWeight.H"

    mixedFvPatchField<vector>::updateCoeffs();
}

// Evaluate the field on the patch
void surfaceWaveVelocityFvPatchVectorField::evaluate()
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    
    Field<vector>::operator=
    (
        this->refValue()
    );

    fvPatchField<vector>::evaluate();
}


// Write
void surfaceWaveVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("phi")
         << "phi" << token::END_STATEMENT << nl;
    this->refValue().writeEntry("value", os);
}

void surfaceWaveVelocityFvPatchVectorField::determineDirections()
{
	zDir = -1 * g_ / mag(g_);

	if (g_[0]==0 && g_[1]==0 && g_[2] !=0)
	{
	zComp = 2;
	}
	else if (g_[1]==0 && g_[2]==0 && g_[0] !=0)
	{
	zComp = 0;
	}
	else if (g_[0]==0 && g_[2]==0 && g_[1] !=0)
	{
	zComp = 1;
	}
	else
	{
        FatalErrorIn("surfaceWaveVelocityVectorField.C : ")
            << "gravity vector must be orthogonal... "
            << abort(FatalError);
	}

	if (xDir[zComp] != 0)
	{
	FatalErrorIn("surfaceWaveVelocityVectorField.C : ")
            << "wave travel direction must be perpendicular to gravity vector... "
            << abort(FatalError);
	}
	else
	{
	xDir /= mag(xDir);
	}
}

void surfaceWaveVelocityFvPatchVectorField::readIrregularCoeffs()
{

}

void surfaceWaveVelocityFvPatchVectorField::solveDispersionRelation()
{
	omega_ = 2*mathematicalConstant::pi / wavePeriod_;
	k_ = sqr(omega_) / mag(g_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, surfaceWaveVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
