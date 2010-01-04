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

Class
    surfaceWaveVelocityFvPatchVectorField

Description

SourceFiles
    surfaceWaveVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceWaveVelocityFvPatchVectorField_H
#define surfaceWaveVelocityFvPatchVectorField_H

#include "mixedFvPatchField.H"
#include "mathematicalConstants.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class inletOutletFvPatch Declaration
\*---------------------------------------------------------------------------*/
//template<class Type>

class surfaceWaveVelocityFvPatchVectorField
:
    public mixedFvPatchField<vector>
{
    // Private data

        //- Value field
//        vectorField value_;

protected:

    // Protected data
        vector g_;
	vector startPoint_;
	vector refU_;
	word waveType_;
        scalarField ampl_;
        scalarField k_;
        scalarField wavePeriod_;
	scalarField phase_;
	vector xDir, zDir;
	int xComp, zComp;
        scalar seaLevel_;
	scalar distance;
	scalarField omega_;
    
public:

    //- Runtime type information
    TypeName("surfaceWaveVelocity");


    // Constructors

        //- Construct from patch and internal field
        surfaceWaveVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surfaceWaveVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given surfaceWaveVelocityFvPatchVectorField onto a new patch
        surfaceWaveVelocityFvPatchVectorField
        (
            const surfaceWaveVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new surfaceWaveVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        surfaceWaveVelocityFvPatchVectorField
        (
            const surfaceWaveVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone(const DimensionedField<vector, volMesh>& iF) const
        {
            return tmp<fvPatchField<vector> >
            (
                new surfaceWaveVelocityFvPatchVectorField(*this, iF)
            );
        }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
        
        //- Evaluate the patch field
        virtual void evaluate();

        //- Write
        virtual void write(Ostream&) const;

	// determine wave elevation direction
	void determineDirections();

	void solveDispersionRelation();

	void readIrregularCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "surfaceWaveVelocityFvPatchVectorField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
