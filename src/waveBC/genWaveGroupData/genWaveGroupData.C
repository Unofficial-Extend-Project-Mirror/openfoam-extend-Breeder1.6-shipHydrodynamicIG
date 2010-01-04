/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

Application
    

Description
    

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "List.H"
#include "vector.H"
#include "OFstream.H"
#include "Random.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createDynamicFvMesh.H"
	scalar pi_ = mathematicalConstant::pi;
	IOdictionary waveProperties
	(
		IOobject
		(
			"waveDict",
			runTime.constant(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	const dictionary& waveDict = waveProperties.subDict("irregularCoeffs");

	fileName waveGroupFileName(waveDict.lookup("waveGroupFileName"));
	word spectrum(waveDict.lookup("spectrum"));
	scalar H_13 = readScalar(waveDict.lookup("H_13"));		//Let op: H1/3 hier gebruitk als Hs!
	scalar T1 = readScalar(waveDict.lookup("T1"));		//Let op: Tp en T1!!!
	scalar n(100);
	scalar maxOmega(2.5);
	scalar dOmega = maxOmega / (n+1);
	scalar gamma(0), sigma(0), Tp(0), omegaP(0), A(0);
	//Prepare Spectrum
	scalarField omega(n,0);
	scalarField S(n,0);
	forAll (omega,i)
	{
		omega[i] = 0.5*dOmega + maxOmega * i / n;
	}
	if (spectrum == "Bretschneider")
	{
		S = 172.8*pow(H_13,2)*pow(T1,-4)*pow(omega,-5)*Foam::exp(-691.2*pow(T1,-4)*pow(omega,-4));
	}
	else if (spectrum == "JONSWAP")
	{	
		Tp = T1*1.200;
		omegaP = 2*pi_/Tp;	
		gamma = 3.3;	
		forAll (omega, i)
		{
			if (omega[i] > 2*pi_/Tp)
			{
				sigma = 0.09;
			}	
			else
			{
				sigma = 0.07;
			}	
//		Info << "i : " << i << "; omega(i) : " << omega[i] << "; sigma : " << sigma << nl;
		A = Foam::exp( - pow(    (omega[i]/omegaP - 1) / (sigma*Foam::sqrt(2.0))    ,2 ));
//	Info << A << nl;
		S[i] = 320*pow(H_13,2)*pow(Tp,-4)*pow(omega[i],-5)*Foam::exp(-1950*pow(Tp,-4)*pow(omega[i],-4))*Foam::pow(gamma,A);
		}
		
	}
	else if (spectrum == "Neuman")
	{

	}
    Random rand(clock::getTime());
	    
    List<vector>waveComponents(n);
        OFstream dataFile(waveGroupFileName);
	dataFile << n << nl << "(" << nl;

    forAll(waveComponents, i)
    {
//	if (Foam::sqr(S[i]) < 1E-20)
//	{
//		amplitude = 1E-5;
//	}
//	else
//	{
//		amplitude = Foam::sqrt(2*S[i]);
//	}
        waveComponents[i] = vector
        (
            S[i],
            omega[i],
            2*pi_*rand.scalar01()
        );
	dataFile << waveComponents[i] << nl;
    }

	dataFile << ")" << nl;
 //   {
 //       OFstream dataFile(waveGroupFileName);
 //       dataFile << waveComponents << nl;
 //   }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
