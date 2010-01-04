const vectorField & c = patch().Cf(); // Er der overhovedet brug for denne?
const fvMesh& mesh = this->dimensionedInternalField().mesh(); // Gives a reference to the mesh
const pointField& points = mesh.points(); // gives a reference to the points in the mesh

const label patchI = mesh.boundaryMesh().findPatchID( this->patch().name()); // ID for the inlet patch
const labelList cNumbers = this->patch().faceCells();  // IDs for the cells on next to the present patch

word faceType;
scalar interfaceRatio(0.0), maxZ(0.0), minZ(0.0);


forAll(cNumbers, cI)
{
	const cell& cProp = mesh.cells()[cNumbers[cI]]; // a cell is a list of faces
	forAll (cProp, cJ)
	{
		// This gives you the patch ID for a face. If the face is internal, you get -1
//		label patchID = mesh.boundaryMesh().whichPatch(cProp[cJ]); 
		if (mesh.boundaryMesh().whichPatch(cProp[cJ]) == patchI)
		{
		const face& fProp = mesh.faces()[cProp[cJ]]; // Gives a list of points, which the face consists of
//		face temp(fProp);
		maxZ = -1e100;
		minZ = 1e100;
		forAll(fProp, fI)
			{
			maxZ = max(maxZ, (points[fProp[fI]]& zDir));
			minZ = min(minZ, (points[fProp[fI]]& zDir));
			}
		}
	}

//Determine perpendicular horizontal distance from "startPoint": distance
//using xDir (which is normalized)
//cos angle = (u dot v)/(mag(u) * mag(v))
vector x_ = c[cI] - startPoint_;		//distance vector between current patch centre and startpoint
x_[zComp] = 0;				//perpendicular distance vector in case sealavel is not zero
distance = ((x_ - refU_ * db().time().value()) & xDir);	// /mag(xDir) which is 1; & = dot product

#include "evalEta.H"  // Her skal der fyldes lidt paa!
if (maxZ <= eta)
	{
	faceType = "WET";
	}
else if (minZ >= eta)
	{
	faceType = "DRY";
	}
else
	{
	faceType = "INTERFACE";
	interfaceRatio = mag((eta - minZ) / (maxZ - minZ));
	}

#   include "evaluateQuantity.H"
}

