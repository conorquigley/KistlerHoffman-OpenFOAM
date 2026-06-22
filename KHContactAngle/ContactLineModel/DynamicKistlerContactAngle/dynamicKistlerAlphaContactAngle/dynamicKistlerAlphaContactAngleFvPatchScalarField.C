/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
      \\/    M anipulation  |
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

#include "dynamicKistlerAlphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "fvPatchFields.H"
#include "volMesh.H"
#include "Time.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "IOdictionary.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const scalar dynamicKistlerAlphaContactAngleFvPatchScalarField::convertToDeg =
    180.0/constant::mathematical::pi;

const scalar dynamicKistlerAlphaContactAngleFvPatchScalarField::convertToRad =
    constant::mathematical::pi/180.0;

const scalar dynamicKistlerAlphaContactAngleFvPatchScalarField::theta0 = 85.5;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicKistlerAlphaContactAngleFvPatchScalarField::
dynamicKistlerAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    dynamicAlphaContactAngleFvPatchScalarField(p, iF),
    thetaA_(0.0),
    thetaR_(0.0),
    muWater_(0.0)
{}

dynamicKistlerAlphaContactAngleFvPatchScalarField::
dynamicKistlerAlphaContactAngleFvPatchScalarField
(
    const dynamicKistlerAlphaContactAngleFvPatchScalarField& acpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    dynamicAlphaContactAngleFvPatchScalarField(acpsf, p, iF, mapper),
    thetaA_(acpsf.thetaA_),
    thetaR_(acpsf.thetaR_),
    muWater_(acpsf.muWater_)
{}


dynamicKistlerAlphaContactAngleFvPatchScalarField::
dynamicKistlerAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    dynamicAlphaContactAngleFvPatchScalarField(p, iF),
    thetaA_(readScalar(dict.lookup("thetaA"))),
    thetaR_(readScalar(dict.lookup("thetaR"))),
    muWater_(readScalar(dict.lookup("muWater")))
{
    evaluate();
}


dynamicKistlerAlphaContactAngleFvPatchScalarField::
dynamicKistlerAlphaContactAngleFvPatchScalarField
(
    const dynamicKistlerAlphaContactAngleFvPatchScalarField& acpsf
)
:
    dynamicAlphaContactAngleFvPatchScalarField(acpsf),
    thetaA_(acpsf.thetaA_),
    thetaR_(acpsf.thetaR_),
    muWater_(acpsf.muWater_)
{}


dynamicKistlerAlphaContactAngleFvPatchScalarField::
dynamicKistlerAlphaContactAngleFvPatchScalarField
(
    const dynamicKistlerAlphaContactAngleFvPatchScalarField& acpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    dynamicAlphaContactAngleFvPatchScalarField(acpsf, iF),
    thetaA_(acpsf.thetaA_),
    thetaR_(acpsf.thetaR_),
    muWater_(acpsf.muWater_)
{}


//* * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> dynamicKistlerAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
  const scalar sigmap =
      this->db().lookupObject<IOdictionary>("transportProperties")
      .get<dimensionedScalar>("sigma").value();
//--------------------------------------------------------------------------------------------------Access the layers above the boundary layer-----------------------------


// Access mesh
const fvMesh& mesh = patch().boundaryMesh().mesh();
const labelUList& faceCells = patch().faceCells(); // cells adjacent to boundary
const vector wallNormalDir(0, 0, 1); // "up" = +Z

// How many layers to skip (1 = first internal cell, 2 = second, etc.)
const label nLayers = 3;

// Build a list of cell indices N layers into the domain
labelList targetCells(patch().size());
forAll(faceCells, pfacei)
{
    label cellI = faceCells[pfacei];
    for (label layer = 0; layer < nLayers; layer++)
    {
        // Find the internal face of cellI pointing away from the wall
        // by picking the face whose normal best aligns with -nf
        const cell& c = mesh.cells()[cellI];
        scalar maxDot = -GREAT;
        label bestFace = -1;
        forAll(c, cFacei)
        {
            label faceI = c[cFacei];
            if (mesh.isInternalFace(faceI))
            {
                vector fNorm = mesh.faceAreas()[faceI] / mag(mesh.faceAreas()[faceI]);
                // nf points outward (into wall), so interior direction is -nf
                scalar d = fNorm & wallNormalDir;
                if (d > maxDot)
                {
                    maxDot = d;
                    bestFace = faceI;
                }
            }
        }
        // Walk to the neighbour cell across bestFace
        if (bestFace != -1)
        {
            label own = mesh.faceOwner()[bestFace];
            label nei = mesh.faceNeighbour()[bestFace];
            cellI = (own == cellI) ? nei : own;
        }
    }
    targetCells[pfacei] = cellI;
}

//----------------------------Deep wall velocity

vectorField nf = patch().nf();
const volVectorField& U = mesh.lookupObject<volVectorField>("U");

vectorField UwallDeep(patch().size());
forAll(targetCells, pfacei)
{
    UwallDeep[pfacei] = U[targetCells[pfacei]];
}

// Remove wall-normal component (same as before)
UwallDeep -= (nf & UwallDeep) * nf;

// X-component (or swap back to nWall dot product if preferred)
scalarField uwall = -UwallDeep.component(vector::X);

// * * * * * * * * * * * Kistler Conact Angle Computation * * * * * * * * //

    //eb - Calculate local Capillary number
    scalarField Ca = muWater_*mag(uwall)/sigmap;

    //eb - Instantiate function object InverseHoffmanFunction for thetaA and thetaR
    dynamicKistlerAlphaContactAngleFvPatchScalarField::InverseHoffmanFunction
    InvHoffFuncThetaA
    (
        convertToRad*thetaA_
    );

    dynamicKistlerAlphaContactAngleFvPatchScalarField::InverseHoffmanFunction
    InvHoffFuncThetaR
    (
        convertToRad*(thetaR_)//calculate the receding contact angl
    );

    //eb - Calculate InverseHoffmanFunction for thetaA and thetaR using
    // RiddersRoot
    RiddersRoot RRInvHoffFuncThetaA(InvHoffFuncThetaA, 1.e-10);
    scalar InvHoffFuncThetaAroot = RRInvHoffFuncThetaA.root(0,150);

    RiddersRoot RRInvHoffFuncThetaR(InvHoffFuncThetaR, 1.e-10);
    scalar InvHoffFuncThetaRroot = RRInvHoffFuncThetaR.root(0,150);

    //eb - Calculate and return the value of contact angle on patch faces,
    //     a general approach: the product of Uwall and nWall is negative
    //     for advancing and positiv for receding motion.
    //     thetaDp is initialized to 90 degrees corresponding to no wall
    //     adhesion

//---------------------------------------Kistler--------------
const scalar uwallTol = 1e-4;

    scalarField thetaDp(patch().size(), convertToRad*theta0);
    forAll(uwall, pfacei)
    {

        if (uwall[pfacei] < -uwallTol)
        {
        // advancing
            thetaDp[pfacei] = HoffmanFunction(  mag(Ca[pfacei])
                                                + InvHoffFuncThetaAroot);
        }
        else if (uwall[pfacei] > uwallTol)
        {
        // receding
            thetaDp[pfacei] =HoffmanFunction( mag(Ca[pfacei]) + InvHoffFuncThetaRroot);
        }
        else
          {
        // pinned / equilibrium
            thetaDp[pfacei] = convertToRad*theta0;
        }
    }
    
    return convertToDeg*thetaDp;
}
//---------------------------------------------------------------------------------Define Hoffman--------------------------------------

scalar dynamicKistlerAlphaContactAngleFvPatchScalarField::HoffmanFunction
(
    const scalar& x
) const
{
    // Original Hoffman formula: 5.16 * (x / (1 + 1.31 * x^0.99))^0.706
    scalar inner = 5.16 * pow(x / (1.0 + 1.31 * pow(mag(x), 0.99)), 0.706);
    scalar tanh_val = tanh(inner);

    //cos(theta) = 1 - 2*tanh(f)
    scalar arg = 1.0 - 2.0 * tanh_val;

//    if (arg > 1.0) arg = 1.0;
//    if (arg < -1.0) arg = -1.0;

    return acos(arg);
}

void dynamicKistlerAlphaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    // Call base class write method to handle base class data
    dynamicAlphaContactAngleFvPatchScalarField::write(os);

    // Write scalar members explicitly
    os.writeKeyword("thetaA") << thetaA_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaR") << thetaR_ << token::END_STATEMENT << nl;

    // Write the dictionary keys for the referenced fields
    os.writeKeyword("muWater") << muWater_ << token::END_STATEMENT << nl;

    // Write the actual field value stored in this patch field object
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    dynamicKistlerAlphaContactAngleFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //


