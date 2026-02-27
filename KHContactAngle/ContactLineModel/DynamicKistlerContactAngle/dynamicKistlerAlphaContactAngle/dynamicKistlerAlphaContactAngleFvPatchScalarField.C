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
#include "Time.H" // Required for this->db() to return Time object, which is the objectRegistry
#include "interfaceProperties.H"
#include "twoPhaseMixture.H" // Often needed to find phase properties
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

const scalar dynamicKistlerAlphaContactAngleFvPatchScalarField::theta0 = 90.0;

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
    muWater_(0.0),
    sigmaName_("undefined")
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
    muWater_(acpsf.muWater_),
    sigmaName_(acpsf.sigmaName_)
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
    muWater_(readScalar(dict.lookup("muWater"))),
    sigmaName_(dict.lookup("sigmaKistler"))
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
    muWater_(acpsf.muWater_),
    sigmaName_(acpsf.sigmaName_)
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
    muWater_(acpsf.muWater_),
    sigmaName_(acpsf.sigmaName_)
{}


//----------------------------------------------------------- * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> dynamicKistlerAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
    //eb - Lookup and return the patchField of dynamic viscosity of mixture
    //     and surface tension
    if((sigmaName_ != "sigmaKistler"))
    {
        FatalErrorIn
        (
            "dynamicKistlerAlphaContactAngleFvPatchScalarField"
        )   << " muEffKistler or sigma set inconsitently"
            << "sigmaKistler = " << sigmaName_ << '.' << nl
            << "     Set both muEffKistler and sigmaKistler according to the "
            << "definition of dynamicKistlerAlphaContactAngle"
            << "\n     on patch " << this->patch().name()
            << exit(FatalError);
    }


    // * * * v7 to v2412 MIGRATION CHANGE: Field Lookup Update * * *
    // Old: patch().lookupPatchField<volScalarField, scalar>(muName_);
    // New: Use the object registry (this->db()) to look up the volField
    // and extract the corresponding patchField.

    const objectRegistry& obr = this->db();

    const fvPatchField<scalar>& sigmap =
        obr.lookupObject<volScalarField>(sigmaName_).boundaryField()[this->patch().index()];

//----------------------------------------------------------------------------------------------------Vectors and Fields ----------------------------------------------------

    vectorField nf = patch().nf();

    // Calculate the component of the velocity parallel to the wall
    vectorField Uwall = Up.patchInternalField(); //- Up; Up.patchInternalField() is the U in the cells adjacent to the bounday.
    Uwall -= (nf & Uwall)*nf;

    // Find the direction of the interface parallel to the wall
    vectorField nWall = nHat - (nf & nHat)*nf;

    // Normalise nWall
    nWall /= (mag(nWall) + SMALL); // SMALL = 1.0e-6; Points into the drop from the interface

    // Calculate Uwall resolved normal to the interface perpendicular to
    // the wall (substrate)
    scalarField uwall = nWall & Uwall; // & is the inner product...... faces into the droplet when receding and outwards when advancing

//----------------------------------------------------------------------------------------------------Average uwall over the triple-Line------------------------------------

//  Initialize sums
    scalar sumUwallArea = 0.0;
    scalar sumArea = 0.0;

    // Access area and alpha fields
    const scalarField& magSf = patch().magSf();
    const scalarField& alpha = *this;

    //Define threshold for "interface" (e.g. 0.01 to 0.99)
    const scalar alphaMin = 0.001;
    const scalar alphaMax = 0.999;

    //Accumulate weighted sums over local faces
    forAll(uwall, pfacei)
    {
        if (alpha[pfacei] > alphaMin && alpha[pfacei] < alphaMax)
        {
            sumUwallArea += uwall[pfacei] * magSf[pfacei];
            sumArea += magSf[pfacei];
        }
    }

    //Parallel reduction (sums values across all processors)
    reduce(sumUwallArea, sumOp<scalar>());
    reduce(sumArea, sumOp<scalar>());

    //Calculate average
    scalar uwallAvg = 0.0;
    if (sumArea > VSMALL)
    {
        uwallAvg = sumUwallArea / sumArea;
        forAll(uwall, pfacei)
        {
            if (alpha[pfacei] > alphaMin && alpha[pfacei] < alphaMax)
            {
                uwall[pfacei] = uwallAvg;
            }
        }
    }
//---------------------------------------------------------------------------------------------------- End uwall Averaging Code
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
    scalar InvHoffFuncThetaAroot = RRInvHoffFuncThetaA.root(0,65);

    RiddersRoot RRInvHoffFuncThetaR(InvHoffFuncThetaR, 1.e-10);
    scalar InvHoffFuncThetaRroot = RRInvHoffFuncThetaR.root(0,65);

    //eb - Calculate and return the value of contact angle on patch faces,
    //     a general approach: the product of Uwall and nWall is negative
    //     for advancing and positiv for receding motion.
    //     thetaDp is initialized to 90 degrees corresponding to no wall
    //     adhesion

//-------------------------------------------------------------------------------Kistler------------------------------------------------------------------------
const scalar uwallTol = 5e-3; // or velocity-based
//const scalar uwallTolR = 1e-6; // Receding condition more agressive to deal with spurious currents along the interface
//
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
            // Symmetric Kistler: 180 - f( Ca + f_inv(180 - thetaR) )
            // As Ca increases, f(...) goes to 180, so thetaDp goes to 0.
            thetaDp[pfacei] = HoffmanFunction( mag(Ca[pfacei]) + InvHoffFuncThetaRroot);
//            thetaDp[pfacei] = constant::mathematical::pi - thetaGasDyn;
        }
        else
          {
        // pinned / equilibrium
            thetaDp[pfacei] = convertToRad*theta0;
        }
   }

//- -- ------------------------------------------------ DATA EXPORT FOR PYTHON ----------------------

    //const scalarField& alphaWall = *this;
//
  //  if (this->db().time().writeTime())
//    {
        // Define path: postProcessing/contactAngleData/<patchName>_time.csv
    //    fileName outputDir = this->db().time().path()/"postProcessing"/"contactAngleData";
  //      mkDir(outputDir);
//
      //  fileName logName = outputDir/(this->patch().name() + "_" + this->db().time().timeName() + ".csv");
//
	// Declare the pointer locally
    //    autoPtr<OFstream> OF_logFilePtr;
  //      OF_logFilePtr.reset(new OFstream(logName));
//
        // Write CSV Header
  //      *OF_logFilePtr << "faceID,uwall,Ca,thetaDeg,alpha" << endl;
//
    //    forAll(uwall, facei)
  //      {
//            *OF_logFilePtr << facei << ","
            //               << uwall[facei] << ","
          //                 << Ca[facei] << ","
        //                   << (thetaDp[facei] * convertToDeg) << ","
      //                     << alphaWall[facei] << endl;
    //    }
  //  }
//
    return convertToDeg*thetaDp;
}
//------------------------------------------------------------------------------------Kistler and Data Export End-------------------------


//---------------------------------------------------------------------------------Define Hoffman--------------------------------------
scalar dynamicKistlerAlphaContactAngleFvPatchScalarField::HoffmanFunction
(
    const scalar& x
) const
{
    // Use mag(x) to prevent domain errors with pow(x, 0.99) if x is slightly negative
    //if (x < 0 ) x = 0;

    // Original Hoffman formula: 5.16 * (x / (1 + 1.31 * x^0.99))^0.706
    scalar inner = 5.16 * pow(x / (1.0 + 1.31 * pow(mag(x), 0.99)), 0.706);
    scalar tanh_val = tanh(inner);

    // Kistler model maps the contact angle as: cos(theta) = 1 - 2*tanh(f)
    scalar arg = 1.0 - 2.0 * tanh_val;

    // Explicitly clamp for acos stability to prevent NaNs from precision drift
    if (arg > 1.0) arg = 1.0;
    if (arg < -1.0) arg = -1.0;

    return acos(arg);
}

// * * * v7 to v2412 MIGRATION CHANGE: Full write function restored * * *
// The original provided code had an incomplete write function active.
// This version is complete and writes all custom parameters.

void dynamicKistlerAlphaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    // Call base class write method to handle base class data
    dynamicAlphaContactAngleFvPatchScalarField::write(os);

    // Write scalar members explicitly
    os.writeKeyword("thetaA") << thetaA_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaR") << thetaR_ << token::END_STATEMENT << nl;

    // Write the dictionary keys for the referenced fields
    os.writeKeyword("muWater") << muWater_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaKistler") << sigmaName_ << token::END_STATEMENT << nl;

    // Write the actual field value stored in this patch field object
    writeEntry("value", os);
}

// * * * END MIGRATION CHANGE * * *


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    dynamicKistlerAlphaContactAngleFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //


