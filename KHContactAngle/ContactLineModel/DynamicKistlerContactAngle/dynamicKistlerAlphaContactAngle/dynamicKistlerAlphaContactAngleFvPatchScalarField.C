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

const scalar dynamicKistlerAlphaContactAngleFvPatchScalarField::theta0 = 110.0;

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
    muName_("undefined"),
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
    muName_(acpsf.muName_),
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
    muName_(dict.lookup("muEffKistler")),
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
    muName_(acpsf.muName_),
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
    muName_(acpsf.muName_),
    sigmaName_(acpsf.sigmaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> dynamicKistlerAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
    //eb - Lookup and return the patchField of dynamic viscosity of mixture
    //     and surface tension
    if((muName_ != "muEffKistler") || (sigmaName_ != "sigmaKistler"))
    {
        FatalErrorIn
        (
            "dynamicKistlerAlphaContactAngleFvPatchScalarField"
        )   << " muEffKistler or sigma set inconsitently, muEffKistler = "
            << muName_ << ", sigmaKistler = " << sigmaName_ << '.' << nl
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

    const fvPatchField<scalar>& mup =
        obr.lookupObject<volScalarField>(muName_).boundaryField()[this->patch().index()];

    const fvPatchField<scalar>& sigmap =
        obr.lookupObject<volScalarField>(sigmaName_).boundaryField()[this->patch().index()];
    // * * * END MIGRATION CHANGE * * *

    vectorField nf = patch().nf();

    // Calculate the component of the velocity parallel to the wall
    vectorField Uwall = Up.patchInternalField(); //- Up; Up.patchInternalField() is the U in the cells adjacent to the bounday.
    Uwall -= (nf & Uwall)*nf;

    // Find the direction of the interface parallel to the wall
    vectorField nWall = nHat - (nf & nHat)*nf; // can this be used to define the receding or advancing contact angle?????????????????????????????????

    // Normalise nWall
    nWall /= (mag(nWall) + SMALL); // SMALL = 1.0e-6; Points into the drop from the interface

    // Calculate Uwall resolved normal to the interface perpendicular to
    // the wall
    scalarField uwall = nWall & Uwall; // & is the inner product......................... faces into the droplet when receding and outwards when advancing

    //eb - Calculate local Capillary number
    scalarField Ca = mup*mag(uwall)/sigmap;

    //eb - Instantiate function object InverseHoffmanFunction for thetaA and thetaR
    dynamicKistlerAlphaContactAngleFvPatchScalarField::InverseHoffmanFunction
    InvHoffFuncThetaA
    (
        convertToRad*thetaA_
    );

    dynamicKistlerAlphaContactAngleFvPatchScalarField::InverseHoffmanFunction
    InvHoffFuncThetaR
    (
        convertToRad*(thetaR_)//calculate the receding contact angle
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

//Original
   // scalarField thetaDp(patch().size(), convertToRad*theta0);
  //  forAll(uwall, pfacei)
//    {
//
        //if(uwall[pfacei] < 0.0)
        //{
      //      thetaDp[pfacei] = HoffmanFunction(  mag(Ca[pfacei])
    //                                            + InvHoffFuncThetaAroot);
  //      }
//	else if (uwall[pfacei] >  0.0)
      //  {
    //        thetaDp[pfacei] = HoffmanFunction(   mag(Ca[pfacei])
  //                                              + InvHoffFuncThetaRroot);
//		//attempted improvement
    //        //scalar thetaGasDyn = HoffmanFunction(Ca[pfacei] + InvHoffFuncThetaRroot);
  //          //thetaDp[pfacei] = constant::mathematical::pi - thetaGasDyn;
//
//        }
//   }
//
//    return convertToDeg*thetaDp;
//}

//improved
const scalar uwallTol = 1e-4; // or velocity-based
//
    scalarField thetaDp(patch().size(), convertToRad*theta0);
    forAll(uwall, pfacei)
    {

        if (uwall[pfacei] < -uwallTol)
        {
        // advancing
            thetaDp[pfacei] = HoffmanFunction(  mag(Ca[pfacei])
                                                + 0.1077076);
        }
        else if (uwall[pfacei] > uwallTol)
        {
        // receding
            // Symmetric Kistler: 180 - f( Ca + f_inv(180 - thetaR) )
            // As Ca increases, f(...) goes to 180, so thetaDp goes to 0.
            thetaDp[pfacei] = HoffmanFunction( mag(Ca[pfacei]) + 0.01797121);
//            thetaDp[pfacei] = constant::mathematical::pi - thetaGasDyn;
        }
        else
          {
        // pinned / equilibrium
            thetaDp[pfacei] = convertToRad*theta0;
        }
   }

    // --  DATA EXPORT FOR PYTHON ---.

	const scalarField& alphaWall = *this;

    if (this->db().time().writeTime())
    {
        // Define path: postProcessing/contactAngleData/<patchName>_time.csv
        fileName outputDir = this->db().time().path()/"postProcessing"/"contactAngleData";
        mkDir(outputDir);

        fileName logName = outputDir/(this->patch().name() + "_" + this->db().time().timeName() + ".csv");

	// Declare the pointer locally
        autoPtr<OFstream> OF_logFilePtr;
        OF_logFilePtr.reset(new OFstream(logName));

        // Write CSV Header
        *OF_logFilePtr << "faceID,uwall,Ca,thetaDeg,alpha" << endl;

        forAll(uwall, facei)
        {
            *OF_logFilePtr << facei << ","
                           << uwall[facei] << ","
                           << Ca[facei] << ","
                           << (thetaDp[facei] * convertToDeg) << ","
                           << alphaWall[facei] << endl;
        }
    }
    // --- END DATA EXPORT ---

    return convertToDeg*thetaDp;
}

//scalar dynamicKistlerAlphaContactAngleFvPatchScalarField::HoffmanFunction
//(
  //  const scalar& x
//) const
//{
//    return acos(1 - 2*tanh(5.16*pow(x/(1+1.31*pow(mag(x),0.99)),0.706)));
//}

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
    os.writeKeyword("muEffKistler") << muName_ << token::END_STATEMENT << nl;
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


