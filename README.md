# KistlerHoffman OpenFOAM

These are two codependate pieces of code. It's been updated from V7 to v2412 of OpenFOAM with some improvements made to fix an issue with the receding phase of the droplet. 

### Contact Line install
For the contact line model simple place it in your OpenFOAM users src/ section, or where ever you deem most suits you. Then use the Allmake.sh file to install it.

### interFOAM update
This is a copy of interFOAM, the only change made was to make the solver write mu and sigma so that the contact angle model could access them. 
Place it in your OpenFOAM user folders, I have mine in a subfolder called applications. Then use the wmake command. 
