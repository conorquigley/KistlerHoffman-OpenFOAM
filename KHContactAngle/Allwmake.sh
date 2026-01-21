#!/bin/bash
echo "Cleaning old build files"
wclean libso DynamicKistlerContactAngle
echo "Building library"
wmake libso DynamicKistlerContactAngle

