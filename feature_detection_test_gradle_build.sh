#!/bin/bash
# Fast development run script for MZmine
# add "clean" after gradlew command for a clean build
echo "Starting EIC Feature Resolver test..."
./gradlew :mzmine-community:runEICResolver --configuration-cache -PskipLicenseTasks=true \
  -x :mzmine-community:checkLicensePreparation \
  -x :mzmine-community:checkLicense \
  -x :mzmine-community:generateLicenseReport \
  -x :mzmine-community:copyLicenseInformationToResources \
  -x :mzmine-community:copyLicenseInformationToBuild \
  -DinputDir=/Users/A14D579/Downloads/trainingsdaten_restriktiv/038_Sa13_Methanol_POS.raw \
  -DoutDir=/Users/A14D579/Downloads/trainingsdaten_restriktiv/038_Sa13_Methanol_POS.raw/resolver_eval_html5 \
  -Dformat=html -Diou=0.3 -DrtTol=0.02 -DpreZeroQ=1 -DbaselineFrac=0 -DminAbs=0 -DsgMin=0 -DmsThresh=0 \
  -DcwtSnr=1 -DcwtMin=0 -DcwtPersist=2