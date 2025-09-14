#!/bin/bash
CP=$(./gradlew -q :mzmine-community:printRuntimeClasspath)
CP="mzmine-community/build/classes/java/main:mzmine-community/build/resources/main:$CP"
export CP

./gradlew :mzmine-community:classes -x test -PskipLicenseTasks=true --daemon --offline

java --enable-preview -cp "$CP" -DinputDir=/Users/A14D579/Downloads/trainingsdaten_restriktiv/038_Sa13_Methanol_POS.raw -DoutDir=/Users/A14D579/Downloads/trainingsdaten_restriktiv/038_Sa13_Methanol_POS.raw/resolver_hparam -Dformat=html -Diou=0.3 -DrtTol=0.02 -DminPoints=3 -DpreZeroQ=1 -DbaselineFrac=0 io.github.mzmine.tools.eicresolver.EICGridRunner