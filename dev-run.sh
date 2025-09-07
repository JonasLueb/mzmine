#!/bin/bash
# Fast development run script for MZmine
# add "clean" after gradlew command for a clean build
echo "Starting MZmine in development mode (fast)..."
export JAVA_TOOL_OPTIONS="$JAVA_TOOL_OPTIONS --add-exports=javafx.base/com.sun.javafx.event=ALL-UNNAMED"
./gradlew :mzmine-community:run -x test -x checkLicense -x generateLicenseReport -x copyLicenseInformationToResources -x copyLicenseInformationToBuild --offline --daemon