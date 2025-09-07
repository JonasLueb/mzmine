#!/bin/bash
CP=$(./gradlew -q :mzmine-community:printRuntimeClasspath)
CP="mzmine-community/build/classes/java/main:mzmine-community/build/resources/main:$CP"
export CP