plugins {
    id("io.github.mzmine.java-common-conv")
}

pluginManager.apply("java-library")
pluginManager.apply("maven-publish")

repositories {
    mavenLocal()
    mavenCentral()
}
