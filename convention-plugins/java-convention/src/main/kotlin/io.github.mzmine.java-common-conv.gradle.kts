/*
 * Copyright (c) 2004-2024 The MZmine Development Team
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

import org.gradle.api.artifacts.VersionCatalogsExtension
import org.gradle.api.plugins.JavaPluginExtension
import org.gradle.api.tasks.JavaExec
import org.gradle.api.tasks.compile.JavaCompile
import org.gradle.api.tasks.testing.Test
import org.gradle.jvm.toolchain.JavaLanguageVersion
import org.gradle.kotlin.dsl.configure
import org.gradle.kotlin.dsl.getByType
import org.gradle.kotlin.dsl.named
import org.gradle.kotlin.dsl.withType

plugins {
    java
}

// default group that all modules will use
// in build.gradle use group = "${group}.subunit" to add another subunit to this group
group = "io.github.mzmine"

// https://github.com/gradle/gradle/issues/15383
val libs = extensions.getByType<VersionCatalogsExtension>().named("libs")

extensions.configure<JavaPluginExtension> {
    toolchain {
        languageVersion.set(JavaLanguageVersion.of(libs.findVersion("jdk").get().strictVersion))
    }
}

configurations.all {
    resolutionStrategy.eachDependency {
        if (requested.group == "org.apache.logging.log4j") {
            useVersion(libs.findVersion("log4j.core").get().preferredVersion)  // strict version
            because("patch all transitive dependencies to latest version")
        }

        /*if (requested.group == "org.apache.xmlgraphics") { // currently not needed, but may be requiredto enforce a single batik version in the future
            if (requested.name != "xmlgraphics-commons") { // different dependency version than batik itself
                val goodVersion = libs.findVersion("batik").get().requiredVersion
                useVersion(goodVersion) // strict version
                because("Replacing batik version " + requested.version + " with " + goodVersion)
            }
        }*/
    }
}

repositories {
    mavenCentral()
}

dependencies {
    // issue type safe access to libs is not implemented
    // https://github.com/gradle/gradle/issues/15383
    add("implementation", libs.findBundle("default-convention").get())
    // tests
    add("testImplementation", libs.findLibrary("junit.jupiter").get())
    add("testRuntimeOnly", libs.findLibrary("junit.platform").get())
    add("testImplementation", libs.findLibrary("mockito").get())
}

tasks.withType<JavaCompile>().configureEach {
    val compilerArgs = options.compilerArgs
    compilerArgs.add("--enable-preview")
    options.encoding = "UTF-8"
}

tasks.withType<JavaExec>().configureEach {
    jvmArgs("--enable-preview")
}

tasks.named<Test>("test") {
    useJUnitPlatform()
    testLogging {
        events("passed", "skipped", "failed")
    }
    jvmArgs("--enable-preview")
}
