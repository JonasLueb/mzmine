/*
 * Copyright (c) 2004-2025 The mzmine Development Team
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

package io.github.mzmine.main;

import io.github.mzmine.gui.DesktopService;
import io.github.mzmine.gui.preferences.MZminePreferences;
import io.github.mzmine.util.StringUtils;
import io.github.mzmine.util.files.FileAndPathUtil;
import io.mzio.mzmine.startup.MzmineCliArgs;
import io.mzio.users.gui.fx.LoginOptions;
import io.mzio.users.gui.fx.UsersController;
import io.mzio.users.user.CurrentUserService;
import java.io.File;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jetbrains.annotations.NotNull;

/**
 * Helper class to apply all parsed args from the parser-agnostic {@link MzmineCliArgs} holder to
 * the {@link ConfigService}.
 */
class ArgsToConfigUtils {

  private static final Logger logger = Logger.getLogger(ArgsToConfigUtils.class.getName());

  /**
   * Applies every relevant CLI argument to {@link ConfigService}. Reads the parser-agnostic
   * {@link MzmineCliArgs} produced by either the legacy Commons CLI parser or the picocli
   * subcommand tree.
   *
   * @param args parsed command-line arguments
   */
  static void applyArgsToConfig(final MzmineCliArgs args) {
    ConfigService.setTdfPseudoProfile(args.loadTdfPseudoProfile());

    checkAndLoadArgsConfiguration(args);
    TmpFileCleanup.runCleanup(); // clean old temp files in old dir

    // parse args temp dir after config was loaded, so we can override
    checkAndOverrideArgsTempDir(args);
    applyTempDirFromConfiguration();
    TmpFileCleanup.runCleanup(); // clean temp files in new dir

    checkAndOverrideArgsUser(args);

    checkAndOverrideArgsMemoryOption(args);

    setNumThreadsOverride(args);

    checkAndHandleArgsUserLoginOptions(args);

    ConfigService.setIgnoreParameterWarningsInBatch(args.ignoreParameterWarnings());
  }

  static void checkAndOverrideArgsTempDir(final @NotNull MzmineCliArgs args) {
    // override temp directory
    final File tempDirectory = args.tempDirectory();
    if (tempDirectory != null) {
      // needs to be accessible
      if (FileAndPathUtil.createDirectory(tempDirectory)) {
        ConfigService.getPreferences().setParameter(MZminePreferences.tempDirectory, tempDirectory);
      } else {
        logger.log(Level.WARNING,
            "Cannot create or access temp file directory that was set via program argument: "
                + tempDirectory.getAbsolutePath());
      }
    }
  }

  static void checkAndHandleArgsUserLoginOptions(final @NotNull MzmineCliArgs args) {
    final boolean isCliBatchProcessing = args.batchFile() != null;

    // login user by cli direct password
    if (args.cliLoginPassword()) {
      if (commandLineLogin(isCliBatchProcessing, LoginOptions.CONSOLE_ENTER_CREDENTIALS)) {
        return;
      }
    }

    // login user if cli option
    if (args.cliLogin()) {
      if (commandLineLogin(isCliBatchProcessing, LoginOptions.CONSOLE)) {
        return;
      }
    }
  }

  static void checkAndOverrideArgsMemoryOption(final @NotNull MzmineCliArgs args) {
    KeepInMemory keepInMemory;
    try {
      final String memory = args.keepInMemory();
      if (StringUtils.hasValue(memory)) {
        keepInMemory = KeepInMemory.parse(memory);

        // set to preferences
        ConfigService.getPreferences().setParameter(MZminePreferences.memoryOption, keepInMemory);
      } else {
        keepInMemory = ConfigService.getPreferences().getParameter(MZminePreferences.memoryOption)
            .getValue();
      }
    } catch (Exception exception) {
      logger.warning("Issue while reading keep in memory option from CLI argument");
      System.exit(1);
      return;
    }

    if (keepInMemory == null) {
      keepInMemory = KeepInMemory.NONE;
    }

    // apply memory management option
    keepInMemory.enforceToMemoryMapping();
  }

  static void checkAndOverrideArgsUser(final @NotNull MzmineCliArgs args) {
    if (args.userFile() == null) {
      // listen for user changes so that the latest user is saved
      final String username = ConfigService.getPreference(MZminePreferences.username);
      // this will set the current user to CurrentUserService
      // loads all users already logged in from the user folder
      if (StringUtils.hasValue(username)) {
        UsersController.getInstance().setCurrentUserByName(username);
      }
    }
  }

  static void checkAndLoadArgsConfiguration(final @NotNull MzmineCliArgs args) {
    // override preferences file by command line argument pref
    final File prefFile = Objects.requireNonNullElse(args.preferencesFile(),
        MZmineConfiguration.CONFIG_FILE);
    if ("null".equals(prefFile.getName())) {
      logger.info("Preference file was set to null, not loading configuration.");
      return;
    }

    // Load configuration
    if (prefFile.exists() && prefFile.canRead()) {
      try {
        ConfigService.getConfiguration().loadConfiguration(prefFile, true);
      } catch (Exception e) {
        logger.log(Level.WARNING, "Error while reading configuration " + prefFile.getAbsolutePath(),
            e);
      }
    } else {
      logger.log(Level.WARNING, "Cannot read configuration " + prefFile.getAbsolutePath());
    }
  }

  /**
   * Set number of cores to automatic or to fixed number
   */
  static void setNumThreadsOverride(final @NotNull MzmineCliArgs args) {
    final String numCores = args.numCores();
    if (numCores != null) {
      // set to preferences
      var parameter = ConfigService.getPreferences().getParameter(MZminePreferences.numOfThreads);
      if (numCores.equalsIgnoreCase("auto") || numCores.equalsIgnoreCase("automatic")) {
        parameter.setAutomatic(true);
      } else {
        try {
          parameter.setValue(Integer.parseInt(numCores));
        } catch (Exception ex) {
          logger.log(Level.SEVERE,
              "Cannot parse command line argument threads (int) set to " + numCores);
          throw new IllegalArgumentException("numCores was set to " + numCores, ex);
        }
      }
    }
  }

  /**
   * @param isCliBatchProcessing
   * @param option
   * @return true if application finished
   */
  static boolean commandLineLogin(final boolean isCliBatchProcessing, LoginOptions option) {
    boolean success = false;
    try {
      logger.info("CLI user login");
      UsersController.getInstance().loginOrRegisterConsoleBlocking(option);
      success = true;
    } catch (Exception ex) {
      DesktopService.getDesktop().displayMessage(
          "Requires user login. Open mzmine GUI and login to a user. Then provide the user file as command line argument -user path/user.mzuser");
      if (!isCliBatchProcessing) {
        System.exit(1);
        return true;
      }
    }
    // if no batch select - that means it was only a login call.
    // save config and close mzmine
    if (success && !isCliBatchProcessing) {
      String currentUserName = CurrentUserService.getUserName().orElse("");
      ConfigService.getPreferences().setParameter(MZminePreferences.username, currentUserName);
      if (!ConfigService.saveUserConfig()) {
        logger.severe(
            "Failed to save user config after login. A solution may be to delete the .mzconfig file in the system user directory /.mzmine/");
        System.exit(1);
        return true;
      } else {
        logger.info("User login successful, user configuration is saved with the new user "
            + currentUserName);
        System.exit(0);
        return true;
      }
    }
    return false;
  }

  static void applyTempDirFromConfiguration() {
    final File tempDir = ConfigService.getConfiguration().getPreferences()
        .getParameter(MZminePreferences.tempDirectory).getValue();
    if (tempDir == null) {
      logger.warning(
          () -> "Invalid temporary directory. Defaulting to system temp directory. %s".formatted(
              System.getProperty("java.io.tmpdir")));
      return;
    }

    if (!tempDir.exists()) {
      if (!tempDir.mkdirs()) {
        logger.warning(
            () -> "Could not create temporary directory %s. Defaulting to system temp directory. %s".formatted(
                tempDir.getAbsolutePath(), System.getProperty("java.io.tmpdir")));
        return;
      }
    }

    if (tempDir.isDirectory()) {
      FileAndPathUtil.setTempDir(tempDir.getAbsoluteFile());
      logger.finest(() -> "Default temporary directory is %s".formatted(
          System.getProperty("java.io.tmpdir")));
      logger.finest(() -> "Working temporary directory is %s".formatted(
          FileAndPathUtil.getTempDir().toString()));
    }
  }
}
