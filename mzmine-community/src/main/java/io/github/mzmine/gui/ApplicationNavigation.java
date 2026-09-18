/*
 * Copyright (c) 2004-2026 The mzmine Development Team
 */
package io.github.mzmine.gui;

import io.github.mzmine.datamodel.features.FeatureList;
import io.github.mzmine.gui.mainwindow.MZmineTab;
import io.github.mzmine.main.MZmineCore;
import io.github.mzmine.modules.tools.batchwizard.BatchWizardTab;
import io.github.mzmine.modules.visualization.featurelisttable_modular.FeatureTableTab;
import io.github.mzmine.util.FeatureTableFXUtil;
import javafx.application.Platform;
import javafx.stage.Stage;
import javafx.stage.Window;
import org.jetbrains.annotations.NotNull;

/** Common, local UI destinations for integrations. All methods must run on the JavaFX thread. */
public final class ApplicationNavigation {

  private ApplicationNavigation() {
  }

  /** Opens the existing processing wizard, or creates one, and brings its window to the front. */
  public static @NotNull BatchWizardTab openOrFocusWizard() {
    requireFxThread();
    final BatchWizardTab wizard = MZmineCore.getDesktop().getAllTabs().stream()
        .filter(BatchWizardTab.class::isInstance).map(BatchWizardTab.class::cast).findFirst()
        .orElseGet(() -> {
          final BatchWizardTab newWizard = new BatchWizardTab();
          MZmineCore.getDesktop().addTab(newWizard);
          return newWizard;
        });
    focus(wizard);
    return wizard;
  }

  /** Focuses a particular wizard, for example the wizard that owns a reviewed proposal. */
  public static void focusWizard(final @NotNull BatchWizardTab wizard) {
    requireFxThread();
    focus(wizard);
  }

  /** Opens an existing table for the feature list, or creates one, and brings it to the front. */
  public static void showFeatureTable(final @NotNull FeatureList featureList) {
    requireFxThread();
    final MZmineTab existing = MZmineCore.getDesktop().getAllTabs().stream()
        .filter(FeatureTableTab.class::isInstance).map(FeatureTableTab.class::cast)
        .filter(tab -> tab.getFeatureList() == featureList).findFirst().orElse(null);
    if (existing == null) {
      FeatureTableFXUtil.addFeatureTableTab(featureList);
      focusMainWindow();
    } else {
      focus(existing);
    }
  }

  public static void focus(final @NotNull MZmineTab tab) {
    requireFxThread();
    if (tab.getTabPane() != null) {
      tab.getTabPane().getSelectionModel().select(tab);
      final Window window = tab.getTabPane().getScene() == null ? null
          : tab.getTabPane().getScene().getWindow();
      if (window != null) {
        if (window instanceof Stage stage) {
          stage.setIconified(false);
          stage.toFront();
        }
        window.requestFocus();
        return;
      }
    }
    focusMainWindow();
  }

  private static void focusMainWindow() {
    final Stage mainWindow = MZmineCore.getDesktop().getMainWindow();
    mainWindow.setIconified(false);
    mainWindow.toFront();
    mainWindow.requestFocus();
  }

  private static void requireFxThread() {
    if (!Platform.isFxApplicationThread()) {
      throw new IllegalStateException(
          "ApplicationNavigation must be called on the JavaFX application thread; use FxThread.runLater");
    }
  }
}
