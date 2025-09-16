/*
 * Copyright (c) 2004-2025 The mzmine Development Team
 */

package io.github.mzmine.modules.dataprocessing.id_ms2deepscore_vectorsearch;

import io.github.mzmine.datamodel.MZmineProject;
import io.github.mzmine.datamodel.features.FeatureList;
import io.github.mzmine.datamodel.features.ModularFeatureList;
import io.github.mzmine.parameters.parametertypes.selectors.FeatureListsSelection;
import io.github.mzmine.project.impl.MZmineProjectImpl;
import io.github.mzmine.util.ExitCode;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Collection;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

class MS2DeepscoreVectorSearchModuleTest {

  @Test
  void testRunModuleCreatesTaskAndReturnsOK() {
    MZmineProject project = new MZmineProjectImpl();
    ModularFeatureList fl = FeatureList.createDummy();
    project.addFeatureList(fl);

    MS2DeepscoreVectorSearchParameters params = new MS2DeepscoreVectorSearchParameters();
    FeatureListsSelection sel = new FeatureListsSelection(fl);
    params.getParameter(MS2DeepscoreVectorSearchParameters.PEAK_LISTS).setValue(sel);

    MS2DeepscoreVectorSearchModule module = new MS2DeepscoreVectorSearchModule();
    Collection<io.github.mzmine.taskcontrol.Task> tasks = new ArrayList<>();
    ExitCode code = module.runModule(project, params, tasks, Instant.now());

    Assertions.assertEquals(ExitCode.OK, code);
    Assertions.assertFalse(tasks.isEmpty());
    Assertions.assertTrue(tasks.iterator().next() instanceof MS2DeepscoreVectorSearchTask);
  }
}


