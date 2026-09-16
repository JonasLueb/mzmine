/*
 * Copyright (c) 2004-2026 The mzmine Development Team
 */
package io.github.mzmine.modules.batchmode;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import io.github.mzmine.datamodel.MZmineProject;
import io.github.mzmine.taskcontrol.TaskStatus;
import java.time.Instant;
import org.junit.jupiter.api.Test;
import org.mockito.Mockito;

class BatchTaskFixedProjectTest {

  @Test
  void fixedProjectTaskCancelsBeforeStartingAgainstAnotherProject() {
    final MZmineProject reviewedProject = Mockito.mock(MZmineProject.class);
    final BatchModeParameters parameters = new BatchModeParameters();
    parameters.getParameter(BatchModeParameters.batchQueue).setValue(new BatchQueue());
    final BatchTask task = BatchTask.forFixedProject(reviewedProject, parameters, Instant.now());

    task.run();

    assertEquals(TaskStatus.CANCELED, task.getStatus());
    assertEquals("The reviewed project changed before the batch could continue.",
        task.getErrorMessage());
    assertTrue(task.getResultFeatureLists().isEmpty());
  }
}
