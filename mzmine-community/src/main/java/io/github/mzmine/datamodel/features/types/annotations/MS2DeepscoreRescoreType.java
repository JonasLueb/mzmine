package io.github.mzmine.datamodel.features.types.annotations;

import io.github.mzmine.datamodel.features.types.numbers.abstr.ScoreType;
import org.jetbrains.annotations.NotNull;

/**
 * Column type exposing the heuristically rescored MS2Deepscore similarity.
 */
public class MS2DeepscoreRescoreType extends ScoreType {

  @NotNull
  @Override
  public String getUniqueID() {
    return "ms2deepscore_rescored_similarity";
  }

  @NotNull
  @Override
  public String getHeaderString() {
    return "MS2Deepscore rescored";
  }

  @Override
  public boolean getDefaultVisibility() {
    return true;
  }

  @Override
  public int getPrefColumnWidth() {
    return 120;
  }
}

