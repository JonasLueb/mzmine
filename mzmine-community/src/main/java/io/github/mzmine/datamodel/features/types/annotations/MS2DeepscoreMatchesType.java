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

package io.github.mzmine.datamodel.features.types.annotations;

import io.github.mzmine.datamodel.features.types.DataType;
import io.github.mzmine.datamodel.features.types.JsonStringType;
import io.github.mzmine.datamodel.features.types.ListWithSubsType;
import io.github.mzmine.datamodel.features.types.annotations.CompoundNameType;
import io.github.mzmine.datamodel.features.types.annotations.InChIStructureType;
import io.github.mzmine.datamodel.features.types.annotations.MolecularStructureType;
import io.github.mzmine.datamodel.features.types.annotations.MS2DeepscoreRescoreType;
import io.github.mzmine.datamodel.features.types.annotations.SmilesStructureType;
import io.github.mzmine.datamodel.features.types.annotations.formula.FormulaType;
import io.github.mzmine.datamodel.features.types.annotations.iin.IonAdductType;
import io.github.mzmine.datamodel.features.types.modifiers.AnnotationType;
import io.github.mzmine.datamodel.features.types.numbers.PrecursorMZType;
import io.github.mzmine.datamodel.features.types.numbers.MatchingSignalsType;
import io.github.mzmine.datamodel.features.types.numbers.scores.ExplainedIntensityPercentType;
import io.github.mzmine.datamodel.features.types.numbers.scores.SimilarityType;
import io.github.mzmine.modules.dataprocessing.id_ms2deepscore_vectorsearch.MS2DeepscoreScoringUtils;
import io.github.mzmine.util.spectraldb.entry.DBEntryField;
import io.github.mzmine.util.spectraldb.entry.SpectralDBAnnotation;
import io.github.mzmine.util.spectraldb.entry.SpectralLibraryEntry;
import java.util.List;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * Separate column group for MS2Deepscore ANN matches to allow side-by-side comparison with
 * classical spectral library search results.
 */
public class MS2DeepscoreMatchesType extends ListWithSubsType<SpectralDBAnnotation> implements
    AnnotationType {

  private static final List<DataType> subTypes = List.of(
      new MS2DeepscoreMatchesType(),
      new CompoundNameType(),
      new SimilarityType(),
      new MS2DeepscoreRescoreType(),
      new MatchingSignalsType(),
      new ExplainedIntensityPercentType(),
      new IonAdductType(),
      new FormulaType(),
      new MolecularStructureType(),
      new SmilesStructureType(),
      new InChIStructureType(),
      new PrecursorMZType(),
      new JsonStringType()
  );

  @NotNull
  @Override
  public String getUniqueID() {
    return "ms2deepscore_matches";
  }

  @NotNull
  @Override
  public List<DataType> getSubDataTypes() {
    return subTypes;
  }

  @Override
  public <K> @Nullable K map(@NotNull final DataType<K> subType, final SpectralDBAnnotation match) {
    final SpectralLibraryEntry entry = match.getEntry();
    return (K) switch (subType) {
      case MS2DeepscoreMatchesType __ -> match;
      case CompoundNameType __ -> {
        final String name = entry.getField(DBEntryField.NAME).orElse("").toString();
        if (name != null && !name.isBlank()) yield name;
        // fallback to DB# if no name present
        yield entry.getField(DBEntryField.ENTRY_ID).orElse("").toString();
      }
      case FormulaType __ -> entry.getField(DBEntryField.FORMULA).orElse("").toString();
      case IonAdductType __ -> entry.getField(DBEntryField.ION_TYPE).orElse("").toString();
      case MolecularStructureType __ -> entry.getStructure();
      case SmilesStructureType __ -> entry.getField(DBEntryField.SMILES).orElse("").toString();
      case InChIStructureType __ -> entry.getField(DBEntryField.INCHI).orElse("").toString();
      case SimilarityType __ -> (float) match.getSimilarity().getScore();
      case MS2DeepscoreRescoreType __ -> (float) MS2DeepscoreScoringUtils.computeRescore(match);
      case MatchingSignalsType __ -> match.getSimilarity().getOverlap();
      case ExplainedIntensityPercentType __ -> match.getSimilarity().getExplainedLibraryIntensity();
      case PrecursorMZType __ -> entry.getField(DBEntryField.PRECURSOR_MZ).orElse(null);
      case JsonStringType __ -> entry.getOrElse(DBEntryField.JSON_STRING, null);
      default -> throw new UnsupportedOperationException(
          "DataType %s is not covered in map".formatted(subType.toString()));
    };
  }

  @NotNull
  @Override
  public String getHeaderString() {
    return "MS2Deepscore match";
  }

  @Override
  public boolean getDefaultVisibility() {
    return true;
  }

  @Override
  public int getPrefColumnWidth() {
    return 150;
  }
}
