package io.github.mzmine.datamodel;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/** Canonical acquisition terms plus native-only declared fields. Never serialize this object to an agent. */
public record AcquisitionMetadata(@NotNull List<Term> terms, @NotNull Map<String, String> localFields) {
  public static final AcquisitionMetadata EMPTY = new AcquisitionMetadata(List.of());
  public enum Field { INSTRUMENT_MODEL, ANALYZER, IONIZATION, DETECTOR, ACQUISITION_METHOD }
  public record Term(@NotNull Field field, @NotNull String accession, @NotNull String label) {}
  private static final Map<String, List<Term>> VOCABULARY = loadVocabulary();

  public AcquisitionMetadata {
    terms = terms.stream().distinct().toList();
    localFields = Map.copyOf(localFields);
  }

  public AcquisitionMetadata(final @NotNull List<Term> terms) { this(terms, Map.of()); }

  /** Exact canonical label match only; unknown vendor model strings remain local. */
  public static @NotNull List<Term> resolveLabel(final @NotNull Field field, final @Nullable String label) {
    if (label == null || label.isBlank()) return List.of();
    return VOCABULARY.values().stream().flatMap(List::stream).filter(term -> term.field() == field
        && term.label().equalsIgnoreCase(label.strip())).toList();
  }

  /** Recognized controlled terms; arbitrary CV names and values are deliberately ignored. */
  public static @NotNull List<Term> resolve(final @Nullable String accession) {
    return accession == null ? List.of() : VOCABULARY.getOrDefault(accession, List.of());
  }

  public @NotNull AcquisitionMetadata plus(final @NotNull AcquisitionMetadata other) {
    final List<Term> merged = new ArrayList<>(terms);
    merged.addAll(other.terms);
    final Map<String, String> fields = new java.util.LinkedHashMap<>(localFields);
    fields.putAll(other.localFields);
    return new AcquisitionMetadata(merged, fields);
  }

  private static @NotNull Map<String, List<Term>> loadVocabulary() {
    try (final var reader = new BufferedReader(new InputStreamReader(Objects.requireNonNull(
        AcquisitionMetadata.class.getResourceAsStream("/acquisition-cv.tsv")), StandardCharsets.UTF_8))) {
      return reader.lines().filter(line -> !line.startsWith("#") && !line.isBlank())
          .map(line -> line.split("\t", 3))
          .map(parts -> new Term(Field.valueOf(parts[0]), parts[1], parts[2]))
          .collect(Collectors.groupingBy(Term::accession, Collectors.toUnmodifiableList()));
    } catch (java.io.IOException exception) {
      throw new IllegalStateException("Cannot read packaged acquisition vocabulary", exception);
    }
  }
}
