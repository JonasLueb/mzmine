# Java formatting

Run formatting through the Gradle wrapper. Spotless **8.10.4** downloads the Eclipse JDT
**4.41** formatter library automatically. No Eclipse or IntelliJ installation, IDE plugin,
or paid license is required. Use the project's Java 26 JDK.

```sh
./gradlew spotlessApply
./gradlew spotlessCheck
```

Both commands cover every Java subproject in this build. To format only one edited file:

```sh
./gradlew :mzmine-community:spotlessApply -PspotlessIdeHook=/absolute/path/File.java
```

Replace `mzmine-community` with the module owning the file. The selected file must still be
within the Java targets and changed from the adoption baseline. Each module's normal
`check` task also depends on `spotlessCheck`; checking never rewrites source.
The dedicated Java formatting workflow runs the root check on pull requests and master pushes.

## Profile and IDE compatibility

[gradle/spotless/eclipse.xml](../gradle/spotless/eclipse.xml) is the formatting authority
for Gradle. It is identical in the community and Pro repositories, with a copy in each
so Pro also works without a sibling community checkout. Java targets are
`src/*/java/**/*.java` inside each Java module; generated files under `build/` and
resource fixtures are outside the targets. No import sorting, unused-import removal,
license-header replacement, or member reordering is enabled.

The profile starts from community's modern Eclipse export and is adjusted toward the
shared IntelliJ `GoogleStyle` profile:

- Two-space indentation, four-space continuation indentation, and a 100-column limit.
- Wrap outer nested expressions and allow throws clauses to wrap.
- Preserve up to two existing blank lines, matching IntelliJ's declaration spacing more closely.
- Preserve comment and license-header text, including Javadocs, to avoid JDT-specific reflow
  of documentation already formatted by the engineers. New Javadocs are not reflowed by Spotless.
- Use explicit modern Java settings for records, switch expressions, and type annotations.

Existing IntelliJ and Eclipse IDE exports are unchanged. Exact agreement with IntelliJ's
native Reformat Code action is **not** guaranteed. Observed differences remain in nested
calls, chained calls, annotations, enum/empty-type spacing, text blocks, and some comments.
Use Gradle as the final formatter before committing; an IDE can invoke the selected-file
command above to use the exact same engine and profile as CI.

## Gradual adoption

The fixed adoption baseline is `f0a1a8301ba0596064613d512692f25ec6660d36`.
Only new or modified Java files relative to that baseline are checked and formatted.
The entire changed file is formatted, not only its edited lines. The baseline is fixed
so formatting stays enforced after merging a branch; it does not depend on a developer's
remote names or current `master` reference. Do not advance it to bypass formatting failures.

CI fetches full history (`fetch-depth: 0`). For a shallow local clone, fetch enough history
to include the baseline before running Spotless. To deliberately format an unchanged file,
make the intended edit first. No repository-wide reformat is part of this integration.

## Comparison with IntelliJ

The reference was IntelliJ IDEA **2026.2.3**, build **IU-262.10968.63**, using the unchanged
[mzmine-intellij-code-formater.xml](../mzmine-intellij-code-formater.xml) export. It was run only as a
comparison tool on temporary copies. It is not a build or CI dependency. This validates
against that version/profile, not an unidentified developer's personal IDE settings.

The initial sample set spans the Java modules, main/test sources, annotations, records,
sealed types, switches, unnamed variables, text blocks, and a synthetic Java 26 preview
fixture. Another 50 files were selected by sorting the SHA-256 hashes of their relative
paths, excluding the initial set and files larger than 30 KB. These additional files were
not used to tune the profile.

| Set | Files | Original Eclipse profile: differing lines | New profile: differing lines | New profile: byte-identical files |
| --- | ---: | ---: | ---: | ---: |
| Initial samples | 18 | 207 | 114 | 3/18 |
| Additional samples | 50 | 758 | 470 | 18/50 |

“Differing lines” counts insertions plus deletions against the independently formatted
IntelliJ output, excluding diff context. It is a churn measure, not a percentage of
semantic or complete formatter compatibility. The comparison shows both profile
differences and remaining formatter-layout differences.

Validation also exercises apply/check, repeat-apply idempotence, a deliberately malformed
formatting case, the selected-file hook, generated-file exclusion, and compilation of the
formatted Java 26 preview fixture. Application behavior and native IDE UI integration are
outside these formatting checks.

References: [Spotless Eclipse support](https://github.com/diffplug/spotless/blob/main/plugin-gradle/README.md#eclipse-jdt),
[gradual enforcement](https://github.com/diffplug/spotless/blob/main/plugin-gradle/README.md#how-can-i-enforce-formatting-gradually-aka-ratchet),
[selected-file integration](https://github.com/diffplug/spotless/blob/main/plugin-gradle/IDE_HOOK.md).
