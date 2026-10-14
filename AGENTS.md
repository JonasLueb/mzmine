# Java formatting

- Format edited Java files using this repository's Gradle wrapper; do not imitate the style by hand.
- For one file, run `./gradlew :mzmine-community:spotlessApply -PspotlessIdeHook=/absolute/path/File.java`, substituting the owning module.
- Before submitting changes, run `./gradlew spotlessCheck` from this repository root.
- See [docs/formatting.md](docs/formatting.md) for formatter versions, profiles, and the Git baseline.
