/*
 * Copyright (c) 2004-2026 The mzmine Development Team
 */
package io.github.mzmine.localbridge;

import java.util.Optional;
import java.util.UUID;
import org.jetbrains.annotations.NotNull;

/** An opaque, validated target received through the mzmine URL scheme. */
public record DeepLink(@NotNull UUID instanceId, @NotNull UUID projectId,
                       @NotNull Destination destination, @NotNull Optional<UUID> resourceId) {

  public DeepLink {
    if ((destination == Destination.WIZARD) != resourceId.isEmpty()) {
      throw new IllegalArgumentException("Only wizard destinations may omit a resource identifier");
    }
  }

  public enum Destination {
    WIZARD("wizard"), PROPOSAL("proposal"), FEATURE_TABLE("feature-table");

    private final String pathSegment;

    Destination(final @NotNull String pathSegment) {
      this.pathSegment = pathSegment;
    }

    static @NotNull Destination fromPathSegment(final @NotNull String pathSegment) {
      for (final Destination value : values()) {
        if (value.pathSegment.equals(pathSegment)) {
          return value;
        }
      }
      throw new IllegalArgumentException("Unknown mzmine deep-link destination: " + pathSegment);
    }
  }
}
