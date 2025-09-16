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

package io.github.mzmine.modules.dataprocessing.id_ms2deepscore_vectorsearch;

import java.io.Closeable;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * Minimal vector DB client abstraction to support ANN queries against MS2Deepscore embeddings.
 * Implementations may wrap Milvus, Qdrant, FAISS, SQLite, etc.
 */
public interface VectorDbClient extends Closeable {

  /**
   * Open or attach to the database given a URI (file path or connection string).
   */
  void connect(@NotNull String uri) throws IOException;

  /**
   * @param collection name of the collection (e.g., ionmode-specific)
   * @param queryVector embedding vector
   * @param topK number of nearest neighbors to return
   * @param minScore minimum cosine similarity to retain
   * @param scalarFilters optional scalar filters (e.g., {"precursor_mz_between": [lo, hi]})
   */
  @NotNull
  List<io.github.mzmine.modules.dataprocessing.id_ms2deepscore_vectorsearch.VectorDbHit> query(@NotNull String collection, float[] queryVector, int topK,
      double minScore, @Nullable Map<String, Object> scalarFilters) throws IOException;

  @Override
  void close() throws IOException;
}


