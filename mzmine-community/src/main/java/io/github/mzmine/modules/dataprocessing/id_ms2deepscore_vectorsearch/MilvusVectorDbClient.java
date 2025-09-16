/*
 * Minimal Milvus client for ANN queries against MS2Deepscore embeddings.
 */
package io.github.mzmine.modules.dataprocessing.id_ms2deepscore_vectorsearch;

import io.milvus.v2.client.ConnectConfig;
import io.milvus.v2.client.MilvusClientV2;
import io.milvus.v2.service.vector.request.SearchReq;
import io.milvus.v2.service.vector.response.SearchResp;
import io.milvus.v2.service.vector.request.data.FloatVec;
import io.milvus.v2.service.vector.request.data.BaseVector;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.Collections;
import org.json.JSONObject;
import java.util.HashMap;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class MilvusVectorDbClient implements VectorDbClient {

  private MilvusClientV2 client;
  private String dbName; // optional
  private String user;
  private String password;
  private String idField = "id";
  private String vectorField = "vector";

  @Override
  public void connect(@NotNull String uri) throws IOException {
    try {
      String normalized = uri;
      if (normalized == null || normalized.isBlank()) {
        normalized = "http://localhost:19530";
      }
      URI u = new URI(normalized);
      if (u.getScheme() == null) {
        // accept bare host:port
        u = new URI("http://" + normalized);
      }
      final String scheme = u.getScheme();
      if (!("http".equalsIgnoreCase(scheme) || "https".equalsIgnoreCase(scheme) || "milvus".equalsIgnoreCase(scheme))) {
        throw new IOException("Unsupported URI scheme for Milvus: " + scheme);
      }
      final String host = (u.getHost() != null) ? u.getHost() : "localhost";
      final int port = (u.getPort() > 0) ? u.getPort() : 19530;
      final String q = u.getQuery();
      if (q != null && !q.isBlank()) {
        for (String kv : q.split("&")) {
          final String[] p = kv.split("=", 2);
          if (p.length == 2) {
            switch (p[0]) {
              case "db" -> dbName = p[1];
              case "user" -> user = p[1];
              case "password" -> password = p[1];
              case "idField" -> idField = p[1];
              case "vectorField" -> vectorField = p[1];
            }
          }
        }
      }
      String tokenParam = null;
      if (q != null && !q.isBlank()) {
        for (String kv : q.split("&")) {
          final String[] p = kv.split("=", 2);
          if (p.length == 2 && p[0].equals("token")) {
            tokenParam = java.net.URLDecoder.decode(p[1], java.nio.charset.StandardCharsets.UTF_8);
          }
        }
      }

      // normalize milvus:// to http:// expected by MilvusClientV2
      final String finalUri = ("milvus".equalsIgnoreCase(scheme) ? "http" : scheme) + "://" + host + ":" + port;
      ConnectConfig connectConfig = ConnectConfig.builder()
          .uri(finalUri)
          .token(tokenParam != null ? tokenParam : (user != null && password != null ? user + ":" + password : null))
          .build();
      client = new MilvusClientV2(connectConfig);
      if (dbName != null && !dbName.isBlank()) {
        client.useDatabase(dbName);
      }
    } catch (URISyntaxException e) {
      throw new IOException("Invalid Milvus URI: " + uri, e);
    } catch (Exception e) {
      throw new IOException("Failed to connect to Milvus: " + e.getMessage(), e);
    }
  }

  @Override
  public @NotNull List<VectorDbHit> query(@NotNull String collection, float[] queryVector, int topK,
      double minScore, @Nullable Map<String, Object> scalarFilters) throws IOException {
    if (client == null) throw new IOException("Not connected");
    try {
      // Build query vector (single query)
      List<BaseVector> data = Collections.singletonList(new FloatVec(queryVector));

      String filter = null;
      if (scalarFilters != null && scalarFilters.containsKey("precursor_mz_between")) {
        final Object filterValue = scalarFilters.get("precursor_mz_between");
        double lo, hi;
        if (filterValue instanceof double[] arr && arr.length == 2) { lo = arr[0]; hi = arr[1]; }
        else if (filterValue instanceof List<?> list && list.size() == 2) {
          lo = Double.parseDouble(String.valueOf(list.get(0)));
          hi = Double.parseDouble(String.valueOf(list.get(1)));
        } else { lo = Double.NaN; hi = Double.NaN; }
        if (!Double.isNaN(lo) && !Double.isNaN(hi)) {
          filter = "precursor_mz >= " + lo + " and precursor_mz <= " + hi;
        }
      }

      SearchReq req = SearchReq.builder()
          .collectionName(collection)
          .data(data)
          .topK(topK)
          .filter(filter)
          .outputFields(List.of("ext_id","precursor_mz","meta"))
          .build();

      SearchResp resp = client.search(req);
      final List<VectorDbHit> out = new ArrayList<>();
      // v2 resp provides nested results per query
      List<List<SearchResp.SearchResult>> results = resp.getSearchResults();
      if (results != null && !results.isEmpty()) {
        for (SearchResp.SearchResult r : results.get(0)) {
          final double score = r.getScore();
          if (Double.isNaN(score) || score < minScore) continue;
          final String idStr = String.valueOf(r.getId());

          // Extract metadata and precursor mz
          Map<String,Object> meta = new HashMap<>();
          Object extId = r.getEntity() != null ? r.getEntity().get("ext_id") : null;
          if (extId != null) meta.put("ext_id", extId);
          Double precMz = null;
          Object prec = r.getEntity() != null ? r.getEntity().get("precursor_mz") : null;
          if (prec instanceof Number n) precMz = n.doubleValue();
          Object metaStr = r.getEntity() != null ? r.getEntity().get("meta") : null;
          if (metaStr instanceof String s && !s.isBlank()) {
            try {
              JSONObject obj = new JSONObject(s);
              // Merge selected keys if present
              for (String k : List.of("name","compound_name","smiles","inchi","inchikey","formula","ionmode","collision_energy","source","db_id","accession","peaks","peaks_unit")) {
                if (obj.has(k) && !obj.isNull(k)) meta.put(k, obj.get(k));
              }
              // If only compound_name exists, mirror it to name for UI
              if (!meta.containsKey("name") && obj.has("compound_name") && !obj.isNull("compound_name")) {
                meta.put("name", obj.get("compound_name"));
              }
            } catch (Exception ignored) {}
          }

          out.add(new VectorDbHit(idStr, score, precMz, meta.isEmpty() ? null : meta));
        }
      }
      return out;
    } catch (Exception e) {
      throw new IOException("Milvus query failed: " + e.getMessage(), e);
    }
  }

  @Override
  public void close() throws IOException {
    try { if (client != null) client.close(); } catch (Exception ignored) {}
  }
}


