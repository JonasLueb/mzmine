/*
 * Hyperparameter grid runner that executes EICResolverRunner multiple times within a single JVM.
 * Avoids repeated Gradle/Java process startup and speeds up sweeps dramatically.
 */
package io.github.mzmine.tools.eicresolver;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

public class EICGridRunner {

  public static void main(String[] args) throws Exception {
    final Path inputDir = Paths.get(System.getProperty("inputDir", "")).toAbsolutePath();
    if (Files.notExists(inputDir)) {
      System.err.println("Missing -DinputDir");
      return;
    }
    final Path outRoot = Paths.get(System.getProperty("outDir",
        inputDir.resolve("resolver_hparam").toString())).toAbsolutePath();
    Files.createDirectories(outRoot);

    final String format = System.getProperty("format", "html");
    final String iou = System.getProperty("iou", "0.30");
    final String rtTol = System.getProperty("rtTol", "0.02");
    final String minPoints = System.getProperty("minPoints", "3");

    // carry-through optional preprocessing overrides if provided
    final String preZeroQ = System.getProperty("preZeroQ", "");
    final String baselineFrac = System.getProperty("baselineFrac", "");

    // Grids (can be overridden with simple comma-separated -D lists)
    final List<String> ms_minAbs = listOrDefault("minAbsGrid", new String[]{"0","10","20"});
    final List<String> ms_thresh = listOrDefault("msThreshGrid", new String[]{"0.2","0.3","0.5"});
    final List<String> ms_ratio = listOrDefault("msRatioGrid", new String[]{"1.3","1.5"});
    final List<String> ms_rel = listOrDefault("msRelGrid", new String[]{"0.0","0.01"});
    final List<String> ms_search = listOrDefault("msSearchGrid", new String[]{"0.03","0.05"});
    final List<String> ms_durMin = listOrDefault("msDurMinGrid", new String[]{"0.0"});
    final List<String> ms_durMax = listOrDefault("msDurMaxGrid", new String[]{"10.0"});

    final List<String> sg_min = listOrDefault("sgMinGrid", new String[]{"0","10","20"});
    final List<String> sg_der = listOrDefault("sgDerivGrid", new String[]{"0.2","0.3","0.5"});
    final List<String> sg_durMin = listOrDefault("sgDurMinGrid", new String[]{"0.0"});
    final List<String> sg_durMax = listOrDefault("sgDurMaxGrid", new String[]{"10.0"});

    final List<String> cwt_snr = listOrDefault("cwtSnrGrid", new String[]{"0.8","1.0","1.5"});
    final List<String> cwt_per = listOrDefault("cwtPersistGrid", new String[]{"1","2","3"});
    final List<String> cwt_min = listOrDefault("cwtMinGrid", new String[]{"0"});
    final List<String> cwt_bound = listOrDefault("cwtBoundaryGrid", new String[]{"0.05"});
    final List<String> cwt_quant = listOrDefault("cwtQuantGrid", new String[]{"LOW_QUANTILE","NONE"});
    final List<String> cwt_q = listOrDefault("cwtQuantQGrid", new String[]{"0.1","0.2"});

    final List<RunCfg> cfgs = new ArrayList<>();
    // MinimumSearch
    for (String a : ms_minAbs) for (String t : ms_thresh) for (String rr: ms_ratio) for (String rel: ms_rel) for (String sr: ms_search) for (String d0: ms_durMin) for (String d1: ms_durMax) {
      cfgs.add(new RunCfg("MS_minAbs"+a+"_thr"+t+"_ratio"+rr+"_rel"+rel+"_sr"+sr+"_dur"+d0+"-"+d1, "MinimumSearch",
          mapOf(new String[][] {{"minAbs",a},{"msThresh",t},{"ratio",rr},{"minRel",rel},{"searchRange",sr},{"msDurMin",d0},{"msDurMax",d1}})));
    }
    // Savitzky-Golay
    for (String m : sg_min) for (String d : sg_der) for (String d0: sg_durMin) for (String d1: sg_durMax) {
      cfgs.add(new RunCfg("SG_min"+m+"_der"+d+"_dur"+d0+"-"+d1, "SavitzkyGolay",
          mapOf(new String[][] {{"sgMin",m},{"sgDeriv",d},{"sgDurMin",d0},{"sgDurMax",d1}})));
    }
    // CWT Ridge
    for (String s : cwt_snr) for (String p : cwt_per) for (String m : cwt_min) for (String b : cwt_bound) for (String qb : cwt_quant) for (String q : cwt_q) {
      cfgs.add(new RunCfg("CWT_snr"+s+"_per"+p+"_min"+m+"_b"+b+"_qb"+qb+"_q"+q, "CWT_Ridge",
          mapOf(new String[][] {{"cwtSnr",s},{"cwtPersist",p},{"cwtMin",m},{"cwtBoundary",b},{"cwtQuant",qb},{"cwtQuantQ",q}})));
    }

    // Optional quick filter: -Dresolver=Name to only run one family
    final String onlyResolver = System.getProperty("resolver", "").trim();

    for (RunCfg cfg : cfgs) {
      if (!onlyResolver.isEmpty() && !cfg.resolver.equalsIgnoreCase(onlyResolver)) continue;
      final Path outDir = outRoot.resolve(cfg.id);
      Files.createDirectories(outDir);

      // backup overridden keys (do not rely on inputDir/outDir system props; pass via args instead)
      final Map<String, String> backup = new LinkedHashMap<>();
      setProp(backup, "format", format);
      setProp(backup, "iou", iou);
      setProp(backup, "rtTol", rtTol);
      setProp(backup, "minPoints", minPoints);
      setProp(backup, "resolver", cfg.resolver);
      setProp(backup, "configId", cfg.id);
      if (!preZeroQ.isEmpty()) setProp(backup, "preZeroQ", preZeroQ);
      if (!baselineFrac.isEmpty()) setProp(backup, "baselineFrac", baselineFrac);
      for (Map.Entry<String,String> e : cfg.params.entrySet()) setProp(backup, e.getKey(), e.getValue());
      // special handling for CWT quant baseline selection
      if (cfg.params.containsKey("cwtQuant")) {
        final String mode = cfg.params.get("cwtQuant");
        // Use embedded parameter names from CwtRidgeResolverParameters via system props in EICResolverRunner
        System.setProperty("cwtQuantSel", mode);
      }

      System.out.printf(Locale.US, "Run %s (%s) -> %s\n", cfg.id, cfg.resolver, outDir);
      // pass input and output via args to avoid interference with system properties
      EICResolverRunner.main(new String[]{inputDir.toString(), outDir.toString()});

      // restore properties
      for (Map.Entry<String,String> e : backup.entrySet()) {
        if (e.getValue() == null) System.clearProperty(e.getKey());
        else System.setProperty(e.getKey(), e.getValue());
      }
    }

    // Aggregate results into search_results.csv and ranked_best.csv and index.html
    aggregate(outRoot);
    System.out.println("Done. Dashboard: " + outRoot.resolve("index.html"));
  }

  private static void aggregate(Path outRoot) throws IOException {
    final Path aggCsv = outRoot.resolve("search_results.csv");
    Files.writeString(aggCsv, "config\tresolver\ttp\tfp\tfn\tprecision\trecall\tf1\tweighted_f1\n", StandardCharsets.UTF_8);

    try (DirectoryStream<Path> ds = Files.newDirectoryStream(outRoot)) {
      for (Path run : ds) {
        if (!Files.isDirectory(run)) continue;
        final Path sum = run.resolve("evaluation_summary.csv");
        if (Files.notExists(sum)) continue;
        // sum across EICs by resolver
        final Map<String, long[]> totals = new LinkedHashMap<>();
        try (BufferedReader br = Files.newBufferedReader(sum, StandardCharsets.UTF_8)) {
          String line; boolean header = true;
          while ((line = br.readLine()) != null) {
            if (header) { header = false; continue; }
            if (line.isBlank()) continue;
            final String[] t = line.split("\t");
            if (t.length < 8) continue;
            final String resolver = t[1];
            final long tp = parseLong(t[2]);
            final long fp = parseLong(t[3]);
            final long fn = parseLong(t[4]);
            totals.computeIfAbsent(resolver, k -> new long[3]);
            totals.get(resolver)[0] += tp;
            totals.get(resolver)[1] += fp;
            totals.get(resolver)[2] += fn;
          }
        }
        for (Map.Entry<String,long[]> e : totals.entrySet()) {
          final long tp = e.getValue()[0], fp = e.getValue()[1], fn = e.getValue()[2];
          final double p = (tp+fp)==0?0:tp/(double)(tp+fp);
          final double r = (tp+fn)==0?0:tp/(double)(tp+fn);
          final double f1 = (p+r)==0?0:2*p*r/(p+r);
          final String row = String.format(Locale.US, "%s\t%s\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\n",
              run.getFileName(), e.getKey(), tp, fp, fn, p, r, f1, f1);
          Files.writeString(aggCsv, row, StandardCharsets.UTF_8, java.nio.file.StandardOpenOption.APPEND);
        }
      }
    }

    // Build ranked best per resolver and simple HTML
    final Path ranked = outRoot.resolve("ranked_best.csv");
    Files.writeString(ranked, "resolver\tconfig\ttp\tfp\tfn\tprecision\trecall\tf1\tweighted_f1\n", StandardCharsets.UTF_8);
    // Load agg and pick best by weighted_f1 (last column)
    final Map<String, String[]> best = new LinkedHashMap<>();
    try (BufferedReader br = Files.newBufferedReader(aggCsv, StandardCharsets.UTF_8)) {
      String line; boolean header = true;
      while ((line = br.readLine()) != null) {
        if (header) { header = false; continue; }
        if (line.isBlank()) continue;
        final String[] t = line.split("\t");
        final String config = t[0], resolver = t[1];
        final double wf = Double.parseDouble(t[8]);
        final String[] cur = best.get(resolver);
        if (cur == null || Double.parseDouble(cur[8]) < wf) best.put(resolver, t);
      }
    }
    for (String[] t : best.values()) {
      final String row = String.join("\t", Arrays.asList(t[1], t[0], t[2], t[3], t[4], t[5], t[6], t[7], t[8])) + "\n";
      Files.writeString(ranked, row, StandardCharsets.UTF_8, java.nio.file.StandardOpenOption.APPEND);
    }

    final StringBuilder html = new StringBuilder();
    html.append("<!doctype html><html><head><meta charset=\"utf-8\"><title>Resolver Hyperparameter Search</title>")
        .append("<style>body{font-family:Sans-Serif} table{border-collapse:collapse} td,th{border:1px solid #ccc;padding:4px 8px}</style>")
        .append("</head><body><h2>Resolver Hyperparameter Search</h2>")
        .append("<h3>Top configuration per resolver</h3>")
        .append("<table><thead><tr><th>Resolver</th><th>Config</th><th>TP</th><th>FP</th><th>FN</th><th>Precision</th><th>Recall</th><th>F1</th></tr></thead><tbody>");
    for (String[] t : best.values()) {
      final String resolver = t[1], cfg = t[0];
      html.append("<tr><td>").append(resolver).append("</td><td><a href=\"")
          .append(cfg).append("/index.html\">").append(cfg).append("</a></td><td>")
          .append(t[2]).append("</td><td>").append(t[3]).append("</td><td>")
          .append(t[4]).append("</td><td>").append(t[5]).append("</td><td>")
          .append(t[6]).append("</td><td>").append(t[7]).append("</td></tr>");
    }
    html.append("</tbody></table>")
        .append("<p>All results: <a href=\"search_results.csv\">search_results.csv</a> &nbsp; <a href=\"ranked_best.csv\">ranked_best.csv</a></p>")
        .append("</body></html>");
    Files.writeString(outRoot.resolve("index.html"), html.toString(), StandardCharsets.UTF_8);
  }

  private static List<String> listOrDefault(String key, String[] def) {
    final String v = System.getProperty(key, "");
    if (v == null || v.isBlank()) return Arrays.asList(def);
    return Arrays.stream(v.split(",")).map(String::trim).filter(s -> !s.isEmpty()).toList();
  }

  private static void setProp(Map<String,String> backup, String k, String v) {
    backup.put(k, System.getProperty(k));
    System.setProperty(k, v);
  }

  private static Map<String,String> mapOf(String[][] pairs) {
    final Map<String,String> m = new LinkedHashMap<>();
    for (String[] p : pairs) {
      if (p.length == 2) m.put(p[0], p[1]);
    }
    return m;
  }

  private static long parseLong(String s) {
    try { return Long.parseLong(s); } catch (Exception e) { return 0L; }
  }

  private record RunCfg(String id, String resolver, Map<String,String> params) {}
}


