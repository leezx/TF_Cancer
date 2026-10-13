#!/usr/bin/env Rscript
# Step 2 design §3.4 (summary), §3.6 (tumor vs normal, bootstrap), §3.7 (pair-level
# recurrence, network), §3.8 (focused Fig 1B), §3.9 (BRCA PAM50), §4 (QC tables).
# Usage (repo root): Rscript scripts/02_tf_coexpression/03_integrate.R
source("scripts/02_tf_coexpression/common2.R")
suppressPackageStartupMessages(library(igraph))
cd <- file.path(res_dir, "coexpr")
cm <- fread("results/01_tcga_tf_de/cohort_manifest.tsv")[, .(cohort, organ_system)]
rd <- function(pat) rbindlist(lapply(cohorts, function(k) { f <- file.path(cd, sprintf(pat, k)); if (file.exists(f)) fread(f) }), fill = TRUE)
cq <- merge(rd("%s_cohort_qc.tsv"), cm, by = "cohort")
aq <- rd("%s_analysis_qc.tsv")
mats <- setNames(lapply(cohorts, function(k) readRDS(file.path(cd, sprintf("%s.rds", k)))), cohorts)
os_of <- setNames(cm$organ_system, cm$cohort)
analyses <- c("primary", "prolif_only", "unadjusted")

# ---- cohort QC ----
w <- dcast(aq, cohort ~ analysis, value.var = c("n_samples", "n_edges_pos", "n_edges_neg", "perm_edges_mean", "perm_edges_max"))
cq <- merge(cq, w, by = "cohort")
cq[, score_MKI67_check := fifelse(rho_score_MKI67 >= P$mki67_check, "pass", "below_0.7_investigate")]
fwrite_tsv(cq[order(match(cohort, cohorts))], file.path(repo_out, "cohort_qc.tsv"))
fwrite_tsv(aq, file.path(repo_out, "qc", "analysis_qc_permutation.tsv"))
fwrite_tsv(rd("%s_adjustment_check.tsv"), file.path(repo_out, "qc", "adjustment_check_E2F1_FOXM1.tsv"))

# ---- proliferation association + summary (§3.4) ----
pa <- rd("%s_prolif_assoc.tsv")
fwrite_tsv(pa, file.path(repo_out, "TF_proliferation_association.tsv"))
ps <- pa[, .(n_evaluable = .N, n_coupled = sum(coupled)), by = .(TF, definition = analysis)]
ps[, frac_coupled := n_coupled / n_evaluable]
ps[, recurrently_proliferation_coupled := frac_coupled >= P$recur_cohort_frac & n_coupled >= P$recur_min_cohorts]
fwrite_tsv(ps[order(definition, -n_coupled, TF)], file.path(repo_out, "TF_proliferation_summary.tsv"))

# ---- positive controls (§4.4) ----
ctl <- rbindlist(lapply(seq_len(nrow(controls)), function(i) {
  m <- mats[[controls$cohort[i]]]$mats; a <- controls$a[i]; b <- controls$b[i]
  rbindlist(lapply(names(m), function(nm) data.table(cohort = controls$cohort[i], pair = paste(a, b, sep = "-"), analysis = nm,
    rho = if (all(c(a, b) %in% rownames(m[[nm]]))) m[[nm]][a, b] else NA_real_)))
}))
ed_all <- rd("%s_edges.tsv.gz")
ctl <- merge(ctl, rbind(ed_all[, .(cohort, analysis, pair = paste(TF_a, TF_b, sep = "-"), q)], ed_all[, .(cohort, analysis, pair = paste(TF_b, TF_a, sep = "-"), q)]),
             by = c("cohort", "pair", "analysis"), all.x = TRUE)
ctl[, is_edge := !is.na(q)]
fwrite_tsv(ctl, file.path(repo_out, "qc", "positive_control_pairs.tsv"))

# ---- pair-level recurrence (§3.7) ----
all_tf <- sort(unique(unlist(lapply(mats, `[[`, "tf"))))
ix <- function(v) match(v, all_tf)
systems <- unique(os_of[cohorts])
rec_list <- list(); posmats <- list()
for (an in analyses) {
  n_eval <- n_pos <- matrix(0L, length(all_tf), length(all_tf), dimnames = list(all_tf, all_tf))
  sys_detail <- list()
  for (sy in systems) {
    ev <- ps_ <- matrix(FALSE, length(all_tf), length(all_tf))
    for (k in names(os_of)[os_of == sy & names(os_of) %in% cohorts]) {
      if (an == "primary" && !mats[[k]]$primary_eligible) next
      t <- ix(mats[[k]]$tf); ev[t, t] <- TRUE
      e <- ed_all[cohort == k & analysis == an & sign == "positive"]
      if (nrow(e)) { ps_[cbind(ix(e$TF_a), ix(e$TF_b))] <- TRUE; ps_[cbind(ix(e$TF_b), ix(e$TF_a))] <- TRUE }
    }
    n_eval <- n_eval + ev; n_pos <- n_pos + ps_
  }
  ut <- which(upper.tri(n_pos) & n_pos > 0, arr.ind = TRUE)
  rec_list[[an]] <- data.table(analysis = an, TF_a = all_tf[ut[, 1]], TF_b = all_tf[ut[, 2]],
                               n_organ_systems_positive = n_pos[ut], n_organ_systems_evaluable = n_eval[ut])
  log_msg("%s: %d pairs positive in >=1 system; %d recurrent (>=%d systems)", an, nrow(ut),
          sum(n_pos[ut] >= P$recur_organ_systems), P$recur_organ_systems)
}
rec <- rbindlist(rec_list)
rec[, frac_positive := n_organ_systems_positive / n_organ_systems_evaluable]
rec[, recurrent := n_organ_systems_positive >= P$recur_organ_systems]
# per-system detail for recurrent pairs: positive/evaluable cohorts and sign agreement (Q4).
# Vectorized: one long table of (pair, cohort) for cohorts where both TFs are analysed.
rr <- rec[recurrent == TRUE, .(analysis, TF_a, TF_b)]
long <- rbindlist(lapply(analyses, function(an) {
  pr_ <- rr[analysis == an]; if (!nrow(pr_)) return(NULL)
  rbindlist(lapply(cohorts, function(k) {
    if (an == "primary" && !mats[[k]]$primary_eligible) return(NULL)
    ok <- pr_$TF_a %in% mats[[k]]$tf & pr_$TF_b %in% mats[[k]]$tf
    if (!any(ok)) return(NULL)
    d <- pr_[ok]
    d[, `:=`(cohort = k, organ_system = os_of[[k]], rho = mats[[k]]$mats[[an]][cbind(TF_a, TF_b)])]
    d
  }))
}))
pos_keys <- ed_all[sign == "positive", paste(analysis, cohort, pmin(TF_a, TF_b), pmax(TF_a, TF_b))]
long[, pos := paste(analysis, cohort, pmin(TF_a, TF_b), pmax(TF_a, TF_b)) %chin% pos_keys]
ds <- long[, .(txt = sprintf("%s:%d/%d(sign+ %.2f)", organ_system[1], sum(pos), .N, mean(rho > 0))), by = .(analysis, TF_a, TF_b, organ_system)]
ds <- ds[, .(per_system = paste(txt, collapse = ";")), by = .(analysis, TF_a, TF_b)]
rec <- merge(rec, ds, by = c("analysis", "TF_a", "TF_b"), all.x = TRUE)
fwrite(rec[order(analysis, -n_organ_systems_positive, TF_a, TF_b)], file.path(res_dir, "coexpr", "TF_pair_recurrence_all.tsv.gz"), sep = "\t")
fwrite_tsv(rec[recurrent == TRUE][order(analysis, -n_organ_systems_positive, TF_a, TF_b)], file.path(res_dir, "report", "TF_pair_recurrence.tsv.gz"))

# ---- tumor vs normal (§3.6): prolif_only both sides ----
tn <- rbindlist(lapply(cohorts, function(k) {
  m <- mats[[k]]$mats; if (is.null(m$normal_prolif_only)) return(NULL)
  e <- ed_all[cohort == k & analysis == "prolif_only"]
  rn <- m$normal_prolif_only; nn <- cq[cohort == k, n_normal]
  ut <- upper.tri(rn); qn <- rn; qn[] <- NA; qn[ut] <- p.adjust(r2p(rn[ut], nn - 3), "BH")
  e[, `:=`(normal_rho = rn[cbind(TF_a, TF_b)], normal_q = pmin(qn[cbind(TF_a, TF_b)], qn[cbind(TF_b, TF_a)], na.rm = TRUE))]
  e[, delta_rho := rho - normal_rho]
  e[, shared_coexpression := abs(normal_rho) >= P$edge_rho & normal_q < P$fdr & sign(normal_rho) == sign(rho)]
  e
}))
# bootstrap CI of delta for primary-recurrent positive pairs (§3.6)
rp <- rec[analysis == "primary" & recurrent == TRUE, .(TF_a, TF_b)]
boot <- rbindlist(lapply(cohorts, function(k) {
  if (is.null(mats[[k]]$mats$normal_prolif_only)) return(NULL)
  pr <- readRDS(file.path(res_dir, "prep", sprintf("%s.rds", k)))
  pp <- rp[TF_a %in% pr$tf_symbols & TF_b %in% pr$tf_symbols]
  if (!nrow(pp)) return(NULL)
  g <- unique(c(pp$TF_a, pp$TF_b))
  Xt <- pr$tumor$X[, g, drop = FALSE]; Xn <- pr$normal$X[, g, drop = FALSE]
  st <- pr$tumor$samples$prolif; sn <- pr$normal$samples$prolif
  pc <- function(X, s) { E <- resid_ranks(X, matrix(s)); r <- suppressWarnings(cor(E)); r[cbind(pp$TF_a, pp$TF_b)] }
  set.seed(P$seed)
  B <- vapply(seq_len(P$n_boot), function(b) {
    it <- sample.int(nrow(Xt), replace = TRUE); inn <- sample.int(nrow(Xn), replace = TRUE)
    pc(Xt[it, , drop = FALSE], st[it]) - pc(Xn[inn, , drop = FALSE], sn[inn])
  }, numeric(nrow(pp)))
  B <- matrix(B, nrow = nrow(pp))
  tr <- pc(Xt, st); nr <- pc(Xn, sn)
  log_msg("bootstrap %s: %d recurrent pairs", k, nrow(pp))
  data.table(cohort = k, pp, tumor_rho_prolif_only = tr, normal_rho_prolif_only = nr, delta_rho = tr - nr,
             delta_ci_low = apply(B, 1, quantile, 0.025, na.rm = TRUE), delta_ci_high = apply(B, 1, quantile, 0.975, na.rm = TRUE))
}))
boot[, tumor_enriched_coexpression := delta_ci_low > 0 & normal_rho_prolif_only < P$tumor_enriched_normal_rho_max]
fwrite_tsv(boot, file.path(res_dir, "report", "TF_tumor_vs_normal_recurrent_pairs.tsv.gz"))

# ---- edges file: full in DATA; repo keeps edges of primary-recurrent pairs ----
ed_all <- merge(ed_all, tn[, .(cohort, analysis, TF_a, TF_b, normal_rho, normal_q, delta_rho, shared_coexpression)],
                by = c("cohort", "analysis", "TF_a", "TF_b"), all.x = TRUE)
fwrite(ed_all, file.path(res_dir, "coexpr", "TF_coexpression_edges_all.tsv.gz"), sep = "\t")
key <- paste(rp$TF_a, rp$TF_b)
fwrite_tsv(ed_all[paste(pmin(TF_a, TF_b), pmax(TF_a, TF_b)) %in% key], file.path(res_dir, "report", "TF_coexpression_edges_recurrent_pairs.tsv.gz"))

# ---- network + modules (§3.7) ----
ne <- rec[analysis == "primary" & recurrent == TRUE]
g <- graph_from_data_frame(ne[, .(TF_a, TF_b, weight = n_organ_systems_positive)], directed = FALSE)
set.seed(P$seed)
cl <- cluster_leiden(g, objective_function = "modularity", resolution = P$leiden_res, n_iterations = P$leiden_iter, weights = E(g)$weight)
s1 <- fread("results/01_tcga_tf_de/TF_recurrent_upregulation.tsv")
tfr <- fread("results/01_tcga_tf_de/TF_reference.tsv")
nodes <- data.table(TF = V(g)$name, module = membership(cl), degree = degree(g), weighted_degree = strength(g))
nodes <- merge(nodes, tfr[, .(TF = hgnc_symbol, DBD)], by = "TF", all.x = TRUE)
nodes <- merge(nodes, s1[, .(TF, n_organ_systems_C1_up, fig1B_primary, pipeline_sensitive)], by = "TF", all.x = TRUE)
nodes <- merge(nodes, ps[definition == "purity_adjusted", .(TF, recurrently_proliferation_coupled)], by = "TF", all.x = TRUE)
msz <- nodes[, .N, by = module][N >= P$module_min, module]
nodes[, module_reported := module %in% msz]
fwrite_tsv(nodes[order(module, -weighted_degree)], file.path(repo_out, "TF_network_modules.tsv"))
el <- as.data.table(as_edgelist(g)); setnames(el, c("a", "b")); el[, `:=`(ma = membership(cl)[a], mb = membership(cl)[b])]
msum <- nodes[module_reported == TRUE, .(n_nodes = .N, top_TFs = paste(head(TF[order(-weighted_degree)], 8), collapse = ","),
                                         n_prolif_coupled = sum(recurrently_proliferation_coupled %in% TRUE),
                                         n_fig1B_primary = sum(fig1B_primary %in% TRUE)), by = module]
msum <- merge(msum, el[ma == mb, .(n_edges = .N), by = .(module = ma)], by = "module", all.x = TRUE)
fwrite_tsv(msum[order(-n_nodes)], file.path(repo_out, "TF_network_module_summary.tsv"))
log_msg("network: %d nodes, %d edges, %d modules >= %d TFs", vcount(g), ecount(g), length(msz), P$module_min)

# SOX9 reference (CRC only, descriptive)
sx <- ed_all[cohort == "COADREAD" & analysis == "primary" & sign == "positive" & (TF_a == "SOX9" | TF_b == "SOX9")]
sx[, partner := fifelse(TF_a == "SOX9", TF_b, TF_a)]
sx <- merge(sx[, .(partner, rho, q)], rec[analysis == "primary", .(partner = fifelse(TF_a == "SOX9", TF_b, TF_a), n_organ_systems_positive, n_organ_systems_evaluable, sox = TF_a == "SOX9" | TF_b == "SOX9")][sox == TRUE, -"sox"],
            by = "partner", all.x = TRUE)
fwrite_tsv(sx[order(-rho)], file.path(repo_out, "SOX9_COADREAD_reference_edges.tsv"))

# ---- focused Fig 1B (§3.8) ----
foc <- merge(s1[fig1B_primary == TRUE, .(TF, n_organ_systems_C1_up, n_organ_systems_C1_tier1, pipeline_sensitive)],
             ps[definition == "purity_adjusted", .(TF, n_evaluable, n_coupled, recurrently_proliferation_coupled)], by = "TF", all.x = TRUE)
foc[is.na(recurrently_proliferation_coupled), recurrently_proliferation_coupled := FALSE]
foc[, in_focused_view := !recurrently_proliferation_coupled]
fwrite_tsv(foc[order(-in_focused_view, -n_organ_systems_C1_up, TF)], file.path(repo_out, "TF_focused_fig1B.tsv"))
log_msg("focused Fig 1B: %d of %d Fig1B-primary TFs not recurrently proliferation-coupled", sum(foc$in_focused_view), nrow(foc))

# ---- BRCA PAM50 sensitivity (§3.9) ----
pr <- readRDS(file.path(res_dir, "prep", "BRCA.rds")); sm <- pr$tumor$samples; X <- pr$tumor$X
base <- ed_all[cohort == "BRCA" & analysis == "primary" & sign == "positive", paste(TF_a, TF_b)]
sens <- list()
ok <- !is.na(sm$purity) & sm$pam50 %in% c("LumA", "LumB", "Her2", "Basal")
Zs <- cbind(as.matrix(sm[ok, .(purity, prolif)]), model.matrix(~ factor(pam50, levels = c("LumA", "LumB", "Her2", "Basal")), sm[ok])[, -1])
sens$pam50_adjusted <- list(X = X[ok, ], Z = Zs)
okl <- !is.na(sm$purity) & sm$pam50 %in% "LumA"
if (sum(okl) >= P$luma_n_min) sens$luma_only <- list(X = X[okl, ], Z = as.matrix(sm[okl, .(purity, prolif)]))
bs <- rbindlist(lapply(names(sens), function(nm) {
  # factor indicators are not ranked: residualize on [ranked purity, ranked prolif, raw indicators]
  Zr <- cbind(1, apply(sens[[nm]]$Z[, 1:2, drop = FALSE], 2, rank), sens[[nm]]$Z[, -(1:2), drop = FALSE])
  E <- qr.resid(qr(Zr), apply(sens[[nm]]$X, 2, rank)); r <- suppressWarnings(cor(E)); df <- nrow(E) - 2 - (ncol(Zr) - 1)
  ut <- which(upper.tri(r), arr.ind = TRUE); rho <- r[ut]; q <- p.adjust(r2p(rho, df), "BH")
  pos <- paste(rownames(r)[ut[, 1]], colnames(r)[ut[, 2]])[rho >= P$edge_rho & q < P$fdr]
  data.table(analysis = nm, n = nrow(E), n_pos_edges = length(pos), n_pos_edges_BRCA_primary = length(base),
             jaccard_vs_primary = length(intersect(pos, base)) / length(union(pos, base)),
             rho_FOXA1_GATA3 = r["FOXA1", "GATA3"])
}))
fwrite_tsv(bs, file.path(repo_out, "BRCA_PAM50_sensitivity.tsv"))
print(bs)

# primary vs unadjusted overlap (§4.7)
jac <- ed_all[sign == "positive", .(set = list(paste(TF_a, TF_b))), by = .(cohort, analysis)]
jj <- jac[, { s <- setNames(set, analysis); data.table(
  jaccard_primary_vs_unadjusted = if (!is.null(s$primary)) length(intersect(s$primary, s$unadjusted)) / length(union(s$primary, s$unadjusted)) else NA_real_,
  jaccard_prolif_only_vs_unadjusted = length(intersect(s$prolif_only, s$unadjusted)) / length(union(s$prolif_only, s$unadjusted))) }, by = cohort]
fwrite_tsv(jj, file.path(repo_out, "qc", "edge_overlap_adjusted_vs_unadjusted.tsv"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "report", "sessionInfo_03_integrate.txt"))
