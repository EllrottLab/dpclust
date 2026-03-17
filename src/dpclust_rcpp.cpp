#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::export]]
int omp_thread_count_cpp(int requested = -1) {
#ifdef _OPENMP
    if (requested > 0) {
        omp_set_num_threads(requested);
        return requested;
    }
    return omp_get_max_threads();
#else
    return 0; // OpenMP not compiled in
#endif
}

// -----------------------------------------------------------------------------
// OPTIMIZATION (10, 11): Inline helpers and simplified math
// -----------------------------------------------------------------------------
inline double clamp_prob(double x) {
    if (x < 1e-6) return 1e-6;
    if (x > 0.999999) return 0.999999;
    return x;
}

double mutationBurdenToMutationCopyNumber(double burden, double totalCopyNumber, double cellularity, double normalCopyNumber) {
    if (std::abs(cellularity) < 1e-9) return 0.0;
    double denom = cellularity * totalCopyNumber + normalCopyNumber * (1.0 - cellularity);
    double mutCopyNumber = burden / cellularity * denom;
    return (R_IsNaN(mutCopyNumber)) ? 0.0 : mutCopyNumber;
}

// -----------------------------------------------------------------------------
// OPTIMIZATION (5): Contiguous data layout for cache-friendly access
// -----------------------------------------------------------------------------
struct MutationRow {
    double mut;
    double wt;
    double mb_unit;
};

// [[Rcpp::export]]
List subclone_dirichlet_gibbs_cpp(NumericMatrix mutCount, NumericMatrix WTCount, 
                                  NumericMatrix totalCopyNumber, NumericMatrix normalCopyNumber, 
                                  NumericMatrix copyNumberAdjustment, 
                                  int C, NumericVector cellularity, int iter, 
                                  double conc_param, double cluster_conc,
                                  bool keep_aux_fields,
                                  int num_threads,
                                  IntegerVector stored_iters,
                                  IntegerVector conflict_i = IntegerVector::create(), 
                                  IntegerVector conflict_j = IntegerVector::create(), 
                                  NumericVector conflict_w = NumericVector::create(),
                                  Function log_func = R_NilValue) {
    
    int num_muts = mutCount.nrow();
    int num_timepoints = mutCount.ncol();
    
    std::vector<std::vector<std::pair<int, double>>> conflicts_adj(num_muts);
    if (conflict_i.size() > 0) {
        for (int k = 0; k < conflict_i.size(); ++k) {
            int u = conflict_i[k] - 1;
            int v = conflict_j[k] - 1;
            if (u >= 0 && u < num_muts && v >= 0 && v < num_muts) {
                if (conflict_w[k] > 1.0) conflicts_adj[u].push_back({v, std::log(conflict_w[k])});
            }
        }
    }

    int active_threads = omp_thread_count_cpp(num_threads);

    std::vector<std::mt19937> thread_rngs(static_cast<size_t>(active_threads));
    for (int t = 0; t < active_threads; ++t) {
        thread_rngs[t].seed(static_cast<uint32_t>(R::runif(0, 4294967295.0)));
    }
    
    double A = 1.0;
    double B = conc_param;
    
    NumericVector pi_h(iter * C * num_timepoints);
    NumericMatrix V_h(iter, C);
    std::fill(V_h.begin(), V_h.end(), 1.0);
    
    std::vector<int> keep_iters_zero_based;
    if (stored_iters.size() > 0) {
        std::vector<bool> seen(iter, false);
        for (int i = 0; i < stored_iters.size(); ++i) {
            int iter_idx = stored_iters[i] - 1;
            if (iter_idx >= 0 && iter_idx < iter && !seen[iter_idx]) {
                seen[iter_idx] = true;
                keep_iters_zero_based.push_back(iter_idx);
            }
        }
        std::sort(keep_iters_zero_based.begin(), keep_iters_zero_based.end());
    }
    bool store_all_iters = keep_iters_zero_based.empty();
    int stored_rows = store_all_iters ? iter : static_cast<int>(keep_iters_zero_based.size());
    IntegerMatrix S_i(stored_rows, num_muts);
    IntegerVector stored_iters_out;
    std::vector<int> iter_to_store_index(iter, -1);
    if (store_all_iters) {
        for (int i = 0; i < iter; ++i) iter_to_store_index[i] = i;
    } else {
        stored_iters_out = IntegerVector(stored_rows);
        for (int i = 0; i < stored_rows; ++i) {
            iter_to_store_index[keep_iters_zero_based[i]] = i;
            stored_iters_out[i] = keep_iters_zero_based[i] + 1;
        }
    }
    NumericVector alpha(iter);
    
    NumericVector lower(num_timepoints);
    NumericVector upper(num_timepoints);
    
    // -------------------------------------------------------------------------
    // OPTIMIZATION (5, 12): Contiguous row-major data layout
    // -------------------------------------------------------------------------
    std::vector<MutationRow> mut_data(static_cast<size_t>(num_muts) * num_timepoints);
    
    for (int t = 0; t < num_timepoints; ++t) {
        const double cell_t = cellularity[t];
        const double* mut_ptr = &mutCount(0, t);
        const double* wt_ptr = &WTCount(0, t);
        const double* tcn_ptr = &totalCopyNumber(0, t);
        const double* ncn_ptr = &normalCopyNumber(0, t);
        const double* cna_ptr = &copyNumberAdjustment(0, t);

        for (int k = 0; k < num_muts; ++k) {
            MutationRow& row = mut_data[static_cast<size_t>(k) * num_timepoints + t];
            row.mut = mut_ptr[k];
            row.wt = wt_ptr[k];
            double denom = cell_t * tcn_ptr[k] + ncn_ptr[k] * (1.0 - cell_t);
            double mb_unit = (std::abs(denom) < 1e-9) ? 0.000001 : (cell_t / denom) * cna_ptr[k];
            row.mb_unit = (std::isfinite(mb_unit) && mb_unit > 1e-6) ? mb_unit : 1e-6;
        }
    }

    if (log_func != R_NilValue) {
        log_func("Starting Gibbs sampler for " + std::to_string(num_muts) + " mutations...");
    } else {
        Rcout << "Starting Gibbs sampler for " << num_muts << " mutations..." << std::endl;
    }
    
    // Initialization
    for (int t = 0; t < num_timepoints; ++t) {
        double min_val = R_PosInf;
        double max_val = R_NegInf;
        for (int k = 0; k < num_muts; ++k) {
            const MutationRow& row = mut_data[static_cast<size_t>(k) * num_timepoints + t];
            double burden = row.mut / (row.mut + row.wt);
            if (!std::isfinite(burden)) burden = 0; 
            double mcn = mutationBurdenToMutationCopyNumber(burden, totalCopyNumber(k, t), cellularity[t], normalCopyNumber(k, t));
            mcn /= copyNumberAdjustment(k, t);
            if (mcn < min_val) min_val = mcn;
            if (mcn > max_val) max_val = mcn;
        }
        lower[t] = min_val; upper[t] = max_val;
        double diff = upper[t] - lower[t];
        lower[t] -= diff / 10.0; upper[t] += diff / 10.0;
        std::uniform_real_distribution<double> start_dist(lower[t], upper[t]);
        for (int c = 0; c < C; ++c) pi_h[0 + iter * c + iter * C * t] = start_dist(thread_rngs[0]);
    }
    for (int c = 0; c < C - 1; ++c) V_h(0, c) = 0.5;
    V_h(0, C - 1) = 1.0;
    std::vector<int> S_curr(num_muts, 1);
    if (iter_to_store_index[0] >= 0) {
        int row = iter_to_store_index[0];
        for (int k = 0; k < num_muts; ++k) S_i(row, k) = S_curr[k];
    }
    alpha[0] = 1.0;
    
    // OPTIMIZATION (1, 2, 13): Cache Line Padding and Buffer Reduction
    int C_padded = ((C + 7) / 8) * 8;
    int ct_padded = (((C * num_timepoints) + 7) / 8) * 8;
    int ct = C * num_timepoints;

    std::vector<double> Pr_S_threads(static_cast<size_t>(active_threads) * C_padded, 0.0);
    std::vector<double> shape_sums_thread(static_cast<size_t>(active_threads) * ct_padded, 0.0);
    std::vector<double> rate_sums_thread(static_cast<size_t>(active_threads) * ct_padded, 0.0);
    
    std::vector<double> shape_sums(ct, 0.0);
    std::vector<double> rate_sums(ct, 0.0);
    std::vector<int> cluster_counts(C, 0);

    for (int m = 1; m < iter; ++m) {
        if ((m + 1) % 100 == 0) {
            if (log_func != R_NilValue) log_func("Iteration " + std::to_string(m + 1) + " / " + std::to_string(iter));
            else Rcout << "Iteration " << m + 1 << " / " << iter << std::endl;
            Rcpp::checkUserInterrupt(); 
        }
        
        std::vector<double> log_prior(C);
        log_prior[0] = std::log(V_h(m - 1, 0));
        double sum_log_1_minus_V = 0.0;
        for (int j = 1; j < C; ++j) {
            sum_log_1_minus_V += std::log1p(-V_h(m - 1, j - 1)); // Suggestion (9)
            log_prior[j] = std::log(V_h(m - 1, j)) + sum_log_1_minus_V;
        }

        // OPTIMIZATION (8): Compact current slice
        std::vector<double> current_ccfs(ct);
        for (int c = 0; c < C; ++c) {
            for (int t = 0; t < num_timepoints; ++t) {
                current_ccfs[c + C * t] = pi_h[(m - 1) + iter * c + iter * C * t];
            }
        }

        // OPTIMIZATION (1): Clear buffers using fast memory functions
        std::fill(shape_sums_thread.begin(), shape_sums_thread.end(), 0.0);
        std::fill(rate_sums_thread.begin(), rate_sums_thread.end(), 0.0);

#ifdef _OPENMP
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            double* Pr_S = &Pr_S_threads[static_cast<size_t>(tid) * C_padded];
            std::mt19937& rng = thread_rngs[tid];
            std::uniform_real_distribution<double> udist(0.0, 1.0);

#pragma omp for schedule(static)
            for (int k = 0; k < num_muts; ++k) {
                for (int c = 0; c < C; ++c) Pr_S[c] = log_prior[c];
                
                size_t mut_offset = static_cast<size_t>(k) * num_timepoints;
                for (int t = 0; t < num_timepoints; ++t) {
                    const MutationRow& row = mut_data[mut_offset + t];
                    const double* ccfs_t = &current_ccfs[C * t];
                    
                    for (int c = 0; c < C; ++c) {
                        double mb = clamp_prob(ccfs_t[c] * row.mb_unit); // Suggestion (10)
                        Pr_S[c] += row.mut * std::log(mb) + row.wt * std::log1p(-mb); // Suggestion (9)
                    }
                }
                
                if (!conflicts_adj[k].empty()) {
                    for (auto& edge : conflicts_adj[k]) {
                        int neighbor_cluster = S_curr[edge.first] - 1;
                        if (neighbor_cluster >= 0 && neighbor_cluster < C) Pr_S[neighbor_cluster] -= edge.second;
                    }
                }
                
                double max_val = Pr_S[0];
                for(int c=1; c<C; ++c) if(Pr_S[c] > max_val) max_val = Pr_S[c];
                double sum_exp = 0.0;
                for(int c=0; c<C; ++c) {
                    double ex = std::exp(std::max(Pr_S[c] - max_val, -700.0));
                    Pr_S[c] = ex;
                    sum_exp += ex;
                }
                
                if (sum_exp <= 0.0 || !std::isfinite(sum_exp)) {
                    sum_exp = (double)C;
                    for (int c = 0; c < C; ++c) Pr_S[c] = 1.0;
                }
                
                double r = udist(rng) * sum_exp;
                double cum_sum = 0.0;
                int picked = C - 1;
                for(int c = 0; c < C; ++c) {
                    cum_sum += Pr_S[c];
                    if (r <= cum_sum) { picked = c; break; }
                }
                S_curr[k] = picked + 1;
            }
        }
#else
        std::vector<double> Pr_S_vec(C);
        double* Pr_S = Pr_S_vec.data();
        for (int k = 0; k < num_muts; ++k) {
            for (int c = 0; c < C; ++c) Pr_S[c] = log_prior[c];
            size_t mut_offset = static_cast<size_t>(k) * num_timepoints;
            for (int t = 0; t < num_timepoints; ++t) {
                const MutationRow& row = mut_data[mut_offset + t];
                const double* ccfs_t = &current_ccfs[C * t];
                for (int c = 0; c < C; ++c) {
                    double mb = clamp_prob(ccfs_t[c] * row.mb_unit);
                    Pr_S[c] += row.mut * std::log(mb) + row.wt * std::log1p(-mb);
                }
            }
            if (!conflicts_adj[k].empty()) {
                for (auto& edge : conflicts_adj[k]) {
                    int neighbor_cluster = S_curr[edge.first] - 1;
                    if (neighbor_cluster >= 0 && neighbor_cluster < C) Pr_S[neighbor_cluster] -= edge.second;
                }
            }
            double max_val = Pr_S[0];
            for(int c=1; c<C; ++c) if(Pr_S[c] > max_val) max_val = Pr_S[c];
            double sum_exp = 0.0;
            for(int c=0; c<C; ++c) { Pr_S[c] = std::exp(std::max(Pr_S[c] - max_val, -700.0)); sum_exp += Pr_S[c]; }
            double r = R::runif(0, 1) * sum_exp;
            double cum_sum = 0.0;
            int picked = C - 1;
            for(int c = 0; c < C; ++c) { cum_sum += Pr_S[c]; if (r <= cum_sum) { picked = c; break; } }
            S_curr[k] = picked + 1;
        }
#endif
        
        std::fill(cluster_counts.begin(), cluster_counts.end(), 0);
        for (int k = 0; k < num_muts; ++k) cluster_counts[S_curr[k] - 1]++;
        
        int cumulative_count = 0;
        for (int c = 0; c < C - 1; ++c) {
            int count_eq = cluster_counts[c];
            cumulative_count += count_eq;
            int count_gt = num_muts - cumulative_count;
            V_h(m, c) = R::rbeta(1.0 + count_eq, alpha[m-1] + count_gt);
            if(V_h(m, c) >= 1.0) V_h(m, c) = 0.999;
        }
        V_h(m, C - 1) = 1.0;
        
        for(int t=0; t<num_timepoints; ++t) {
             for(int c=0; c<C; ++c) pi_h[m + iter * c + iter * C * t] = R::runif(lower[t], upper[t]);
        }
        
#ifdef _OPENMP
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            double* sh_loc = &shape_sums_thread[static_cast<size_t>(tid) * ct_padded];
            double* rt_loc = &rate_sums_thread[static_cast<size_t>(tid) * ct_padded];
#pragma omp for schedule(static)
            for (int k = 0; k < num_muts; ++k) {
                int c = S_curr[k] - 1;
                size_t mut_offset = static_cast<size_t>(k) * num_timepoints;
                for (int t = 0; t < num_timepoints; ++t) {
                    const MutationRow& row = mut_data[mut_offset + t];
                    sh_loc[c + C * t] += row.mut;
                    rt_loc[c + C * t] += (row.mut + row.wt) * row.mb_unit;
                }
            }
        }
        // OPTIMIZATION (2): Parallel reduction over indices
#pragma omp parallel for schedule(static)
        for (int idx = 0; idx < ct; ++idx) {
            int t = idx / C; int c = idx % C;
            double s = 0, r = 0;
            for (int tid = 0; tid < active_threads; ++tid) {
                s += shape_sums_thread[tid * ct_padded + (c + C * t)];
                r += rate_sums_thread[tid * ct_padded + (c + C * t)];
            }
            shape_sums[idx] = s;
            rate_sums[idx] = r;
        }
#else
        std::fill(shape_sums.begin(), shape_sums.end(), 0.0);
        std::fill(rate_sums.begin(), rate_sums.end(), 0.0);
        for (int k = 0; k < num_muts; ++k) {
            int c = S_curr[k] - 1;
            size_t mut_offset = static_cast<size_t>(k) * num_timepoints;
            for (int t = 0; t < num_timepoints; ++t) {
                const MutationRow& row = mut_data[mut_offset + t];
                shape_sums[c + C * t] += row.mut;
                rate_sums[c + C * t] += (row.mut + row.wt) * row.mb_unit;
            }
        }
#endif

        for (int c = 0; c < C; ++c) {
            if (cluster_counts[c] == 0) continue;
            for (int t = 0; t < num_timepoints; ++t) {
                int idx = c + C * t;
                double s = shape_sums[idx]; double r = rate_sums[idx];
                pi_h[m + iter * c + iter * C * t] = (r == 0) ? 0.0 : R::rgamma(s, 1.0 / r);
            }
        }

        if (iter_to_store_index[m] >= 0) {
            int row = iter_to_store_index[m];
            for (int k = 0; k < num_muts; ++k) S_i(row, k) = S_curr[k];
        }
        double sum_log = 0.0;
        for(int j=0; j<C-1; ++j) sum_log += std::log1p(-V_h(m, j));
        alpha[m] = R::rgamma(C + A - 1, 1.0 / (B - sum_log));
    }
    
    pi_h.attr("dim") = Dimension(iter, C, num_timepoints);
    return List::create(Named("S.i") = S_i, Named("V.h") = V_h, Named("pi.h") = pi_h, Named("stored_iters") = (store_all_iters ? R_NilValue : Rcpp::wrap(stored_iters_out)), Named("alpha") = alpha);
}

// [[Rcpp::export]]
NumericMatrix assign_mutations_1d_cpp(IntegerMatrix S_i, NumericMatrix pi_h, NumericVector boundary, IntegerVector sampledIters_pi, IntegerVector sampledIters_state, int num_threads = -1) {
    if (sampledIters_pi.size() != sampledIters_state.size()) stop("sampledIters_pi and sampledIters_state must have equal length.");
    int num_muts = S_i.ncol(); int num_optima = boundary.size() + 1; int num_sampled = sampledIters_pi.size();
    NumericMatrix mutation_preferences(num_muts, num_optima);
    omp_thread_count_cpp(num_threads);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < num_muts; ++k) {
        for (int s_idx = 0; s_idx < num_sampled; ++s_idx) {
            int s_pi = sampledIters_pi[s_idx] - 1, s_state = sampledIters_state[s_idx] - 1;
            int c = S_i(s_state, k); if (c <= 0) continue;
            double val = pi_h(s_pi, c - 1); int opt = 0;
            for (int b = 0; b < boundary.size(); ++b) { if (val > boundary[b]) opt++; }
            mutation_preferences(k, opt) += 1.0;
        }
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < num_muts; ++k) { for (int o = 0; o < num_optima; ++o) mutation_preferences(k, o) /= num_sampled; }
    return mutation_preferences;
}

// [[Rcpp::export]]
NumericMatrix assign_mutations_nd_cpp(IntegerMatrix S_i, NumericVector pi_h_flat, IntegerVector pi_h_dims, NumericMatrix boundary, NumericVector plane_vector_flat, NumericMatrix vector_length, LogicalMatrix vector_direction, IntegerVector sampledIters_pi, IntegerVector sampledIters_state, int num_threads = -1) {
    if (sampledIters_pi.size() != sampledIters_state.size()) stop("sampledIters_pi and sampledIters_state must have equal length.");
    const int num_muts = S_i.ncol(), no_iters = pi_h_dims[0], C_total = pi_h_dims[1], no_subsamples = pi_h_dims[2], no_optima = vector_length.nrow(), num_sampled = sampledIters_pi.size();
    NumericMatrix mutation_preferences(num_muts, no_optima);
    omp_thread_count_cpp(num_threads);
#ifdef _OPENMP
#pragma omp parallel
    {
        std::vector<double> votes(no_optima, 0.0);
        std::vector<double> pi_slice(no_subsamples); // Suggestion (14)
#pragma omp for schedule(static)
        for (int k = 0; k < num_muts; ++k) {
            for (int s_idx = 0; s_idx < num_sampled; ++s_idx) {
                int s_pi = sampledIters_pi[s_idx] - 1, s_state = sampledIters_state[s_idx] - 1;
                int c_1based = S_i(s_state, k); if (c_1based <= 0) continue;
                int c = c_1based - 1;
                
                // PERFORMANCE (14): Gather pi slice once
                for (int t = 0; t < no_subsamples; ++t) {
                    pi_slice[t] = pi_h_flat[(size_t)s_pi + (size_t)no_iters * c + (size_t)no_iters * C_total * t];
                }
                
                std::fill(votes.begin(), votes.end(), 0.0);
                for (int i = 0; i < no_optima - 1; ++i) {
                    for (int j = i + 1; j < no_optima; ++j) {
                        double dist = 0;
                        size_t plane_offset = (size_t)i + (size_t)no_optima * j;
                        size_t plane_stride = (size_t)no_optima * no_optima;
                        for (int t = 0; t < no_subsamples; ++t) {
                            dist += pi_slice[t] * plane_vector_flat[plane_offset + (size_t)t * plane_stride];
                        }
                        dist /= vector_length(i, j);
                        bool lead_to_i = (dist <= boundary(i, j)) ? vector_direction(i, j) : !vector_direction(i, j);
                        if (lead_to_i) votes[i]++; else votes[j]++;
                    }
                }
                int best_opt = 0; double max_v = -1;
                for (int o = 0; o < no_optima; ++o) { if (votes[o] > max_v) { max_v = votes[o]; best_opt = o; } }
                mutation_preferences(k, best_opt) += 1.0;
            }
        }
    }
#else
    std::vector<double> votes(no_optima, 0.0);
    std::vector<double> pi_slice(no_subsamples);
    for (int k = 0; k < num_muts; ++k) {
        for (int s_idx = 0; s_idx < num_sampled; ++s_idx) {
            int s_pi = sampledIters_pi[s_idx] - 1, s_state = sampledIters_state[s_idx] - 1;
            int c_1based = S_i(s_state, k); if (c_1based <= 0) continue;
            int c = c_1based - 1;
            for (int t = 0; t < no_subsamples; ++t) pi_slice[t] = pi_h_flat[(size_t)s_pi + (size_t)no_iters * c + (size_t)no_iters * C_total * t];
            std::fill(votes.begin(), votes.end(), 0.0);
            for (int i = 0; i < no_optima - 1; ++i) {
                for (int j = i + 1; j < no_optima; ++j) {
                    double dist = 0;
                    size_t plane_offset = (size_t)i + (size_t)no_optima * j;
                    size_t plane_stride = (size_t)no_optima * no_optima;
                    for (int t = 0; t < no_subsamples; ++t) {
                        dist += pi_slice[t] * plane_vector_flat[plane_offset + (size_t)t * plane_stride];
                    }
                    dist /= vector_length(i, j);
                    bool lead_to_i = (dist <= boundary(i, j)) ? vector_direction(i, j) : !vector_direction(i, j);
                    if (lead_to_i) votes[i]++; else votes[j]++;
                }
            }
            int best_opt = 0; double max_v = -1;
            for (int o = 0; o < no_optima; ++o) { if (votes[o] > max_v) { max_v = votes[o]; best_opt = o; } }
            mutation_preferences(k, best_opt) += 1.0;
        }
    }
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < num_muts; ++k) { for (int o = 0; o < no_optima; ++o) mutation_preferences(k, o) /= num_sampled; }
    return mutation_preferences;
}

// [[Rcpp::export]]
NumericVector get_snv_assignment_ccfs_cpp(NumericVector pi_h_flat, IntegerVector pi_h_dims, IntegerMatrix S_i, int no_iters_burn_in, int num_threads = -1) {
    int no_iters = pi_h_dims[0], C_total = pi_h_dims[1], no_timepoints = pi_h_dims[2], no_muts = S_i.ncol(), no_iters_post_burnin = no_iters - no_iters_burn_in;
    NumericVector snv_ccfs(no_iters_post_burnin * no_muts * no_timepoints);
    omp_thread_count_cpp(num_threads);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif
    for (int t = 0; t < no_timepoints; ++t) {
        for (int i = 0; i < no_muts; ++i) {
            for (int j = 0; j < no_iters_post_burnin; ++j) {
                int it_idx = no_iters_burn_in + j;
                int cluster_id = S_i(it_idx, i);
                if (cluster_id <= 0 || cluster_id > C_total) { snv_ccfs[j + no_iters_post_burnin * i + no_iters_post_burnin * no_muts * t] = R_NaN; continue; }
                snv_ccfs[j + no_iters_post_burnin * i + no_iters_post_burnin * no_muts * t] = pi_h_flat[static_cast<size_t>(it_idx) + no_iters * (cluster_id - 1) + no_iters * C_total * t];
            }
        }
    }
    snv_ccfs.attr("dim") = IntegerVector::create(no_iters_post_burnin, no_muts, no_timepoints);
    return snv_ccfs;
}
