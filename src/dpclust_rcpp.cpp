#include <Rcpp.h>
#include <algorithm>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Helper functions for burden conversion
double mutationBurdenToMutationCopyNumber(double burden, double totalCopyNumber, double cellularity, double normalCopyNumber) {
    if (std::abs(cellularity) < 1e-9) return 0.0; // Avoid division by zero
    double denom = cellularity * totalCopyNumber + normalCopyNumber * (1.0 - cellularity);
    double mutCopyNumber = burden / cellularity * denom;
    if (R_IsNaN(mutCopyNumber)) return 0.0;
    return mutCopyNumber;
}

double mutationCopyNumberToMutationBurden(double copyNumber, double totalCopyNumber, double cellularity, double normalCopyNumber) {
    double denom = cellularity * totalCopyNumber + normalCopyNumber * (1.0 - cellularity);
    if (std::abs(denom) < 1e-9) return 0.000001; 
    double burden = copyNumber * cellularity / denom;
    if (R_IsNaN(burden) || burden < 0.000001) return 0.000001;
    if (burden > 0.999999) return 0.999999;
    return burden;
}

// [[Rcpp::export]]
List subclone_dirichlet_gibbs_cpp(NumericMatrix mutCount, NumericMatrix WTCount, 
                                  NumericMatrix totalCopyNumber, NumericMatrix normalCopyNumber, 
                                  NumericMatrix copyNumberAdjustment, 
                                  int C, NumericVector cellularity, int iter, 
                                  double conc_param, double cluster_conc,
                                  bool keep_aux_fields,
                                  int num_threads,
                                  IntegerVector stored_iters,
                                  Function log_func = R_NilValue) {
    
    int num_muts = mutCount.nrow();
    int num_timepoints = mutCount.ncol();
    int active_threads = 1;
#ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
        active_threads = num_threads;
    } else {
        active_threads = omp_get_max_threads();
    }
#endif
    
    double A = 1.0;
    double B = conc_param;
    
    // Output structures
    // pi.h: iter x C x num_timepoints
    // Dimensions: [iter * C * num_timepoints]
    NumericVector pi_h(iter * C * num_timepoints);
    // Index (m, c, t) (0-based) = m + iter * c + iter * C * t
    
    NumericMatrix V_h(iter, C);
    std::fill(V_h.begin(), V_h.end(), 1.0);
    
    std::vector<int> keep_iters_zero_based;
    keep_iters_zero_based.reserve(stored_iters.size());
    if (stored_iters.size() > 0) {
        std::vector<bool> seen(iter, false);
        for (int i = 0; i < stored_iters.size(); ++i) {
            int iter_idx = stored_iters[i] - 1; // 1-based from R
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
        for (int i = 0; i < iter; ++i) {
            iter_to_store_index[i] = i;
        }
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
    
    // Pre-calculate denominator for burden conversion (constant across iterations)
    NumericMatrix burden_denom_inv(num_muts, num_timepoints);
    for (int t = 0; t < num_timepoints; ++t) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int k = 0; k < num_muts; ++k) {
            double denom = cellularity[t] * totalCopyNumber(k, t) + normalCopyNumber(k, t) * (1.0 - cellularity[t]);
            if (std::abs(denom) < 1e-9) {
                burden_denom_inv(k, t) = 0.0;
            } else {
                burden_denom_inv(k, t) = cellularity[t] / denom;
            }
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
            double burden = mutCount(k, t) / (mutCount(k, t) + WTCount(k, t));
            if (R_IsNaN(burden)) burden = 0; 
            
            double mcn = mutationBurdenToMutationCopyNumber(burden, totalCopyNumber(k, t), cellularity[t], normalCopyNumber(k, t));
            mcn /= copyNumberAdjustment(k, t);
            
            if (mcn < min_val) min_val = mcn;
            if (mcn > max_val) max_val = mcn;
        }
        
        lower[t] = min_val;
        upper[t] = max_val;
        double diff = upper[t] - lower[t];
        lower[t] -= diff / 10.0;
        upper[t] += diff / 10.0;
        
        // Randomise starting positions (iteration 0 - which is 1 in R)
        for (int c = 0; c < C; ++c) {
            double val = R::runif(lower[t], upper[t]);
            // pi.h[1, c, t] -> m=0
            pi_h[0 + iter * c + iter * C * t] = val;
        }
    }
    
    // V.h[1, ] init
    for (int c = 0; c < C - 1; ++c) V_h(0, c) = 0.5;
    V_h(0, C - 1) = 1.0;
    
    // S.i[1, ] init: all mutations in cluster 1.
    std::vector<int> S_curr(num_muts, 1);
    if (iter_to_store_index[0] >= 0) {
        int row = iter_to_store_index[0];
        for (int k = 0; k < num_muts; ++k) {
            S_i(row, k) = S_curr[k];
        }
    }
    
    alpha[0] = 1.0;
    
    // Pre-allocate thread-local buffers for allocation step
    std::vector<double> Pr_S_threads(static_cast<size_t>(active_threads) * C);
    std::vector<double> sampled_uniforms(num_muts, 0.0);

    // Build stick-breaking priors once per iteration and reuse for each mutation.
    std::vector<double> log_prior(C, 0.0);

    // Pre-allocate thread-local buffers for shape/rate updates
    int ct = C * num_timepoints;
    std::vector<double> shape_sums_thread(static_cast<size_t>(active_threads) * ct, 0.0);
    std::vector<double> rate_sums_thread(static_cast<size_t>(active_threads) * ct, 0.0);
    std::vector<double> shape_sums(ct, 0.0);
    std::vector<double> rate_sums(ct, 0.0);

    // MCMC Loop
    for (int m = 1; m < iter; ++m) {
        if ((m + 1) % 100 == 0) {
            if (log_func != R_NilValue) {
                log_func("Iteration " + std::to_string(m + 1) + " / " + std::to_string(iter));
            } else {
                Rcout << "Iteration " << m + 1 << " / " << iter << std::endl;
            }
            Rcpp::checkUserInterrupt(); 
        }
        
        // Build stick-breaking priors once per iteration and reuse for each mutation.
        std::vector<double> log_prior(C, 0.0);
        log_prior[0] = std::log(V_h(m - 1, 0));
        double sum_log_1_minus_V = 0.0;
        for (int j = 1; j < C; ++j) {
            sum_log_1_minus_V += std::log(1.0 - V_h(m - 1, j - 1));
            log_prior[j] = std::log(V_h(m - 1, j)) + sum_log_1_minus_V;
        }

        std::vector<double> sampled_uniforms(num_muts, 0.0);
        for (int k = 0; k < num_muts; ++k) {
            sampled_uniforms[k] = R::runif(0.0, 1.0);
        }

#ifdef _OPENMP
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            double* Pr_S = &Pr_S_threads[static_cast<size_t>(tid) * C];
#pragma omp for schedule(static)
#else
        std::vector<double> Pr_S_vec(C);
        double* Pr_S = Pr_S_vec.data();
#endif
        for (int k = 0; k < num_muts; ++k) {
            for (int c = 0; c < C; ++c) {
                Pr_S[c] = log_prior[c];
            }
            
            // Add Likelihood
            for (int t = 0; t < num_timepoints; ++t) {
                for (int c = 0; c < C; ++c) {
                    double val = pi_h[(m - 1) + iter * c + iter * C * t];
                    double mb = val * copyNumberAdjustment(k, t) * burden_denom_inv(k, t);
                    if (R_IsNaN(mb) || mb < 0.000001) mb = 0.000001;
                    if (mb > 0.999999) mb = 0.999999;
                    Pr_S[c] += mutCount(k, t) * std::log(mb) + WTCount(k, t) * std::log(1.0 - mb);
                }
            }
            
            // Normalize in log space then exp
            double max_val = Pr_S[0];
            for(int c=1; c<C; ++c) if(Pr_S[c] > max_val) max_val = Pr_S[c];
            
            double sum_exp = 0.0;
            for(int c=0; c<C; ++c) {
                Pr_S[c] -= max_val;
                if (R_IsNaN(Pr_S[c])) Pr_S[c] = -700; // Small log prob
                double ex = std::exp(Pr_S[c]);
                Pr_S[c] = ex;
                sum_exp += ex;
            }
            
            if (sum_exp <= 0.0 || !std::isfinite(sum_exp)) {
                sum_exp = static_cast<double>(C);
                for (int c = 0; c < C; ++c) {
                    Pr_S[c] = 1.0;
                }
            }
            for(int c=0; c<C; ++c) Pr_S[c] /= sum_exp;
            
            // Multinomial sampling
            double r = sampled_uniforms[k];
            double cum_sum = 0.0;
            int picked = C - 1;
            for(int c = 0; c < C; ++c) {
                cum_sum += Pr_S[c];
                if (r <= cum_sum) {
                    picked = c;
                    break;
                }
            }
            S_curr[k] = picked + 1; // 1-based
        }
#ifdef _OPENMP
        } // end parallel
#endif
        
        std::vector<int> cluster_counts(C, 0);
        for (int k = 0; k < num_muts; ++k) {
            int c_idx = S_curr[k] - 1;
            if (c_idx >= 0 && c_idx < C) {
                cluster_counts[c_idx]++;
            }
        }

        // Update stick-breaking weights
        int cumulative_count = 0;
        for (int c = 0; c < C - 1; ++c) {
            double count_eq = static_cast<double>(cluster_counts[c]);
            cumulative_count += cluster_counts[c];
            double count_gt = static_cast<double>(num_muts - cumulative_count);
            
            V_h(m, c) = R::rbeta(1.0 + count_eq, alpha[m-1] + count_gt);
            if(V_h(m, c) == 1.0) V_h(m, c) = 0.999;
        }
        V_h(m, C - 1) = 1.0;
        
       
        // Randomise unused pi.h
        for(int t=0; t<num_timepoints; ++t) {
             for(int c=0; c<C; ++c) {
                 pi_h[m + iter * c + iter * C * t] = R::runif(lower[t], upper[t]);
             }
        }
        
        // Update populated pi.h using one-pass cluster/timepoint aggregates.
        std::fill(shape_sums_thread.begin(), shape_sums_thread.end(), 0.0);
        std::fill(rate_sums_thread.begin(), rate_sums_thread.end(), 0.0);
        std::fill(shape_sums.begin(), shape_sums.end(), 0.0);
        std::fill(rate_sums.begin(), rate_sums.end(), 0.0);

#ifdef _OPENMP
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            double* shape_local = &shape_sums_thread[static_cast<size_t>(tid) * ct];
            double* rate_local = &rate_sums_thread[static_cast<size_t>(tid) * ct];
#pragma omp for schedule(static)
            for (int k = 0; k < num_muts; ++k) {
                int c = S_curr[k] - 1;
                if (c < 0 || c >= C) continue;
                for (int t = 0; t < num_timepoints; ++t) {
                    int idx = c + C * t;
                    shape_local[idx] += mutCount(k, t);
                    double mb_unit = copyNumberAdjustment(k, t) * burden_denom_inv(k, t);
                    if (R_IsNaN(mb_unit) || mb_unit < 0.000001) mb_unit = 0.000001;
                    if (mb_unit > 0.999999) mb_unit = 0.999999;
                    rate_local[idx] += (mutCount(k, t) + WTCount(k, t)) * mb_unit;
                }
            }
        }
        for (int tid = 0; tid < active_threads; ++tid) {
            const double* shape_local = &shape_sums_thread[static_cast<size_t>(tid) * ct];
            const double* rate_local = &rate_sums_thread[static_cast<size_t>(tid) * ct];
            for (int idx = 0; idx < ct; ++idx) {
                shape_sums[idx] += shape_local[idx];
                rate_sums[idx] += rate_local[idx];
            }
        }
#else
        for (int k = 0; k < num_muts; ++k) {
            int c = S_curr[k] - 1;
            if (c < 0 || c >= C) continue;
            for (int t = 0; t < num_timepoints; ++t) {
                int idx = c + C * t;
                shape_sums[idx] += mutCount(k, t);
                double mb_unit = copyNumberAdjustment(k, t) * burden_denom_inv(k, t);
                if (R_IsNaN(mb_unit) || mb_unit < 0.000001) mb_unit = 0.000001;
                if (mb_unit > 0.999999) mb_unit = 0.999999;
                rate_sums[idx] += (mutCount(k, t) + WTCount(k, t)) * mb_unit;
            }
        }
#endif
        for (int c = 0; c < C; ++c) {
            if (cluster_counts[c] == 0) continue;
            for (int t = 0; t < num_timepoints; ++t) {
                int idx = c + C * t;
                double shape_sum = shape_sums[idx];
                double rate_sum = rate_sums[idx];
                if (rate_sum == 0.0) {
                    pi_h[m + iter * c + iter * C * t] = 0.0;
                } else {
                    pi_h[m + iter * c + iter * C * t] = R::rgamma(shape_sum, 1.0 / rate_sum);
                }
            }
        }
        if (iter_to_store_index[m] >= 0) {
            int row = iter_to_store_index[m];
            for (int k = 0; k < num_muts; ++k) {
                S_i(row, k) = S_curr[k];
            }
        }
        // Update alpha
        double sum_log = 0.0;
        for(int c=0; c<C-1; ++c) sum_log += std::log(1.0 - V_h(m, c));
        alpha[m] = R::rgamma(C + A - 1, 1.0 / (B - sum_log));
    }
    
    // Set dimensions
    pi_h.attr("dim") = Dimension(iter, C, num_timepoints);
    SEXP stored_iters_sexp = R_NilValue;
    if (!store_all_iters) {
        stored_iters_sexp = Rcpp::wrap(stored_iters_out);
    }
    
    SEXP y1_out = keep_aux_fields ? Rcpp::wrap(mutCount) : R_NilValue;
    SEXP n1_out = keep_aux_fields ? Rcpp::wrap(mutCount + WTCount) : R_NilValue;

    return List::create(Named("S.i") = S_i, 
                        Named("V.h") = V_h, 
                        Named("pi.h") = pi_h, 
                        Named("stored_iters") = stored_iters_sexp,
                        Named("mutBurdens") = R_NilValue,
                        Named("alpha") = alpha, 
                        Named("y1") = y1_out, 
                        Named("N1") = n1_out);
}

// [[Rcpp::export]]
NumericMatrix assign_mutations_1d_cpp(IntegerMatrix S_i, NumericMatrix pi_h, NumericVector boundary, IntegerVector sampledIters_pi, IntegerVector sampledIters_state) {
    if (sampledIters_pi.size() != sampledIters_state.size()) {
        stop("sampledIters_pi and sampledIters_state must have equal length.");
    }
    int num_muts = S_i.ncol();
    int num_optima = boundary.size() + 1;
    int num_sampled = sampledIters_pi.size();
    NumericMatrix mutation_preferences(num_muts, num_optima);

    for (int s_idx = 0; s_idx < num_sampled; ++s_idx) {
        int s_pi = sampledIters_pi[s_idx] - 1; // 0-based
        int s_state = sampledIters_state[s_idx] - 1; // 0-based
        if (s_pi < 0 || s_pi >= pi_h.nrow()) {
            stop("sampledIters_pi contains out-of-range index.");
        }
        if (s_state < 0 || s_state >= S_i.nrow()) {
            stop("sampledIters_state contains out-of-range index.");
        }
        
        int max_c = 0;
        for(int k=0; k<num_muts; ++k) if(S_i(s_state, k) > max_c) max_c = S_i(s_state, k);
        
        std::vector<int> c_to_opt(max_c + 1, 0);
        for(int k=0; k<num_muts; ++k) {
            int c = S_i(s_state, k);
            if (c <= 0) continue;
            if(c_to_opt[c] == 0) {
                double val = pi_h(s_pi, c-1); 
                int opt = 0;
                for(int b=0; b<boundary.size(); ++b) {
                    if(val > boundary[b]) opt++;
                }
                c_to_opt[c] = opt + 1;
            }
            mutation_preferences(k, c_to_opt[c]-1) += 1.0;
        }
    }

    for(int k=0; k<num_muts; ++k) {
        for(int o=0; o<num_optima; ++o) {
            mutation_preferences(k, o) /= num_sampled;
        }
    }
    
    return mutation_preferences;
}

// [[Rcpp::export]]
NumericMatrix assign_mutations_nd_cpp(IntegerMatrix S_i, NumericVector pi_h_flat, IntegerVector pi_h_dims, NumericMatrix boundary, NumericVector plane_vector_flat, NumericMatrix vector_length, LogicalMatrix vector_direction, IntegerVector sampledIters_pi, IntegerVector sampledIters_state) {
    if (sampledIters_pi.size() != sampledIters_state.size()) {
        stop("sampledIters_pi and sampledIters_state must have equal length.");
    }
    int num_muts = S_i.ncol();
    int no_iters = pi_h_dims[0];
    int C_total = pi_h_dims[1];
    int no_subsamples = pi_h_dims[2];
    int no_optima = vector_length.nrow(); // boundary might be smaller in some edge cases but length is safe
    int num_sampled = sampledIters_pi.size();
    
    NumericMatrix mutation_preferences(num_muts, no_optima);

    for (int s_idx = 0; s_idx < num_sampled; ++s_idx) {
        int s_pi = sampledIters_pi[s_idx] - 1;
        int s_state = sampledIters_state[s_idx] - 1;
        if (s_pi < 0 || s_pi >= no_iters) {
            stop("sampledIters_pi contains out-of-range index.");
        }
        if (s_state < 0 || s_state >= S_i.nrow()) {
            stop("sampledIters_state contains out-of-range index.");
        }
        
        int max_c = 0;
        for(int k=0; k<num_muts; ++k) if(S_i(s_state, k) > max_c) max_c = S_i(s_state, k);
        std::vector<int> c_to_opt(max_c + 1, -1);
        
        for(int k=0; k<num_muts; ++k) {
            int c_1based = S_i(s_state, k);
            if (c_1based <= 0) continue;
            int c = c_1based - 1;
            if(c_to_opt[c_1based] == -1) {
                std::vector<double> votes(no_optima, 0.0);
                for(int i=0; i<no_optima-1; ++i) {
                    for(int j=i+1; j<no_optima; ++j) {
                        double distance_from_plane = 0;
                        for(int t=0; t<no_subsamples; ++t) {
                            // plane_vector_flat is [no_optima, no_optima, no_subsamples+1]
                            // Index: i + no_optima * j + no_optima * no_optima * t
                            distance_from_plane += pi_h_flat[s_pi + no_iters * c + no_iters * C_total * t] * plane_vector_flat[i + no_optima * j + no_optima * no_optima * t];
                        }
                        distance_from_plane /= vector_length(i, j);
                        
                        bool lead_to_i;
                        if(distance_from_plane <= boundary(i, j)) {
                            lead_to_i = vector_direction(i, j);
                        } else {
                            lead_to_i = !vector_direction(i, j);
                        }
                        
                        if(lead_to_i) votes[i]++; else votes[j]++;
                    }
                }
                double max_votes = -1;
                int best_opt = 0;
                for(int o=0; o<no_optima; ++o) {
                    if(votes[o] > max_votes) {
                        max_votes = votes[o];
                        best_opt = o;
                    }
                }
                c_to_opt[c_1based] = best_opt;
            }
            mutation_preferences(k, c_to_opt[c_1based]) += 1.0;
        }
    }

    for(int k=0; k<num_muts; ++k) {
        for(int o=0; o<no_optima; ++o) {
            mutation_preferences(k, o) /= num_sampled;
        }
    }
    
    return mutation_preferences;
}

// [[Rcpp::export]]
NumericVector get_snv_assignment_ccfs_cpp(NumericVector pi_h_flat, IntegerVector pi_h_dims, IntegerMatrix S_i, int no_iters_burn_in) {
    int no_iters = pi_h_dims[0];
    int C_total = pi_h_dims[1];
    int no_timepoints = pi_h_dims[2];
    int no_muts = S_i.ncol();
    int no_iters_post_burnin = no_iters - no_iters_burn_in;
    
    NumericVector snv_ccfs(no_iters_post_burnin * no_muts * no_timepoints);
    
    for (int t = 0; t < no_timepoints; ++t) {
        for (int i = 0; i < no_muts; ++i) {
            for (int j = 0; j < no_iters_post_burnin; ++j) {
                int it_idx = no_iters_burn_in + j;
                int cluster_id = S_i(it_idx, i); // 1-based
                if (cluster_id <= 0 || cluster_id > C_total) {
                    snv_ccfs[j + no_iters_post_burnin * i + no_iters_post_burnin * no_muts * t] = R_NaN;
                    continue;
                }
                double val = pi_h_flat[it_idx + no_iters * (cluster_id - 1) + no_iters * C_total * t];
                snv_ccfs[j + no_iters_post_burnin * i + no_iters_post_burnin * no_muts * t] = val;
            }
        }
    }
    
    snv_ccfs.attr("dim") = IntegerVector::create(no_iters_post_burnin, no_muts, no_timepoints);
    return snv_ccfs;
}
