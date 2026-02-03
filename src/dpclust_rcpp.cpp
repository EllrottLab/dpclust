#include <Rcpp.h>
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
                                  double conc_param, double cluster_conc) {
    
    int num_muts = mutCount.nrow();
    int num_timepoints = mutCount.ncol();
    
    double A = 1.0;
    double B = conc_param;
    
    // Output structures
    // pi.h: iter x C x num_timepoints
    // Dimensions: [iter * C * num_timepoints]
    NumericVector pi_h(iter * C * num_timepoints);
    // Index (m, c, t) (0-based) = m + iter * c + iter * C * t
    
    // mutBurdens: C x num_timepoints x num_muts
    // dim = c(C, num_timepoints, num_muts)
    NumericVector mutBurdens(C * num_timepoints * num_muts);
    // Index (c, t, k) = c + C * t + C * num_timepoints * k
    
    NumericMatrix V_h(iter, C);
    std::fill(V_h.begin(), V_h.end(), 1.0);
    
    IntegerMatrix S_i(iter, num_muts);
    NumericMatrix Pr_S(num_muts, C);
    NumericVector alpha(iter);
    
    NumericVector lower(num_timepoints);
    NumericVector upper(num_timepoints);
    
    NumericMatrix mutCopyNum(num_muts, num_timepoints);
    
    // Pre-calculate denominator for burden conversion (constant across iterations)
    NumericMatrix burden_denom_inv(num_muts, num_timepoints);
    for (int t = 0; t < num_timepoints; ++t) {
        for (int k = 0; k < num_muts; ++k) {
            double denom = cellularity[t] * totalCopyNumber(k, t) + normalCopyNumber(k, t) * (1.0 - cellularity[t]);
            if (std::abs(denom) < 1e-9) {
                burden_denom_inv(k, t) = 0.0;
            } else {
                burden_denom_inv(k, t) = cellularity[t] / denom;
            }
        }
    }

    Rcout << "Starting Gibbs sampler for " << num_muts << " mutations..." << std::endl;
    
    // Initialization
    for (int t = 0; t < num_timepoints; ++t) {
        double min_val = R_PosInf;
        double max_val = R_NegInf;
        
        for (int k = 0; k < num_muts; ++k) {
            double burden = mutCount(k, t) / (mutCount(k, t) + WTCount(k, t));
            if (R_IsNaN(burden)) burden = 0; 
            
            double mcn = mutationBurdenToMutationCopyNumber(burden, totalCopyNumber(k, t), cellularity[t], normalCopyNumber(k, t));
            mcn /= copyNumberAdjustment(k, t);
            mutCopyNum(k, t) = mcn;
            
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
            
            for (int k = 0; k < num_muts; ++k) {
                double burden = val * copyNumberAdjustment(k, t) * burden_denom_inv(k, t);
                if (R_IsNaN(burden) || burden < 0.000001) burden = 0.000001;
                if (burden > 0.999999) burden = 0.999999;
                mutBurdens[c + C * t + C * num_timepoints * k] = burden;
            }
        }
    }
    
    // V.h[1, ] init
    for (int c = 0; c < C - 1; ++c) V_h(0, c) = 0.5;
    V_h(0, C - 1) = 1.0;
    
    // S.i[1, ] init
    for (int k = 0; k < num_muts; ++k) {
        S_i(0, k) = (k == 0) ? 1 : 1; // R code: c(1, rep(0, ...)). 0 is invalid 1-based index usually. Setting all to 1 is safer.
        // Actually, let's replicate the logic of 1 and 0 just in case downstream code relies on 0 for "unassigned".
        if (k > 0) S_i(0, k) = 1; // Changed to 1 to differ from 0 but keep valid. R code is weird.
        // Actually, since loop starts at m=1 (2nd iter), S_i[0] matters little unless accessed.
        // I will set all to 1.
    }
    
    alpha[0] = 1.0;
    
    // MCMC Loop
    for (int m = 1; m < iter; ++m) {
        if ((m + 1) % 100 == 0) {
            Rcout << "Iteration " << m + 1 << " / " << iter << std::endl;
            Rcpp::checkUserInterrupt(); 
        }
        
        // Update cluster allocation
        for (int k = 0; k < num_muts; ++k) {
            
            // Calculate prior Pr.S
            double log_V_h_0 = std::log(V_h(m - 1, 0));
            Pr_S(k, 0) = log_V_h_0;
            
            double sum_log_1_minus_V = 0.0;
            for (int j = 1; j < C; ++j) {
                sum_log_1_minus_V += std::log(1.0 - V_h(m - 1, j - 1));
                Pr_S(k, j) = std::log(V_h(m - 1, j)) + sum_log_1_minus_V;
            }
            
            // Add Likelihood
            for (int t = 0; t < num_timepoints; ++t) {
                for (int c = 0; c < C; ++c) {
                    double mb = mutBurdens[c + C * t + C * num_timepoints * k];
                    Pr_S(k, c) += mutCount(k, t) * std::log(mb) + WTCount(k, t) * std::log(1.0 - mb);
                }
            }
            
            // Normalize in log space then exp
            double max_val = Pr_S(k, 0);
            for(int c=1; c<C; ++c) if(Pr_S(k, c) > max_val) max_val = Pr_S(k, c);
            
            double sum_exp = 0.0;
            for(int c=0; c<C; ++c) {
                Pr_S(k, c) -= max_val;
                if (R_IsNaN(Pr_S(k, c))) Pr_S(k, c) = -700; // Small log prob
                double ex = std::exp(Pr_S(k, c));
                Pr_S(k, c) = ex;
                sum_exp += ex;
            }
            
            for(int c=0; c<C; ++c) Pr_S(k, c) /= sum_exp;
            
            // Multinomial sampling
            double r = R::runif(0.0, 1.0);
            double cum_sum = 0.0;
            int picked = C - 1;
            for(int c = 0; c < C; ++c) {
                cum_sum += Pr_S(k, c);
                if (r <= cum_sum) {
                    picked = c;
                    break;
                }
            }
            S_i(m, k) = picked + 1; // 1-based
        }
        
        // Update stick-breaking weights
        for (int c = 0; c < C - 1; ++c) {
            double count_eq = 0;
            double count_gt = 0;
            int cluster_id = c + 1;
            for(int k=0; k<num_muts; ++k) {
                if(S_i(m, k) == cluster_id) count_eq++;
                if(S_i(m, k) > cluster_id) count_gt++;
            }
            
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
        
        // Update populated pi.h
        std::set<int> unique_clusters; 
        for(int k=0; k<num_muts; ++k) unique_clusters.insert(S_i(m, k));
        
        for(int c_1based : unique_clusters) {
            int c = c_1based - 1;
            for(int t=0; t<num_timepoints; ++t) {
                 double shape_sum = 0.0;
                 double rate_sum = 0.0;
                 
                 for(int k=0; k<num_muts; ++k) {
                     if(S_i(m, k) == c_1based) {
                         shape_sum += mutCount(k, t);
                         double mb_unit = copyNumberAdjustment(k, t) * burden_denom_inv(k, t);
                         if (R_IsNaN(mb_unit) || mb_unit < 0.000001) mb_unit = 0.000001;
                         if (mb_unit > 0.999999) mb_unit = 0.999999;
                         rate_sum += (mutCount(k, t) + WTCount(k, t)) * mb_unit;
                     }
                 }
                 
                 if (rate_sum == 0.0) {
                     pi_h[m + iter * c + iter * C * t] = 0.0;
                 } else {
                     pi_h[m + iter * c + iter * C * t] = R::rgamma(shape_sum, 1.0/rate_sum);
                 }
            }
        }
        
        // Update mutBurdens
        for(int t=0; t<num_timepoints; ++t) {
            for(int c=0; c<C; ++c) {
                double val = pi_h[m + iter * c + iter * C * t];
                for(int k=0; k<num_muts; ++k) {
                    double burden = val * copyNumberAdjustment(k, t) * burden_denom_inv(k, t);
                    if (R_IsNaN(burden) || burden < 0.000001) burden = 0.000001;
                    if (burden > 0.999999) burden = 0.999999;
                    mutBurdens[c + C * t + C * num_timepoints * k] = burden;
                }
            }
        }
        
        // Update alpha
        double sum_log = 0.0;
        for(int c=0; c<C-1; ++c) sum_log += std::log(1.0 - V_h(m, c));
        alpha[m] = R::rgamma(C + A - 1, 1.0 / (B - sum_log));
    }
    
    // Set dimensions
    pi_h.attr("dim") = Dimension(iter, C, num_timepoints);
    mutBurdens.attr("dim") = Dimension(C, num_timepoints, num_muts);
    
    return List::create(Named("S.i") = S_i, 
                        Named("V.h") = V_h, 
                        Named("pi.h") = pi_h, 
                        Named("mutBurdens") = mutBurdens, 
                        Named("alpha") = alpha, 
                        Named("y1") = mutCount, 
                        Named("N1") = mutCount + WTCount);
}

// [[Rcpp::export]]
NumericMatrix assign_mutations_1d_cpp(IntegerMatrix S_i, NumericMatrix pi_h, NumericVector boundary, IntegerVector sampledIters) {
    int num_muts = S_i.ncol();
    int num_optima = boundary.size() + 1;
    int num_sampled = sampledIters.size();
    NumericMatrix mutation_preferences(num_muts, num_optima);

    for (int s_idx = 0; s_idx < num_sampled; ++s_idx) {
        int s = sampledIters[s_idx] - 1; // 0-based
        NumericMatrix temp_preferences(num_muts, num_optima);
        
        int max_c = 0;
        for(int k=0; k<num_muts; ++k) if(S_i(s, k) > max_c) max_c = S_i(s, k);
        
        std::vector<int> c_to_opt(max_c + 1, 0);
        for(int k=0; k<num_muts; ++k) {
            int c = S_i(s, k);
            if(c_to_opt[c] == 0) {
                double val = pi_h(s, c-1); 
                int opt = 0;
                for(int b=0; b<boundary.size(); ++b) {
                    if(val > boundary[b]) opt++;
                }
                c_to_opt[c] = opt + 1;
            }
            temp_preferences(k, c_to_opt[c]-1)++;
        }
        
        for(int k=0; k<num_muts; ++k) {
            double max_val = 0;
            int count_max = 0;
            for(int o=0; o<num_optima; ++o) {
                if(temp_preferences(k, o) > max_val) {
                    max_val = temp_preferences(k, o);
                    count_max = 1;
                } else if(temp_preferences(k, o) == max_val) {
                    count_max++;
                }
            }
            if(max_val > 0) {
                for(int o=0; o<num_optima; ++o) {
                    if(temp_preferences(k, o) == max_val) {
                        mutation_preferences(k, o) += 1.0 / count_max;
                    }
                }
            }
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
NumericMatrix assign_mutations_nd_cpp(IntegerMatrix S_i, NumericVector pi_h_flat, IntegerVector pi_h_dims, NumericMatrix boundary, NumericVector plane_vector_flat, NumericMatrix vector_length, LogicalMatrix vector_direction, IntegerVector sampledIters) {
    int num_muts = S_i.ncol();
    int no_iters = pi_h_dims[0];
    int C_total = pi_h_dims[1];
    int no_subsamples = pi_h_dims[2];
    int no_optima = vector_length.nrow(); // boundary might be smaller in some edge cases but length is safe
    int num_sampled = sampledIters.size();
    
    NumericMatrix mutation_preferences(num_muts, no_optima);

    for (int s_idx = 0; s_idx < num_sampled; ++s_idx) {
        int s = sampledIters[s_idx] - 1;
        NumericMatrix temp_preferences(num_muts, no_optima);
        
        int max_c = 0;
        for(int k=0; k<num_muts; ++k) if(S_i(s, k) > max_c) max_c = S_i(s, k);
        std::vector<int> c_to_opt(max_c + 1, -1);
        
        for(int k=0; k<num_muts; ++k) {
            int c_1based = S_i(s, k);
            int c = c_1based - 1;
            if(c_to_opt[c_1based] == -1) {
                std::vector<double> votes(no_optima, 0.0);
                for(int i=0; i<no_optima-1; ++i) {
                    for(int j=i+1; j<no_optima; ++j) {
                        double distance_from_plane = 0;
                        for(int t=0; t<no_subsamples; ++t) {
                            // plane_vector_flat is [no_optima, no_optima, no_subsamples+1]
                            // Index: i + no_optima * j + no_optima * no_optima * t
                            distance_from_plane += pi_h_flat[s + no_iters * c + no_iters * C_total * t] * plane_vector_flat[i + no_optima * j + no_optima * no_optima * t];
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
            temp_preferences(k, c_to_opt[c_1based])++;
        }
        
        for(int k=0; k<num_muts; ++k) {
            double max_val = 0;
            int count_max = 0;
            for(int o=0; o<no_optima; ++o) {
                if(temp_preferences(k, o) > max_val) {
                    max_val = temp_preferences(k, o);
                    count_max = 1;
                } else if(temp_preferences(k, o) == max_val) {
                    count_max++;
                }
            }
            if(max_val > 0) {
                for(int o=0; o<no_optima; ++o) {
                    if(temp_preferences(k, o) == max_val) {
                        mutation_preferences(k, o) += 1.0 / count_max;
                    }
                }
            }
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
