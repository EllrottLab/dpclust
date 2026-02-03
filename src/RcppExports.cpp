#include <Rcpp.h>
using namespace Rcpp;

// subclone_dirichlet_gibbs_cpp
List subclone_dirichlet_gibbs_cpp(NumericMatrix mutCount, NumericMatrix WTCount, NumericMatrix totalCopyNumber, NumericMatrix normalCopyNumber, NumericMatrix copyNumberAdjustment, int C, NumericVector cellularity, int iter, double conc_param, double cluster_conc);
RcppExport SEXP _DPClust_subclone_dirichlet_gibbs_cpp(SEXP mutCountSEXP, SEXP WTCountSEXP, SEXP totalCopyNumberSEXP, SEXP normalCopyNumberSEXP, SEXP copyNumberAdjustmentSEXP, SEXP CSEXP, SEXP cellularitySEXP, SEXP iterSEXP, SEXP conc_paramSEXP, SEXP cluster_concSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mutCount(mutCountSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type WTCount(WTCountSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type totalCopyNumber(totalCopyNumberSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type normalCopyNumber(normalCopyNumberSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type copyNumberAdjustment(copyNumberAdjustmentSEXP);
    Rcpp::traits::input_parameter< int >::type C(CSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cellularity(cellularitySEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type conc_param(conc_paramSEXP);
    Rcpp::traits::input_parameter< double >::type cluster_conc(cluster_concSEXP);
    rcpp_result_gen = Rcpp::wrap(subclone_dirichlet_gibbs_cpp(mutCount, WTCount, totalCopyNumber, normalCopyNumber, copyNumberAdjustment, C, cellularity, iter, conc_param, cluster_conc));
    return rcpp_result_gen;
END_RCPP
}

// assign_mutations_1d_cpp
NumericMatrix assign_mutations_1d_cpp(IntegerMatrix S_i, NumericMatrix pi_h, NumericVector boundary, IntegerVector sampledIters);
RcppExport SEXP _DPClust_assign_mutations_1d_cpp(SEXP S_iSEXP, SEXP pi_hSEXP, SEXP boundarySEXP, SEXP sampledItersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type S_i(S_iSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pi_h(pi_hSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sampledIters(sampledItersSEXP);
    rcpp_result_gen = Rcpp::wrap(assign_mutations_1d_cpp(S_i, pi_h, boundary, sampledIters));
    return rcpp_result_gen;
END_RCPP
}

// assign_mutations_nd_cpp
NumericMatrix assign_mutations_nd_cpp(IntegerMatrix S_i, NumericVector pi_h_flat, IntegerVector pi_h_dims, NumericMatrix boundary, NumericVector plane_vector_flat, NumericMatrix vector_length, LogicalMatrix vector_direction, IntegerVector sampledIters);
RcppExport SEXP _DPClust_assign_mutations_nd_cpp(SEXP S_iSEXP, SEXP pi_h_flatSEXP, SEXP pi_h_dimsSEXP, SEXP boundarySEXP, SEXP plane_vector_flatSEXP, SEXP vector_lengthSEXP, SEXP vector_directionSEXP, SEXP sampledItersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type S_i(S_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi_h_flat(pi_h_flatSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pi_h_dims(pi_h_dimsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type plane_vector_flat(plane_vector_flatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type vector_length(vector_lengthSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type vector_direction(vector_directionSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sampledIters(sampledItersSEXP);
    rcpp_result_gen = Rcpp::wrap(assign_mutations_nd_cpp(S_i, pi_h_flat, pi_h_dims, boundary, plane_vector_flat, vector_length, vector_direction, sampledIters));
    return rcpp_result_gen;
END_RCPP
}

// get_snv_assignment_ccfs_cpp
NumericVector get_snv_assignment_ccfs_cpp(NumericVector pi_h_flat, IntegerVector pi_h_dims, IntegerMatrix S_i, int no_iters_burn_in);
RcppExport SEXP _DPClust_get_snv_assignment_ccfs_cpp(SEXP pi_h_flatSEXP, SEXP pi_h_dimsSEXP, SEXP S_iSEXP, SEXP no_iters_burn_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pi_h_flat(pi_h_flatSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pi_h_dims(pi_h_dimsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type S_i(S_iSEXP);
    Rcpp::traits::input_parameter< int >::type no_iters_burn_in(no_iters_burn_inSEXP);
    rcpp_result_gen = Rcpp::wrap(get_snv_assignment_ccfs_cpp(pi_h_flat, pi_h_dims, S_i, no_iters_burn_in));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DPClust_subclone_dirichlet_gibbs_cpp", (DL_FUNC) &_DPClust_subclone_dirichlet_gibbs_cpp, 10},
    {"_DPClust_assign_mutations_1d_cpp", (DL_FUNC) &_DPClust_assign_mutations_1d_cpp, 4},
    {"_DPClust_assign_mutations_nd_cpp", (DL_FUNC) &_DPClust_assign_mutations_nd_cpp, 8},
    {"_DPClust_get_snv_assignment_ccfs_cpp", (DL_FUNC) &_DPClust_get_snv_assignment_ccfs_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_DPClust(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
