// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// tableC
arma::vec tableC(arma::vec x);
RcppExport SEXP PPtreeExt_tableC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(tableC(x));
    return rcpp_result_gen;
END_RCPP
}
// roundme
double roundme(double x);
RcppExport SEXP PPtreeExt_roundme(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(roundme(x));
    return rcpp_result_gen;
END_RCPP
}
// LDAindex2
arma::vec LDAindex2(arma::vec origclass, arma::mat origdata, arma::mat proj, bool weight);
RcppExport SEXP PPtreeExt_LDAindex2(SEXP origclassSEXP, SEXP origdataSEXP, SEXP projSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type origclass(origclassSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type origdata(origdataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type proj(projSEXP);
    Rcpp::traits::input_parameter< bool >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(LDAindex2(origclass, origdata, proj, weight));
    return rcpp_result_gen;
END_RCPP
}
// signC
double signC(double x);
RcppExport SEXP PPtreeExt_signC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(signC(x));
    return rcpp_result_gen;
END_RCPP
}
// LDAopt
arma::vec LDAopt(arma::vec origclass, arma::mat origdata, int q, std::string PPmethod, bool weight);
RcppExport SEXP PPtreeExt_LDAopt(SEXP origclassSEXP, SEXP origdataSEXP, SEXP qSEXP, SEXP PPmethodSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type origclass(origclassSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type origdata(origdataSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< std::string >::type PPmethod(PPmethodSEXP);
    Rcpp::traits::input_parameter< bool >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(LDAopt(origclass, origdata, q, PPmethod, weight));
    return rcpp_result_gen;
END_RCPP
}
// PDAindex2
double PDAindex2(arma::vec origclass, arma::mat origdata, arma::mat proj, bool weight, double lambda);
RcppExport SEXP PPtreeExt_PDAindex2(SEXP origclassSEXP, SEXP origdataSEXP, SEXP projSEXP, SEXP weightSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type origclass(origclassSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type origdata(origdataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type proj(projSEXP);
    Rcpp::traits::input_parameter< bool >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(PDAindex2(origclass, origdata, proj, weight, lambda));
    return rcpp_result_gen;
END_RCPP
}
// PDAopt
arma::vec PDAopt(arma::vec origclass, arma::mat origdata, int q, std::string PPmethod, bool weight, double lambda);
RcppExport SEXP PPtreeExt_PDAopt(SEXP origclassSEXP, SEXP origdataSEXP, SEXP qSEXP, SEXP PPmethodSEXP, SEXP weightSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type origclass(origclassSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type origdata(origdataSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< std::string >::type PPmethod(PPmethodSEXP);
    Rcpp::traits::input_parameter< bool >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(PDAopt(origclass, origdata, q, PPmethod, weight, lambda));
    return rcpp_result_gen;
END_RCPP
}
// varselect
arma::uvec varselect(int p, int s);
RcppExport SEXP PPtreeExt_varselect(SEXP pSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(varselect(p, s));
    return rcpp_result_gen;
END_RCPP
}
// datanode
List datanode(arma::mat origdata, double sizep);
RcppExport SEXP PPtreeExt_datanode(SEXP origdataSEXP, SEXP sizepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type origdata(origdataSEXP);
    Rcpp::traits::input_parameter< double >::type sizep(sizepSEXP);
    rcpp_result_gen = Rcpp::wrap(datanode(origdata, sizep));
    return rcpp_result_gen;
END_RCPP
}
// entropy
double entropy(arma::vec origclass);
RcppExport SEXP PPtreeExt_entropy(SEXP origclassSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type origclass(origclassSEXP);
    rcpp_result_gen = Rcpp::wrap(entropy(origclass));
    return rcpp_result_gen;
END_RCPP
}
// split_entro
double split_entro(arma::vec origclass, arma::colvec projdata);
RcppExport SEXP PPtreeExt_split_entro(SEXP origclassSEXP, SEXP projdataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type origclass(origclassSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type projdata(projdataSEXP);
    rcpp_result_gen = Rcpp::wrap(split_entro(origclass, projdata));
    return rcpp_result_gen;
END_RCPP
}
// split_relMOD
List split_relMOD(arma::vec origclass, arma::colvec projdata, bool entro, bool entroindiv);
RcppExport SEXP PPtreeExt_split_relMOD(SEXP origclassSEXP, SEXP projdataSEXP, SEXP entroSEXP, SEXP entroindivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type origclass(origclassSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type projdata(projdataSEXP);
    Rcpp::traits::input_parameter< bool >::type entro(entroSEXP);
    Rcpp::traits::input_parameter< bool >::type entroindiv(entroindivSEXP);
    rcpp_result_gen = Rcpp::wrap(split_relMOD(origclass, projdata, entro, entroindiv));
    return rcpp_result_gen;
END_RCPP
}
// findprojMOD
List findprojMOD(arma::vec origclass, arma::mat origdata, std::string PPmethod, double lambda, bool entro, bool entroindiv);
RcppExport SEXP PPtreeExt_findprojMOD(SEXP origclassSEXP, SEXP origdataSEXP, SEXP PPmethodSEXP, SEXP lambdaSEXP, SEXP entroSEXP, SEXP entroindivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type origclass(origclassSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type origdata(origdataSEXP);
    Rcpp::traits::input_parameter< std::string >::type PPmethod(PPmethodSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type entro(entroSEXP);
    Rcpp::traits::input_parameter< bool >::type entroindiv(entroindivSEXP);
    rcpp_result_gen = Rcpp::wrap(findprojMOD(origclass, origdata, PPmethod, lambda, entro, entroindiv));
    return rcpp_result_gen;
END_RCPP
}
// findproj1D
List findproj1D(arma::vec origclass, arma::mat origdata, std::string PPmethod, double lambda, bool entro, bool entroindiv);
RcppExport SEXP PPtreeExt_findproj1D(SEXP origclassSEXP, SEXP origdataSEXP, SEXP PPmethodSEXP, SEXP lambdaSEXP, SEXP entroSEXP, SEXP entroindivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type origclass(origclassSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type origdata(origdataSEXP);
    Rcpp::traits::input_parameter< std::string >::type PPmethod(PPmethodSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type entro(entroSEXP);
    Rcpp::traits::input_parameter< bool >::type entroindiv(entroindivSEXP);
    rcpp_result_gen = Rcpp::wrap(findproj1D(origclass, origdata, PPmethod, lambda, entro, entroindiv));
    return rcpp_result_gen;
END_RCPP
}
// arma_sub_cond
arma::uvec arma_sub_cond(arma::vec x, int val);
RcppExport SEXP PPtreeExt_arma_sub_cond(SEXP xSEXP, SEXP valSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type val(valSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_sub_cond(x, val));
    return rcpp_result_gen;
END_RCPP
}
// quantileCpp
double quantileCpp(arma::vec x, double probs);
RcppExport SEXP PPtreeExt_quantileCpp(SEXP xSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(quantileCpp(x, probs));
    return rcpp_result_gen;
END_RCPP
}
// quant
NumericVector quant(NumericVector x, NumericVector q);
RcppExport SEXP PPtreeExt_quant(SEXP xSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(quant(x, q));
    return rcpp_result_gen;
END_RCPP
}
// nodestr
arma::vec nodestr(arma::vec classe, arma::vec projdata);
RcppExport SEXP PPtreeExt_nodestr(SEXP classeSEXP, SEXP projdataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type classe(classeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type projdata(projdataSEXP);
    rcpp_result_gen = Rcpp::wrap(nodestr(classe, projdata));
    return rcpp_result_gen;
END_RCPP
}
// findprojwrapMOD
List findprojwrapMOD(arma::vec origclass, arma::mat origdata, std::string PPmethod, double sizep, double lambda, bool entro, bool entroindiv);
RcppExport SEXP PPtreeExt_findprojwrapMOD(SEXP origclassSEXP, SEXP origdataSEXP, SEXP PPmethodSEXP, SEXP sizepSEXP, SEXP lambdaSEXP, SEXP entroSEXP, SEXP entroindivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type origclass(origclassSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type origdata(origdataSEXP);
    Rcpp::traits::input_parameter< std::string >::type PPmethod(PPmethodSEXP);
    Rcpp::traits::input_parameter< double >::type sizep(sizepSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type entro(entroSEXP);
    Rcpp::traits::input_parameter< bool >::type entroindiv(entroindivSEXP);
    rcpp_result_gen = Rcpp::wrap(findprojwrapMOD(origclass, origdata, PPmethod, sizep, lambda, entro, entroindiv));
    return rcpp_result_gen;
END_RCPP
}
// treeconstructMOD
List treeconstructMOD(arma::vec origclass, arma::mat origdata, arma::mat Treestruct, int id, int rep, int rep1, int rep2, arma::mat projbestnode, arma::mat splitCutoffnode, std::string PPmethod, double lambda, double sizep, bool entro, bool entroindiv);
RcppExport SEXP PPtreeExt_treeconstructMOD(SEXP origclassSEXP, SEXP origdataSEXP, SEXP TreestructSEXP, SEXP idSEXP, SEXP repSEXP, SEXP rep1SEXP, SEXP rep2SEXP, SEXP projbestnodeSEXP, SEXP splitCutoffnodeSEXP, SEXP PPmethodSEXP, SEXP lambdaSEXP, SEXP sizepSEXP, SEXP entroSEXP, SEXP entroindivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type origclass(origclassSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type origdata(origdataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Treestruct(TreestructSEXP);
    Rcpp::traits::input_parameter< int >::type id(idSEXP);
    Rcpp::traits::input_parameter< int >::type rep(repSEXP);
    Rcpp::traits::input_parameter< int >::type rep1(rep1SEXP);
    Rcpp::traits::input_parameter< int >::type rep2(rep2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type projbestnode(projbestnodeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type splitCutoffnode(splitCutoffnodeSEXP);
    Rcpp::traits::input_parameter< std::string >::type PPmethod(PPmethodSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type sizep(sizepSEXP);
    Rcpp::traits::input_parameter< bool >::type entro(entroSEXP);
    Rcpp::traits::input_parameter< bool >::type entroindiv(entroindivSEXP);
    rcpp_result_gen = Rcpp::wrap(treeconstructMOD(origclass, origdata, Treestruct, id, rep, rep1, rep2, projbestnode, splitCutoffnode, PPmethod, lambda, sizep, entro, entroindiv));
    return rcpp_result_gen;
END_RCPP
}
// treeconstructIND
List treeconstructIND(arma::vec origclass, arma::mat origdata, arma::mat Treestruct, int id, int rep, int rep1, int rep2, arma::mat projbestnode, arma::mat splitCutoffnode, std::string PPmethod, double lambda, double sizep, bool entro, bool entroindiv, int tot, int iter);
RcppExport SEXP PPtreeExt_treeconstructIND(SEXP origclassSEXP, SEXP origdataSEXP, SEXP TreestructSEXP, SEXP idSEXP, SEXP repSEXP, SEXP rep1SEXP, SEXP rep2SEXP, SEXP projbestnodeSEXP, SEXP splitCutoffnodeSEXP, SEXP PPmethodSEXP, SEXP lambdaSEXP, SEXP sizepSEXP, SEXP entroSEXP, SEXP entroindivSEXP, SEXP totSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type origclass(origclassSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type origdata(origdataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Treestruct(TreestructSEXP);
    Rcpp::traits::input_parameter< int >::type id(idSEXP);
    Rcpp::traits::input_parameter< int >::type rep(repSEXP);
    Rcpp::traits::input_parameter< int >::type rep1(rep1SEXP);
    Rcpp::traits::input_parameter< int >::type rep2(rep2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type projbestnode(projbestnodeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type splitCutoffnode(splitCutoffnodeSEXP);
    Rcpp::traits::input_parameter< std::string >::type PPmethod(PPmethodSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type sizep(sizepSEXP);
    Rcpp::traits::input_parameter< bool >::type entro(entroSEXP);
    Rcpp::traits::input_parameter< bool >::type entroindiv(entroindivSEXP);
    Rcpp::traits::input_parameter< int >::type tot(totSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(treeconstructIND(origclass, origdata, Treestruct, id, rep, rep1, rep2, projbestnode, splitCutoffnode, PPmethod, lambda, sizep, entro, entroindiv, tot, iter));
    return rcpp_result_gen;
END_RCPP
}
// csample_num
arma::vec csample_num(arma::vec x, int size, bool replace, arma::vec prob);
RcppExport SEXP PPtreeExt_csample_num(SEXP xSEXP, SEXP sizeSEXP, SEXP replaceSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type replace(replaceSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(csample_num(x, size, replace, prob));
    return rcpp_result_gen;
END_RCPP
}
// boot
arma::vec boot(arma::mat origclass, arma::mat origdata);
RcppExport SEXP PPtreeExt_boot(SEXP origclassSEXP, SEXP origdataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type origclass(origclassSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type origdata(origdataSEXP);
    rcpp_result_gen = Rcpp::wrap(boot(origclass, origdata));
    return rcpp_result_gen;
END_RCPP
}
// trainfn
arma::vec trainfn(arma::mat origclass, arma::mat origdata, double sizetr);
RcppExport SEXP PPtreeExt_trainfn(SEXP origclassSEXP, SEXP origdataSEXP, SEXP sizetrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type origclass(origclassSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type origdata(origdataSEXP);
    Rcpp::traits::input_parameter< double >::type sizetr(sizetrSEXP);
    rcpp_result_gen = Rcpp::wrap(trainfn(origclass, origdata, sizetr));
    return rcpp_result_gen;
END_RCPP
}
// proximi
arma::mat proximi(arma::mat predtrnt, int m);
RcppExport SEXP PPtreeExt_proximi(SEXP predtrntSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type predtrnt(predtrntSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(proximi(predtrnt, m));
    return rcpp_result_gen;
END_RCPP
}
// mvote
arma::vec mvote(arma::mat votes);
RcppExport SEXP PPtreeExt_mvote(SEXP votesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type votes(votesSEXP);
    rcpp_result_gen = Rcpp::wrap(mvote(votes));
    return rcpp_result_gen;
END_RCPP
}
// oobindex
NumericMatrix oobindex(List datab, int m);
RcppExport SEXP PPtreeExt_oobindex(SEXP databSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type datab(databSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(oobindex(datab, m));
    return rcpp_result_gen;
END_RCPP
}
// oobobs
arma::mat oobobs(arma::mat index);
RcppExport SEXP PPtreeExt_oobobs(SEXP indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type index(indexSEXP);
    rcpp_result_gen = Rcpp::wrap(oobobs(index));
    return rcpp_result_gen;
END_RCPP
}
// mvoteoob
arma::mat mvoteoob(arma::mat votes, arma::mat oobobs);
RcppExport SEXP PPtreeExt_mvoteoob(SEXP votesSEXP, SEXP oobobsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type votes(votesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type oobobs(oobobsSEXP);
    rcpp_result_gen = Rcpp::wrap(mvoteoob(votes, oobobs));
    return rcpp_result_gen;
END_RCPP
}
// ooberrortree
arma::vec ooberrortree(arma::mat votes, arma::mat oobobs, arma::vec classe, int m);
RcppExport SEXP PPtreeExt_ooberrortree(SEXP votesSEXP, SEXP oobobsSEXP, SEXP classeSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type votes(votesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type oobobs(oobobsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type classe(classeSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(ooberrortree(votes, oobobs, classe, m));
    return rcpp_result_gen;
END_RCPP
}
// PPclassification
List PPclassification(arma::mat Treestruct, arma::mat testclassindex, arma::vec IOindex, arma::vec testclass, int id, int rep);
RcppExport SEXP PPtreeExt_PPclassification(SEXP TreestructSEXP, SEXP testclassindexSEXP, SEXP IOindexSEXP, SEXP testclassSEXP, SEXP idSEXP, SEXP repSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Treestruct(TreestructSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type testclassindex(testclassindexSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type IOindex(IOindexSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type testclass(testclassSEXP);
    Rcpp::traits::input_parameter< int >::type id(idSEXP);
    Rcpp::traits::input_parameter< int >::type rep(repSEXP);
    rcpp_result_gen = Rcpp::wrap(PPclassification(Treestruct, testclassindex, IOindex, testclass, id, rep));
    return rcpp_result_gen;
END_RCPP
}
// PPclassindex
List PPclassindex(arma::vec classtemp, arma::mat testclassindex, arma::mat testdata, arma::mat Treestruct, arma::mat AlphaKeep, arma::mat CKeep, int id, int Rule);
RcppExport SEXP PPtreeExt_PPclassindex(SEXP classtempSEXP, SEXP testclassindexSEXP, SEXP testdataSEXP, SEXP TreestructSEXP, SEXP AlphaKeepSEXP, SEXP CKeepSEXP, SEXP idSEXP, SEXP RuleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type classtemp(classtempSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type testclassindex(testclassindexSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type testdata(testdataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Treestruct(TreestructSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type AlphaKeep(AlphaKeepSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type CKeep(CKeepSEXP);
    Rcpp::traits::input_parameter< int >::type id(idSEXP);
    Rcpp::traits::input_parameter< int >::type Rule(RuleSEXP);
    rcpp_result_gen = Rcpp::wrap(PPclassindex(classtemp, testclassindex, testdata, Treestruct, AlphaKeep, CKeep, id, Rule));
    return rcpp_result_gen;
END_RCPP
}
// PPpred
NumericVector PPpred(NumericMatrix TRstr, NumericMatrix TRprnode, NumericMatrix TRspl, NumericMatrix testdata);
RcppExport SEXP PPtreeExt_PPpred(SEXP TRstrSEXP, SEXP TRprnodeSEXP, SEXP TRsplSEXP, SEXP testdataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type TRstr(TRstrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type TRprnode(TRprnodeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type TRspl(TRsplSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type testdata(testdataSEXP);
    rcpp_result_gen = Rcpp::wrap(PPpred(TRstr, TRprnode, TRspl, testdata));
    return rcpp_result_gen;
END_RCPP
}
// imposoon
NumericMatrix imposoon(NumericMatrix train, NumericVector classes, List oobid, List permute, List trees, IntegerVector noob, List TRstrL, List TRsplL, List TRprnodeL);
RcppExport SEXP PPtreeExt_imposoon(SEXP trainSEXP, SEXP classesSEXP, SEXP oobidSEXP, SEXP permuteSEXP, SEXP treesSEXP, SEXP noobSEXP, SEXP TRstrLSEXP, SEXP TRsplLSEXP, SEXP TRprnodeLSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type train(trainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type classes(classesSEXP);
    Rcpp::traits::input_parameter< List >::type oobid(oobidSEXP);
    Rcpp::traits::input_parameter< List >::type permute(permuteSEXP);
    Rcpp::traits::input_parameter< List >::type trees(treesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type noob(noobSEXP);
    Rcpp::traits::input_parameter< List >::type TRstrL(TRstrLSEXP);
    Rcpp::traits::input_parameter< List >::type TRsplL(TRsplLSEXP);
    Rcpp::traits::input_parameter< List >::type TRprnodeL(TRprnodeLSEXP);
    rcpp_result_gen = Rcpp::wrap(imposoon(train, classes, oobid, permute, trees, noob, TRstrL, TRsplL, TRprnodeL));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"PPtreeExt_tableC", (DL_FUNC) &PPtreeExt_tableC, 1},
    {"PPtreeExt_roundme", (DL_FUNC) &PPtreeExt_roundme, 1},
    {"PPtreeExt_LDAindex2", (DL_FUNC) &PPtreeExt_LDAindex2, 4},
    {"PPtreeExt_signC", (DL_FUNC) &PPtreeExt_signC, 1},
    {"PPtreeExt_LDAopt", (DL_FUNC) &PPtreeExt_LDAopt, 5},
    {"PPtreeExt_PDAindex2", (DL_FUNC) &PPtreeExt_PDAindex2, 5},
    {"PPtreeExt_PDAopt", (DL_FUNC) &PPtreeExt_PDAopt, 6},
    {"PPtreeExt_varselect", (DL_FUNC) &PPtreeExt_varselect, 2},
    {"PPtreeExt_datanode", (DL_FUNC) &PPtreeExt_datanode, 2},
    {"PPtreeExt_entropy", (DL_FUNC) &PPtreeExt_entropy, 1},
    {"PPtreeExt_split_entro", (DL_FUNC) &PPtreeExt_split_entro, 2},
    {"PPtreeExt_split_relMOD", (DL_FUNC) &PPtreeExt_split_relMOD, 4},
    {"PPtreeExt_findprojMOD", (DL_FUNC) &PPtreeExt_findprojMOD, 6},
    {"PPtreeExt_findproj1D", (DL_FUNC) &PPtreeExt_findproj1D, 6},
    {"PPtreeExt_arma_sub_cond", (DL_FUNC) &PPtreeExt_arma_sub_cond, 2},
    {"PPtreeExt_quantileCpp", (DL_FUNC) &PPtreeExt_quantileCpp, 2},
    {"PPtreeExt_quant", (DL_FUNC) &PPtreeExt_quant, 2},
    {"PPtreeExt_nodestr", (DL_FUNC) &PPtreeExt_nodestr, 2},
    {"PPtreeExt_findprojwrapMOD", (DL_FUNC) &PPtreeExt_findprojwrapMOD, 7},
    {"PPtreeExt_treeconstructMOD", (DL_FUNC) &PPtreeExt_treeconstructMOD, 14},
    {"PPtreeExt_treeconstructIND", (DL_FUNC) &PPtreeExt_treeconstructIND, 16},
    {"PPtreeExt_csample_num", (DL_FUNC) &PPtreeExt_csample_num, 4},
    {"PPtreeExt_boot", (DL_FUNC) &PPtreeExt_boot, 2},
    {"PPtreeExt_trainfn", (DL_FUNC) &PPtreeExt_trainfn, 3},
    {"PPtreeExt_proximi", (DL_FUNC) &PPtreeExt_proximi, 2},
    {"PPtreeExt_mvote", (DL_FUNC) &PPtreeExt_mvote, 1},
    {"PPtreeExt_oobindex", (DL_FUNC) &PPtreeExt_oobindex, 2},
    {"PPtreeExt_oobobs", (DL_FUNC) &PPtreeExt_oobobs, 1},
    {"PPtreeExt_mvoteoob", (DL_FUNC) &PPtreeExt_mvoteoob, 2},
    {"PPtreeExt_ooberrortree", (DL_FUNC) &PPtreeExt_ooberrortree, 4},
    {"PPtreeExt_PPclassification", (DL_FUNC) &PPtreeExt_PPclassification, 6},
    {"PPtreeExt_PPclassindex", (DL_FUNC) &PPtreeExt_PPclassindex, 8},
    {"PPtreeExt_PPpred", (DL_FUNC) &PPtreeExt_PPpred, 4},
    {"PPtreeExt_imposoon", (DL_FUNC) &PPtreeExt_imposoon, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_PPtreeExt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
