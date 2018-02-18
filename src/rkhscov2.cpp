#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* *************************************************************
 * **** rewritting the rkhscov by using svec transformation ****
 * *************************************************************/


/* **************
 * **** misc ****
 * **************/

// [[Rcpp::export]]
vec K_sob_cpp(const vec & s, const vec & t){
  vec k1s = s-0.5;
  vec k1t = t-0.5;
  vec k1abst = abs(s-t)-0.5;

  vec k2s = (k1s % k1s - 1.0/12)/2;
  vec k2t = (k1t % k1t - 1.0/12)/2;

  vec k4abst = (pow(k1abst,4) - k1abst%k1abst/2 + 7.0/240)/24;

 return 1.0 + k1s % k1t + k2s % k2t - k4abst;
}

vec K_sob_cpp2(subview_col<double>  s, double t){
  vec k1s = s-0.5;
  double k1t = t-0.5;
  vec k1abst = abs(s-t)-0.5;

  vec k2s = (k1s % k1s - 1.0/12)/2;
  double k2t = (k1t * k1t - 1.0/12)/2;

  vec k4abst = (pow(k1abst,4) - k1abst%k1abst/2 + 7.0/240)/24;

 return 1.0 + k1s * k1t + k2s * k2t - k4abst;
}

vec K_sob_cpp3(const vec & s, double t){
  vec k1s = s-0.5;
  double k1t = t-0.5;
  vec k1abst = abs(s-t)-0.5;

  vec k2s = (k1s % k1s - 1.0/12)/2;
  double k2t = (k1t * k1t - 1.0/12)/2;

  vec k4abst = (pow(k1abst,4) - k1abst%k1abst/2 + 7.0/240)/24;

 return 1.0 + k1s * k1t + k2s * k2t - k4abst;
}

// [[Rcpp::export]]
mat getK(const vec & t){
  int n = t.n_elem, i;
  mat out(n,n);
  for (i=0; i<n; i++){
    out.submat(i, i, n-1, i) = K_sob_cpp2(t.subvec(i, n-1), t[i]);
  }
  return symmatl(out);
}


// [[Rcpp::export]]
mat getM(const vec & t, int k=150, int p=20, double tol=1e-8){
  int n = t.n_elem, i;
  vec eigval;
  mat X(n,n), Q, R, eigvec;

  for (i=0; i<n; i++){
    X.submat(i, i, n-1, i) = K_sob_cpp2(t.subvec(i, n-1), t[i]);
  }
  X = symmatl(X);

  if (k>n){
    k = n;
  }
  mat y = X * randu(n, k+p); // RNG seems to align with R
  qr_econ(Q, R, y);
  eig_sym(eigval, eigvec, Q.t() * X * Q);
  eigvec = Q * eigvec;
  uvec ind = find((eigval/eigval[eigval.n_elem-1])>tol, 1, "first");
  if (ind.n_elem == 0){
    return zeros(X.n_rows, X.n_cols);
  } else {
    int ii = fmax(ind[0], eigvec.n_cols - k);
    return eigvec.cols(ii, eigvec.n_cols-1) *
      diagmat(sqrt(eigval.subvec(ii, eigval.n_elem-1)));
  }
}


/* ***************************
 * **** utility functions ****
 * ***************************/
// svec(X) = [X_{11}, sqrt(2) X_{21}, sqrt(2) X_{31},...., X_{22}, sqrt(2) X_{32}...]

void svec(const mat & B, vec & Bv){
  int n = B.n_cols, i, j, ind=0;
  for (j=0; j < n; ++j)
    for (i=j; i < n; ++i){
      if (i==j) Bv[ind] = B(i,j);
      else Bv[ind] = B(i,j) * std::sqrt(2.0);
      ind += 1;
    }
}

void svec_inv(mat & B, const vec & Bv){
  int n = B.n_cols, i, j, ind=0;
  for (j=0; j < n; ++j)
    for (i=j; i < n; ++i){
      if (i==j) B(i,j) = Bv[ind];
      else {
        B(i,j) = Bv[ind] / std::sqrt(2.0);
        B(j,i) = B(i,j);
      }
      ind += 1;
    }
}

// [[Rcpp::export]]
vec svec_cpp(const mat & B){
  double r = B.n_cols;
  vec Bv(r * (r+1)/2);
  svec(B, Bv);
  return Bv;
}

// not very efficient
void smat(mat & R, const mat & P, int r){
  int rr = R.n_cols, i, j, ind=0;
  // negative n means unsupplied r
  if (r < 0) r = (int) ((-1.0 + std::sqrt(1.0 + 4.0 * 2.0 * rr)) / 2.0);
  mat Q(P);
  for (j=0; j < r; ++j)
    for (i=j+1; i < r; ++i)
      Q.col(i+r*j) = (P.col(i+r*j) + P.col(j+r*i))/sqrt(2.0);

  for (j=0; j < r; ++j)
    for (i=j+1; i < r; ++i)
      Q.row(i+r*j) = (Q.row(i+r*j) + Q.row(j+r*i))/sqrt(2.0);

  uvec index(rr);
  for (j=0; j < r; ++j)
    for (i=j; i<r; ++i){
      index(ind) = i + j*r;
      ind += 1;
    }
  R = Q.submat(index, index);
}

// **** eigendecomposition of svec form ****
vec eig_svec(const vec & Bv, int r){
  // negative n means unsupplied r
  if (r < 0) r = (int) ((-1.0 + std::sqrt(1.0 + 4.0 * 2.0 * Bv.n_elem)) / 2.0);
  mat B(r, r);
  svec_inv(B, Bv);
  return eig_sym(B);
}

// **** eigendecomposition of svec form ****
void eig_svec(vec & eigval, mat & eigvec, const vec & Bv, int r){
  // negative n means unsupplied r
  if (r < 0) r = (int) ((-1.0 + std::sqrt(1.0 + 4.0 * 2.0 * Bv.n_elem)) / 2.0);
  mat B(r, r);
  svec_inv(B, Bv);
  eig_sym(eigval, eigvec, B);
}

/* ******************************
 * **** compute loss objects ****
 * ******************************/
// Further integration with prep (basic2.R)

void Qlossprep(mat & R, vec & Qv, double * c, List Xs, List Ms, const ivec & include){
  // R, Qv and c are outputs
   
  int *ii = INTEGER(Rf_getAttrib(Ms[0], R_DimSymbol));
  int n = include.n_elem, i, j;
  double lconst;
  mat P = zeros(ii[1]*ii[1], ii[1]*ii[1]);
  mat Q = zeros(ii[1], ii[1]);
  *c = 0.0;

  for (j = 0; j < include.n_elem; ++j)
  {
    i = include[j] - 1; // p.s. R indexing to C indexing
    /* avoid copying */
    int* Mdim = INTEGER(Rf_getAttrib(Ms[i], R_DimSymbol));
    mat M(REAL(Ms[i]), Mdim[0], Mdim[1], false);
    int nX = Rf_length(Xs[i]);
    vec x(REAL(Xs[i]), nX, false);
    
    lconst = 2.0 / (n * nX * (nX - 1.0) );
    mat MM = kron(M, M);
    P += MM.t() * diagmat((vectorise(ones(nX,nX) - eye(nX, nX))) *lconst) * MM;
    mat Z = x * x.t();
    Q -= M.t() * (Z - diagmat(Z)) * M * lconst;
    *c += pow(norm(Z - diagmat(Z), "fro"), 2.0) * lconst / 2.0;
  }

  svec(Q, Qv);
  smat(R, P, Q.n_cols);
}

// [[Rcpp::export]]
List Qlossprep_cpp(List Xs, List Ms, const ivec & include){
  // for R
  int *ii = INTEGER(Rf_getAttrib(Ms[0], R_DimSymbol));
  vec Qv(ii[1] *(ii[1] + 1) / 2);
  mat R(Qv.n_elem, Qv.n_elem);
  double c;
  
  Qlossprep(R, Qv, &c, Xs, Ms, include);
  
  List out;
  out["R"] = wrap(R);
  out["Qv"] = wrap(Qv);
  out["c"] = wrap(c);
  return out;
}

// Qloss class
class Qloss {
  public:
    mat R;
    vec Qv;
    double c;
    int r, rr;
    Qloss(SEXP, SEXP, double, int);
    Qloss(mat &, vec &, double, int);
    double loss(const vec &);
    void grad(vec &, const vec &);
};

// setup loss function objects: P, Qv (no copying)
Qloss::Qloss(SEXP RR, SEXP RQv, double c1, int rr1) : R(REAL(RR), rr1, rr1, false),
  Qv(REAL(RQv), rr1, false) {
  c = c1;
  rr = rr1;
  r = (int) ((-1.0 + std::sqrt(1.0 + 8.0 * rr)) / 2.0);
}


// setup loss function objects: P, Qv (no copying)
Qloss::Qloss(mat & R1, vec & Qv1, double c1, int rr1) : R(R1.memptr(), rr1, rr1, false),
  Qv(Qv1.memptr(), rr1, false) {
  c = c1;
  rr = rr1;
  r = (int) ((-1.0 + std::sqrt(1.0 + 8.0 * rr)) / 2.0);
}

// loss function
double Qloss::loss(const vec & Bv){
  return 0.5 * as_scalar(Bv.t() * R * Bv) + sum(Qv % Bv) + c;
}

// grad function
void Qloss::grad(vec & gBv, const vec & Bv){
  gBv = R * Bv + Qv;
}

// [[Rcpp::export]]
List testQloss(vec Bv, SEXP RR, SEXP RQv, double c){
  int rr = INTEGER(Rf_getAttrib(RR, R_DimSymbol))[0];
  Qloss Q(RR, RQv, c, rr);
  vec gBv(Q.rr);
  Q.grad(gBv, Bv);
  List out;
  out["r"] = wrap(Q.r);
  out["rr"] = wrap(Q.rr);
  out["R"] = wrap(Q.R);
  out["Qv"] = wrap(Q.Qv);
  out["c"] = wrap(Q.c);
  out["lossB"] = wrap(Q.loss(Bv));
  out["grad"] = wrap(gBv);
  return out;
}


double objective(const vec & Bv, double lam1, double lam2, const vec & weight,
    Qloss * Qp, bool pos){
  vec val = eig_svec(Bv, Qp->r);
  if (!pos || all(val>=-1e-8) || (val[0]/val[Qp->r-1]>= -1e-6)){ // middle case to safeguard all zeros
    return Qp->loss(Bv) + lam1 * sum(weight % val) + lam2 * sum(square(val));
  } else {
    return INFINITY;
  }
}

/* **************************
 * **** proximal mapping ****
 * **************************/

// [[Rcpp::export]]
mat prox_trace_hs(const vec & Bv, double s1, double s2, const vec & weight, int r, int
    max_rank, bool pos=true){
  // lam1 sum_j w_j |eig_j| + lam2 sum_j eig_j^2 if pos=false
  vec eigval;
  mat eigvec;
  eig_svec(eigval, eigvec, Bv, r);
  if (pos){
    eigval = (eigval - s1 * weight) / (1.0 + 2.0 * s2);
    uvec ind = find(eigval>0, 1, "first");
    if (ind.n_elem == 0){
      return zeros(r, r);
    } else {
      int temp = fmax(ind[0], eigval.n_elem - max_rank);
      return eigvec.cols(temp, eigvec.n_cols-1) *
        diagmat(eigval.subvec(temp, eigval.n_elem-1)) * (eigvec.cols(temp,
              eigvec.n_cols-1)).t();
    }
  } else {
    // reorder weight;
    uvec oo = sort_index(abs(eigval));
    vec weight1(r);
    for (int k=0; k<r; k++){
      weight1[oo[k]] = weight[k];
    }
    
    vec aeigval = (abs(eigval) - s1 * weight1) / (1.0 + 2.0 * s2);
    eigval = sign(eigval) % aeigval;
    vec aeigval1 = sort(aeigval);
    double val = fmax(aeigval1[fmax(aeigval1.n_elem-max_rank, 0)], 0.0);
    uvec ind = find(aeigval >= val);
    if (ind.n_elem == 0){
      return zeros(r, r);
    } else {
      return eigvec.cols(ind) * diagmat(eigval.elem(ind)) * eigvec.cols(ind).t();
    }
  }
}

void prox_trace_hs_e(vec & out_eigval, mat & out_eigvec, uvec & out_ind,
    const vec & Bv, double s1, double s2, const vec & weight, int r, int max_rank, bool pos=true){
  vec eigval;
  mat eigvec;
  eig_svec(eigval, eigvec, Bv, r);
  if (pos){
    eigval = (eigval - s1 * weight) / (1.0 + 2.0 * s2);
    vec eigval1 = sort(eigval);
    double val = fmax(eigval1[fmax(eigval1.n_elem-max_rank, 0)], 0.0);
    out_ind = find(eigval >= val);
    out_eigval = eigval;
    out_eigvec = eigvec;
  } else {
    // reorder weight;
    uvec oo = sort_index(abs(eigval));
    vec weight1(r);
    for (int k=0; k<r; k++){
      weight1[oo[k]] = weight[k];
    }

    vec aeigval = (abs(eigval) - s1 * weight1) / (1.0 + 2.0 * s2);
    eigval = sign(eigval) % aeigval;
    vec aeigval1 = sort(aeigval);
    double val = fmax(aeigval1[fmax(aeigval1.n_elem-max_rank, 0)], 0.0);
    out_ind = find(aeigval >= val);
    out_eigval = eigval;
    out_eigvec = eigvec;
  }
}

// [[Rcpp::export]]
vec prox(const vec & Bv, double s1, double s2, const vec & weight, int r, int max_rank,
    bool pos=true){
  mat B1;
  vec B1v(Bv.n_elem);
  B1 = prox_trace_hs(Bv, s1, s2, weight, r, max_rank, pos);
  svec(B1, B1v);
  return B1v;
}

void prox_e(vec & out_eigval, mat & out_eigvec, uvec & out_ind,
    const vec & Bv, double s1, double s2, const vec & weight, int r, int max_rank, bool
    pos=true){
  prox_trace_hs_e(out_eigval, out_eigvec, out_ind, Bv, s1, s2, weight, r, max_rank, pos);
}


/* ********************
 * **** Variations ****
 * ********************/

// Auslender-Teboulle06  (see TFOCS Becker-Candes-Grant10)
void steps_AT(vec & B1v, vec & barB1v, const vec & gGv,
    const vec & barB0v, const vec & B0v, const vec & weight, double L1,
    double theta1, double lam1, double lam2, int r, int max_rank, bool pos)
{
  // B1v, barB1v are outputs
  vec tGv = barB0v-gGv / (L1 * theta1);
  barB1v = prox(tGv, lam1/(L1*theta1), lam2/(L1*theta1), weight, r, max_rank, pos);
  B1v = (1.0 - theta1) * B0v + theta1 * barB1v;
}

// Neterov's 1983 method  (see TFOCS Becker-Candes-Grant10)
// or FISTA (with step size adaption from Becker-Candes-Grant10)
void steps_N(vec & B1v, vec & barB1v, const vec & Gv, const vec & gGv,
    const vec & B0v, const vec & weight, double L1, double theta1, double lam1, double lam2,
    int r, int max_rank, bool pos)
{
  // B1v, barB1v are outputs
  vec tGv = Gv-gGv/L1;
  B1v = prox(tGv, lam1/L1, lam2/L1, weight, r, max_rank, pos);
  barB1v = (1.0 / theta1) * B1v - ((1.0 - theta1) / theta1) * B0v;
}

// Classical projected gradient generalization or Proximal gradient descent
// (see TFOCS Becker-Candes-Grant10)
void steps_GRA(vec & B1v, vec & barB1v, const vec & Gv, const vec & gGv,
    const vec & B0v, const vec & weight, double L1, double theta1, double lam1, double lam2,
    int r, int max_rank, bool pos) 
{
  // B1v, barB1v are outputs
  vec tGv = Gv-gGv/L1;
  B1v = prox(tGv, lam1/L1, lam2/L1, weight, r, max_rank, pos);
  barB1v = B1v;
}

/* *********************
 * **** L condition ****
 * *********************/

double lip_const(const vec & Bv, const vec & gBv, const vec & Gv,
    const vec & gGv, Qloss * Qp){
  //double out = (Qp->loss(Gv));
  vec Av = Bv - Gv;
  double res = 2.0 * std::max( (Qp->loss(Bv)) - (Qp->loss(Gv)) - sum(Av % gGv), 0.0) / sum(square(Av));
  if (isnan(res)){
    return 0.0;
  } else {
    return res;
  }
}

double lip2_const(const vec & Bv, const vec & gBv, const vec & Gv, const vec & gGv){
  vec Av = Bv - Gv;
  double res = 2.0 * std::abs(sum(Av % (gBv-gGv))) / sum(square(Av));
  if (isnan(res)){
    return 0.0;
  } else {
    return res;
  }
}

/***************************************
 **** accelarated proximal gradient ****
 ***************************************/

// **** check the rank of B0v ****
// output the objective
double checkB0v(vec & B0v, Qloss * Qp, double lam1, double lam2, const vec &
    weight, int * max_rank, bool pos){
  if (*max_rank<0){
    *max_rank = Qp->r;
  } else if (*max_rank > (Qp->r)) {
    Rcout << "max_rank > number of rows of B0: forcing max_rank = number of rows of B0\n";
    *max_rank = Qp->r;
  } else {
    vec B0eigval;
    mat B0eigvec;
    eig_svec(B0eigval, B0eigvec, B0v, Qp->r);
    int temp = B0eigvec.n_cols - *max_rank;
    mat B0 = B0eigvec.cols(temp, B0eigvec.n_cols-1) *
        diagmat(B0eigval.subvec(temp, B0eigval.n_elem-1)) *
        (B0eigvec.cols(temp, B0eigvec.n_cols-1)).t();
    svec(B0, B0v);
  }
  double out = objective(B0v, lam1, lam2, weight, Qp, pos);
  if (out==INFINITY) Rf_error("B is not positive semi-definite!");
  return out;
}

int rkhscov_pg(vec & eigval, mat & eigvec, uvec & ind, vec & vals, vec & dBs, mat & conv,
    Qloss * Qp, double lam, double gam, vec B0v, const vec & weight, double L=1.0,
    double eta=2.0, double alpha=0.9, int maxit=10000, bool traceit=true, double tol=1e-8,
    int max_rank=-1, bool pos=true, int variant=2, int cond=2) {
  // outputs: eigval, eigvec, ind, vals, dBs, conv
  // return iter
  // proximal gradient methods
  // loss = 0.5 * Bv^T R Bv + Qv^T Bv + c (note Bv and Qv are svec version of B and Q)
  // penalty:
  // lam * (gam * sum_j w_j |eig_j| + (1-gam)/2 * sum_j eig_j^2) if pos=false
  // note: no weight on eig_j^2

  // setup lam1 and lam2
  double lam1 = lam * gam;
  double lam2 = lam * (1.0 - gam) / 2.0;

  if ((vals.n_elem < maxit) || (dBs.n_elem < maxit) || (conv.n_rows < maxit)){
    Rf_error("rkhscov_pg: (vals.n_elem < maxit) || (dBs.n_elem < maxit) || (conv.n_rows < maxit)");
  }

  // setup Qloss object
  int rr = Qp->rr, r = Qp->r;
  
  // checking
  if (L<0) Rf_error("L need to be >0!");
  vals[0] = checkB0v(B0v, Qp, lam1, lam2, weight, &max_rank, pos);

  // initialization
  dBs[0] = 0.0;
  conv(0,0) = 0.0; conv(0,1) = 0.0;
  int count = 0, iter;
  vec B1v = B0v, barB0v = B0v;
  vec barB1v(rr), Gv(rr), gGv(rr), tGv(rr), gB1v(rr), tempeigval(r);
  double theta0 = 1e9, theta1, tval, L1, L0 = L, hatL, lossB;

  for (iter=1; iter < maxit; ++iter){
    if (traceit) Rcout << "iter " << iter-1 << "\t obj val: " << vals[iter-1] <<
      "\t rel. changes (obj val, Bv): " << conv(iter-1,1) << "," << conv(iter-1,1);
    L1 = L0 * alpha;
    while (1){
      theta1 = 2.0 / (1.0 + sqrt( 1.0 + 4.0 * L1 / (L0 * std::pow(theta0, 2.0))));
      Gv = (1.0 - theta1) * B0v + theta1 * barB0v;
      Qp->grad(gGv, Gv); // compute gGv

      if (variant==1){
        // Auslender-Teboulle06  (see TFOCS Becker-Candes-Grant10)
        steps_AT(B1v, barB1v, gGv, barB0v, B0v, weight, L1, theta1, lam1, lam2, r, max_rank,
            pos);
      } else if (variant==2) {
        // Neterov's 1983 method  (see TFOCS Becker-Candes-Grant10)
        // or FISTA (with step size adaption from Becker-Candes-Grant10)
        steps_N(B1v, barB1v, Gv, gGv, B0v, weight, L1, theta1, lam1, lam2, r, max_rank, pos);
      } else if (variant==3) {
        // Classical projected gradient generalization or Proximal gradient descent
        // (see TFOCS Becker-Candes-Grant10)
        steps_GRA(B1v, barB1v, Gv, gGv, B0v, weight, L1, theta1, lam1, lam2, r, max_rank, pos);
      } else {
        Rf_error("No such variant!");
      }

      Qp->grad(gB1v, B1v); // compute gB1v
      lossB = (Qp->loss(B1v));
    
      if (cond==1){
        // typical
        L = lip_const(B1v, gB1v, Gv, gGv, Qp);
      } else if (cond==2) {
        // suggested by TFOCS Becker-Candes-Grant10
        L = lip2_const(B1v, gB1v, Gv, gGv);
      } else {
        Rf_error("No such cond!");
      }
      
      if (L1 >= L) break;
      L1 = std::max(eta * L1, L);
      if (traceit) Rcout << ".";
    }
    if (traceit) Rcout << "\n";

    
    tempeigval = eig_svec(B1v, r);
    tval = lossB + lam1 * sum(weight % tempeigval) + lam2 * sum(square(tempeigval));

    vals[iter] = tval;
    dBs[iter] = sum(square(B1v-B0v));
    conv(iter,0) = std::abs(vals[iter]-vals[iter-1])/vals[iter-1];
    conv(iter,1) = dBs[iter] / sum(square(B0v));
    
    B0v = B1v;
    barB0v = barB1v;
    L0 = L1;
    theta0 = theta1;

    if (((conv(iter,0) < tol)&&(conv(iter,1) < tol))||(norm(B1v)==0)){
      count += 1;
    } else
      count = 0;
    if (count >= 10){
      iter += 1;
      break;
    }

  }

  // final proximal step to get thresholded solution

  Gv = B1v;
  Qp->grad(gGv, Gv); // compute gB1v
  tGv = Gv-gGv/L1;
  prox_e(eigval, eigvec, ind, tGv, lam1/L1, lam2/L1, weight, r, max_rank, pos);

  return iter;
}



// [[Rcpp::export]]
List rkhscov_pg_cpp(SEXP RR, SEXP RQv, double c, double lam, double gam, vec
    B0v, vec weight, double L=1.0, double eta=2.0, double alpha=0.9, int
    maxit=10000, bool traceit=true, double tol=1e-8, int max_rank=-1, bool
    pos=true, int variant=2, int cond=2) {
  // R wrapper of rkhscov_pg
  vec vals(maxit), dBs(maxit), eigval;
  mat conv(maxit, 2), eigvec;
  uvec ind;
  int iter;
  
  int rr = INTEGER(Rf_getAttrib(RR, R_DimSymbol))[0];
  Qloss Q(RR, RQv, c, rr);

  iter = rkhscov_pg(eigval, eigvec, ind, vals, dBs, conv, &Q, lam, gam, B0v, weight, L,
     eta, alpha, maxit, traceit, tol, max_rank, pos, variant, cond);

  List outlist, e;
  /* currently, wrap is not defined for subview_col: need to copy to eigval1, eigvec1 */
  vec eigval1;
  mat eigvec1;
  vec vals1;
  vec dBs1;
  mat conv1;
  vals1 = vals.subvec(0, iter-1);
  dBs1 = dBs.subvec(0, iter-1);
  conv1 = conv.rows(0, iter-1);
  if (ind.n_elem == 0){
    eigval1 = zeros(1);
    eigvec1 = zeros(eigvec.n_rows, 1);
  } else {
    eigval1 = eigval.elem(ind);
    eigvec1 = eigvec.cols(ind);
  }
  e["values"] = wrap(eigval1);
  Rf_setAttrib(e["values"], R_DimSymbol, Rf_ScalarReal(eigval1.n_elem)); // flatten
  e["vectors"] = wrap(eigvec1);

  outlist["obj"] = wrap(vals1);
  Rf_setAttrib(outlist["obj"], R_DimSymbol, Rf_ScalarReal(vals1.n_elem)); // flatten
  outlist["e"] = wrap(e);
  outlist["dBs"] = wrap(dBs1);
  outlist["conv"] = wrap(conv1);
  outlist["weight"] = wrap(weight);
  return outlist;
}



/* *********************************
 * **** k-fold cross-validation ****
 * *********************************/

// [[Rcpp::export]]

double start_lambda(SEXP RR, SEXP RQv, double c, vec weight, double gam, bool
    pos=true, double rtol=1e-4) {
  // loss = 0.5 * Bv^T R Bv + Qv^T Bv + c (note Bv and Qv are svec version of B and Q)
  // determine the largest lambda in the grid
  // find a lambda such that prox( B-gB, lambda/L) = 0 when B=0
  // for gam=0 (<1e-20), we cannot shrink all to zero with a finite lambda, hence
  // we try to find a lambda such that the largest (abs) eigenvalue is rtol times itself (if !pos)

  // setup Qloss object
  int rr = INTEGER(Rf_getAttrib(RR, R_DimSymbol))[0];
  Qloss Q(RR, RQv, c, rr);
  int r = Q.r;
  double L = (eig_sym(Q.R))[rr-1]; // lipschitz constant of the gradient

  vec Gv(rr, fill::zeros);
  vec gGv(rr);
  Q.grad(gGv, Gv);
  vec tGv = -gGv/L;
  vec eigval = eig_svec(tGv, r);
  if (gam>1e-20){ // L is not needed
    if (pos)
      return std::max(0.0, eigval[r-1]) * L / (gam * weight[r-1]);
    else
      return std::max(std::abs(eigval[0]), std::abs(eigval[r-1])) * L / (gam * weight[r-1]);
  } else {
    // when gam==0: only l2 penalty
    return (1.0/rtol -1.0) * 0.5 * L; // relative tol
  }
}

// [[Rcpp::export]]
double sumvaliderr(SEXP RR, SEXP RQv, double c, const vec & Bv){
  int rr = INTEGER(Rf_getAttrib(RR, R_DimSymbol))[0];
  Qloss Q(RR, RQv, c, rr);
  return Q.loss(Bv);
}


// [[Rcpp::export]]
mat rkhscov_pg_cvj_cpp(List Xs, List Ms, ivec traingroup, ivec testgroup, const
    vec & lams, double gam,
    vec B0v, vec weight, double L=1.0, double eta=2.0, double alpha=0.9, int
    maxit=1000, bool traceit=true, double tol=1e-8, int max_rank=-1, bool
    pos=true, int variant=2, int cond=2) {
  // loss = 0.5 * Bv^T R Bv + Qv^T Bv + c (note Bv and Qv are svec version of B and Q)

  int n = Xs.size(), nj = testgroup.n_elem, nlam = lams.n_elem;
  int *ii = INTEGER(Rf_getAttrib(Ms[0], R_DimSymbol));
  int r = ii[1], rr = ii[1] *(ii[1] + 1) / 2;
  if (B0v.n_elem != rr) Rf_error("length of B0v does not match!");

  // prepare train and test Qloss object
  mat Rtrain(rr, rr), Rtest(rr, rr);
  vec Qvtrain(rr), Qvtest(rr);
  double ctrain, ctest, err;

  Qlossprep(Rtrain, Qvtrain, &ctrain, Xs, Ms, traingroup);
  Qlossprep(Rtest, Qvtest, &ctest, Xs, Ms, testgroup);
  Qloss Qtrain(Rtrain, Qvtrain, ctrain, rr), Qtest(Rtest, Qvtest, ctest, rr);

  // containers
  vec vals(maxit), dBs(maxit), eigval(r);
  mat conv(maxit, 2), eigvec(r, r), B0(r, r), cvj(nlam, 2);
  uvec ind;

  for (int i = nlam-1; i >= 0; --i){
    rkhscov_pg(eigval, eigvec, ind, vals, dBs, conv, &Qtrain, lams[i], gam, B0v, weight, L, eta,
        alpha, maxit, traceit, tol, max_rank, pos, variant, cond);
    if (ind.n_elem == 0)
      B0 = zeros(r, r);
    else
      B0 = eigvec.cols(ind) * diagmat(eigval.elem(ind)) * eigvec.cols(ind).t();
    svec(B0, B0v);
    err = Qtest.loss(B0v);
    cvj(i,0) = err; cvj(i,1) = err * nj;
/*    if (pos && (min(eigval)<= -1e-8)){
      B0 = eye(r,r) * 1e-8;
      svec(B0, B0v);
    }*/
    if (max(abs(eigval)) <= 1e-15){
      B0 = eye(r,r) * 1e-8;
      svec(B0, B0v);
    }
  }
  return cvj;
}


