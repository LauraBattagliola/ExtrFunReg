#include <RcppArmadillo.h>
#include <Rinternals.h>




using namespace Rcpp;



extern "C" {
  void metrl2_(double *A, double *B, int *M, int *N, int *D, double *METR);
  void metrl2b_(double *A, double *B, int *M, int *N, int *D, double *METR);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix metrl2_cpp(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B) {
  int M = A.nrow();
  int N = B.nrow();
  int D = A.ncol();
  Rcpp::NumericMatrix METR(M, N);

  metrl2_(&A[0], &B[0], &M, &N, &D, &METR[0]);

  return METR;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix metrl2b_cpp(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B) {
  int M = A.nrow();
  int N = B.nrow();
  int D = A.ncol();
  Rcpp::NumericMatrix METR(M, N);

  metrl2b_(&A[0], &B[0], &M, &N, &D, &METR[0]);

  return METR;
}



// [[Rcpp::export]]
double kern(double t, int k) {
  // positive kernel functions
  // k is the switch between kernels
  double result;
  double pi = atan(1)*4;

  result = 0;
  if(k==1){ // uniform kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = 1;
  }
  if(k==2){ // triangular kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = (1-t)*2;
  }
  if(k==3){ // Epanechnikov kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = (1-pow(t,2))*3/2;
  }
  if(k==4){ // biweight kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = pow(1-pow(t,2),2)*15/8;
  }
  if(k==5){ // triweight kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = pow(1-pow(t,2),3)*35/16;
  }
  if(k==6){ // Gaussian kernel
    if (t < 0)
      result = 0;
    else
      result = sqrt(2/pi)*exp(-pow(t,2)/2);
  }
  return result;
}

// [[Rcpp::export]]
double J_tau(double tau, double t)
{
  double res;
  res=1.0;
  if (tau == 1.0){
    return 1.0;
  }
  if (tau > 0.0 && tau <= 0.5) {
    double s_tau;
    s_tau = log(0.5)/log(1.0-tau);
    res=s_tau*pow((1-t),(s_tau-1.0));
  }
  if (tau > 0.5 && tau < 1.0) {
    double r_tau;
    r_tau = log(0.5)/log(tau);
    res=r_tau*pow(t,(r_tau-1.0));
  }
  return res;
}


// [[Rcpp::export]]
double kern_CDF(double t, int k) {
  // positive kernel functions
  // k is the switch between kernels
  double result;
  double pi = atan(1)*4;

  result = 0;

  if(k==2){ // triangular kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = (1-t)*2;
  }
  if(k==3){ // Epanechnikov kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = (1-pow(t,2))*3/2;
  }
  if(k==4){ // biweight kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = pow(1-pow(t,2),2)*15/8;
  }
  if(k==5){ // triweight kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = pow(1-pow(t,2),3)*35/16;
  }

  return result;
}

// [[Rcpp::export]]
double cdf_ksm_single(Rcpp::NumericVector D, Rcpp::NumericVector Y, double y, Rcpp::NumericVector h, int n, int kernCDFI) {
  double res = 0;
  double den = 0;
  double indicator, weight;

  for (int i = 0; i < n; i++) {
    indicator = (Y[i] <= y) ? 1.0 : 0.0;  // Indicator function
    weight = kern_CDF(D[i] / h[i], kernCDFI);
    res += indicator * weight;
    den += weight;
  }

  return den == 0 ? 0 : res / den;  // Return the computed CDF value for point y
}




// [[Rcpp::export]]
Rcpp::NumericMatrix multiplyMatrices(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B) {
  // Convert Rcpp::NumericMatrix to arma::mat
  arma::mat armaA = Rcpp::as<arma::mat>(A);
  arma::mat armaB = Rcpp::as<arma::mat>(B);

  // Perform matrix multiplication using Armadillo
  arma::mat result = armaA * armaB;

  // Convert the result back to Rcpp::NumericMatrix
  return Rcpp::wrap(result);
}


// [[Rcpp::export]]
Rcpp::NumericVector multiplyMatrixVector(Rcpp::NumericMatrix M, Rcpp::NumericVector v) {
  // Convert Rcpp::NumericMatrix and Rcpp::NumericVector to arma::mat and arma::vec
  arma::mat armaM = Rcpp::as<arma::mat>(M);
  arma::vec armaV = Rcpp::as<arma::vec>(v);

  // Perform matrix-vector multiplication using Armadillo
  arma::vec result = armaM * armaV;

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(result);
}



// [[Rcpp::export]]
Rcpp::List solveLinearSystem(Rcpp::NumericMatrix A, Rcpp::NumericVector b) {
  // Convert Rcpp::NumericMatrix and Rcpp::NumericVector to arma::mat and arma::vec
  arma::mat armaA = Rcpp::as<arma::mat>(A);
  arma::vec armaB = Rcpp::as<arma::vec>(b);

  arma::vec x;
  bool success = true;

  // Regularization parameter (small value)
  double lambda = 0.000000000001;

  // Check for non-singularity and ill-conditioning using condition number
  double cond = arma::cond(armaA);
  if (cond > 1e10 || std::isinf(cond) || std::isnan(cond)) {  // Adjust threshold as appropriate
    // Regularize the matrix
    armaA.diag() += lambda;
  }


  // Solve the system using regular solve (which uses LAPACK)
  success = arma::solve(x, armaA, armaB);


  // Prepare the return value
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("solution") = Rcpp::wrap(x),
    Rcpp::Named("success") = success
  );

  return result;
}




/////////////////////////////////////////////////////////////////////////////////
// fllr procedures used for CV to find local bandwith h_i
// [[Rcpp::export]]
Rcpp::NumericMatrix llsm_single_cpp(Rcpp::NumericMatrix C, Rcpp::NumericMatrix Cnew, Rcpp::NumericMatrix D, Rcpp::NumericVector Y, Rcpp::NumericVector h, int n, int m, int J, int kernI) {
  Rcpp::NumericMatrix C1(n, J), C2(J, n), C3(n, J), M1(J, J),  res(m, J);
  Rcpp::NumericVector restemp(J);
  //int info;

  for(int k = 0; k < m; k++) { // for each function
    // create matrix with rows coefficients corresponding to X - x
    for(int i = 0; i < n; i++) {
      C1(i, 0) = 1; // C1[i,1] = 1
      for(int j = 1; j < J; j++) {
        C1(i, j) = C(i, j) - Cnew(k, j); // C1[i,j] = C[i,j] - Cnew[k,j], i>1
      }
    }

    restemp.fill(0.0);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < J; j++) {
        // C2 is C1 transposed, and cut at jj
        // C3 is C1, cut at jj
        // double kern_val = kern(abs(D(i,k) / h[k]), kernI)
        double kern_val = kern((D(k,i) / h[k]), kernI);;
        C2(j, i) = C1(i, j) * (h[k] == -1 ? 1 : kern_val);
        C3(i, j) = C1(i, j);
      }
    }


    M1 = multiplyMatrices(C2,C3);
    restemp = multiplyMatrixVector(C2,Y);
    Rcpp::List inv = solveLinearSystem(M1, restemp);

    // Check the success flag and process the solution
    if(!Rcpp::as<bool>(inv["success"])) {
      for(int j = 0; j < J; j++) {
        res(k, j) = R_PosInf;
      }
    } else {
      Rcpp::NumericVector solution = Rcpp::as<Rcpp::NumericVector>(inv["solution"]);
      for(int j = 0; j < J; j++) {
        res(k, j) = solution[j];
      }
    }


  }


  return res;
}



// [[Rcpp::export]]
void llsm_cv_single_cpp(Rcpp::NumericMatrix C, Rcpp::NumericMatrix D, Rcpp::NumericVector Y, int n, int J, Rcpp::NumericMatrix H, int nH, Rcpp::NumericVector CV, Rcpp::NumericVector CVB, int nCV, int kernI) {
  int i, j, k, l, hi, p;
  Rcpp::NumericVector Yi(n - 1), h_vec(n - 1);
  Rcpp::NumericMatrix Ci(n - 1, J), Cinew(1, J), Di(1,n - 1),  Yt(1,J);
  double h;

  // Initialize CV and CVB
  std::fill(CV.begin(), CV.end(), 0.0);
  std::fill(CVB.begin(), CVB.end(), 0.0);

  // for each bandwidth h leave-one-out cross-validation
  for (hi = 0; hi < nH; hi++) {
    for (i = 0; i < nCV; i++) {
      // cross-validation for first nCV functions only in a leave-one-out sense
      for (j = 0; j < J; j++) {
        Cinew(0, j) = C(i, j);
      }

      k = 0;
      for (j = 0; j < n; j++) {
        if (j != i) {
          Di(0,k) = D(j, i);
          Yi[k] = Y[j];
          for (l = 0; l < J; l++) {
            Ci(k, l) = C(j, l);
          }
          k++;
        }
      }

    // h = H(hi, i);
    h = H(i, hi);

      for (p = 0; p < (n-1); p++) {
        h_vec[p]=h;
      }


      Yt = llsm_single_cpp(Ci, Cinew, Di, Yi,  h_vec, n-1, 1, J, kernI); // Adjusted call to llsm_single


      //
      // cross-validate
      CV[(J - 1) * nH + hi] += pow(Yt[0] - Y[i], 2);  // MSE
      CVB[(J - 1) * nH + hi] += Yt[0] - Y[i];        // bias
    }


     CV[(J - 1) * nH + hi] /= nCV;
    CVB[(J - 1) * nH + hi] = CV[(J - 1) * nH + hi] - pow(CVB[(J - 1) * nH + hi] / nCV, 2); // variance (MSE - bias^2)

  }


}



// [[Rcpp::export]]
Rcpp::NumericMatrix efr_single(Rcpp::NumericMatrix C, Rcpp::NumericMatrix Cnew, Rcpp::NumericMatrix D, Rcpp::NumericVector Y, Rcpp::NumericVector h, Rcpp::NumericVector hF, int n, int m, int J, int kernI, double tau, int kernCDFI) {

  Rcpp::NumericMatrix C1(n, J), C2(J, n), C3(n, J), M1(J, J), Fhat(m,n) ,res(m, J);
  Rcpp::NumericVector restemp(J);

  for (int i = 0; i < m; i++) {
    Rcpp::NumericVector Dd = D.row(i);
    for (int j = 0; j < n; j++) {
      Fhat(i, j) = cdf_ksm_single(Dd, Y, Y[j], hF, n, kernCDFI);
    }
  }

  for(int k = 0; k < m; k++) { // for each function
    // create matrix with rows coefficients corresponding to X - x
    for(int i = 0; i < n; i++) {
      C1(i, 0) = 1; // C1[i,1] = 1
      for(int j = 1; j < J; j++) {
        C1(i, j) = C(i, j) - Cnew(k, j); // C1[i,j] = C[i,j] - Cnew[k,j], i>1
      }
    }

    restemp.fill(0.0);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < J; j++) {
        double kern_val = kern((D(k,i) / h[k]), kernI);
        double tau_weight = J_tau(tau,Fhat(k,i));

        C2(j, i) = C1(i, j) * (h[k] == -1 ? 1 : kern_val) * tau_weight ;
        C3(i, j) = C1(i, j);

      }
    }




    M1 = multiplyMatrices(C2,C3);
    restemp = multiplyMatrixVector(C2,Y);
    Rcpp::List inv = solveLinearSystem(M1, restemp);

    // Check the success flag and process the solution
    if(!Rcpp::as<bool>(inv["success"])) {
      for(int j = 0; j < J; j++) {
        res(k, j) = R_PosInf;
      }
    } else {
      Rcpp::NumericVector solution = Rcpp::as<Rcpp::NumericVector>(inv["solution"]);
      for(int j = 0; j < J; j++) {
        res(k, j) = solution[j];
      }
    }
  }




  return res;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix efr_single_leave(Rcpp::NumericMatrix C, Rcpp::NumericMatrix Cnew, Rcpp::NumericMatrix D, Rcpp::NumericVector Y, Rcpp::NumericVector h, Rcpp::NumericVector hF, int n, int m, int J, int kernI, double tau, int kernCDFI) {

  Rcpp::NumericMatrix C1(n, J), C2(J, n), C3(n, J), M1(J, J),res(m, J), Fhat(m,n);
  Rcpp::NumericVector restemp(J), tauweights(m);//, Fhat(m);
  double cnst;



  for (int i = 0; i < m; i++) {
    Rcpp::NumericVector Dd = D.row(i);
    for (int j = 0; j < n; j++) {
      Fhat(i, j) = cdf_ksm_single(Dd, Y, Y[j], hF, n, kernCDFI);
    }
  }


  for(int k = 0; k < m; k++) { // for each function
    // create matrix with rows coefficients corresponding to X - x
    for(int i = 0; i < n; i++) {
      C1(i, 0) = 1; // C1[i,1] = 1
      for(int j = 1; j < J; j++) {
        C1(i, j) = C(i, j) - Cnew(k, j); // C1[i,j] = C[i,j] - Cnew[k,j], i>1
      }
    }

    restemp.fill(0.0);
    for(int i = 0; i < n; i++) {
      if(i==k) cnst=0.0; else cnst=1.0;
      for(int j = 0; j < J; j++) {

        double kern_val = kern((D(k,i) / h[k]), kernI);
        double tau_weight = J_tau(tau,Fhat(k,i));
        C2(j, i) = C1(i, j) * (h[k] == -1 ? 1 : kern_val) * tau_weight*cnst ;
        C3(i, j) = C1(i, j);

      }
    }



    M1 = multiplyMatrices(C2,C3);
    restemp = multiplyMatrixVector(C2,Y);
    Rcpp::List inv = solveLinearSystem(M1, restemp);

    // Check the success flag and process the solution
    if(!Rcpp::as<bool>(inv["success"])) {
      for(int j = 0; j < J; j++) {
        res(k, j) = R_PosInf;
      }
    } else {
      Rcpp::NumericVector solution = Rcpp::as<Rcpp::NumericVector>(inv["solution"]);
      for(int j = 0; j < J; j++) {
        res(k, j) = solution[j];
      }
    }
  }




  return res;
}

