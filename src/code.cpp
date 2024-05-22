#include <Rcpp.h>
//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Eigen/Cholesky> 
using namespace Rcpp;
using namespace Eigen; 


/////////////////////////////////////////////////////////////////////
/// @brief  Compute (X'X+a P)^(-1/2) for Gibbs sampler where P = B'B
/// @param XtX - XtX is the XtX in the design matrix assuming precision of 1
/// @param B   - B is the Choleskey of the precision matrix for the 
///              i.e., the precision matrix is t(B)*B
/// @param a   - a is the muliplier to t(B)*B
/// @return returns (XtX + t(B)*B)^(-1/2)
inline Eigen::MatrixXd XtX_plusP_chol(Eigen::MatrixXd XtX, Eigen::MatrixXd B, double a) {
  XtX.selfadjointView<Eigen::Lower>().rankUpdate(B.adjoint(),a);  
  Eigen::LDLT<Eigen::MatrixXd> myldl(XtX);
  return myldl.solve(Eigen::MatrixXd::Identity(XtX.cols(),XtX.cols())).llt().matrixL();
} 

/// @brief  sample_beta is a function that does one Gibbs sample of each gene's
///         dose response curve given taus, PenMat and phi. The latent correlation is 
///         assumed to be removed in Y. 
/// @param Y - (n x m) Data matrix. Where n represents the number of genes
/// and m is the number of animals. 
/// @param X - (m x r) Regression matrix of spline coefficients where r 
/// is the number of betas that need to be estimated.  It is the same
/// X matrix for each regression
/// @param taus   - (n x 1) matrix of precision matrices for each gene
/// @param PenMat - Penalty matrix for each regression
/// @param phi    - (n x 1) matrix of smoothing coeffiecients 
/// @return betas - beta 
// [[Rcpp::export]]
Eigen::MatrixXd sample_beta (Eigen::MatrixXd Y,
                             Eigen::MatrixXd X,
                             Eigen::MatrixXd taus,
                             Eigen::MatrixXd PenMat,
                             Eigen::MatrixXd phi){
                Y.transposeInPlace();
                Eigen::MatrixXd betas(X.cols(),Y.cols()); 
                
                Eigen::MatrixXd  XtX = Eigen::MatrixXd(X.cols(), X.cols()).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint()); // should only need to do this once
                                                       // accross all  
                X.transposeInPlace();
                 
                /////////////////////////////////////////////////
                // do the number of regressions which is Y.cols()
                #pragma omp parallel
                {
                #pragma omp for
                for (int i = 0; i < Y.cols(); i++){
                    Eigen::MatrixXd tY = Y.col(i); 
                    Eigen::MatrixXd cholTV = pow(1/taus(i,0),0.5)*XtX_plusP_chol(XtX, PenMat, phi(i,0)/taus(i,0));//Eigen::Inverse<Eigen::LDLT> tV(lu);  
                    Eigen::MatrixXd tM = cholTV*cholTV.transpose()*taus(i,0)*X*Y.col(i); 
                    Eigen::MatrixXd ranNorm(X.rows(),1); // =  rnorm(X.cols(), 0, 1); 
                    for (int i = 0; i < X.rows(); i++){
                      ranNorm(i,0) = R::rnorm(0,1); 
                    }
                    betas.col(i) = tM + cholTV*ranNorm; 
                }
                }
                return betas; 

}

///////////////////////////////////////////////
/// @brief given the the latent etas Y compute the correlation component
///        Assumes Y is mean centered, and tau is the gene variance
/// @param Y      - A (g x n)  untrended matrix of observations, where n is the number of 
///                 observations and g is the number of genes                    
/// @param etas   - An (m x n) matrix of latent random effects weighting the correlation for
///                 the given observation, where m is the number of latent factors
/// @param taus   - Precision of the de-meaned Y's 
/// @param priorP - A (m x 1) Prior precision vector for each column of the Lambdas
/// @return       - A (g x m) Gibbs sample estimating Lambda given Y,etas,taus and priorP
// [[Rcpp::export]]
Eigen::MatrixXd compute_lambda(Eigen::MatrixXd Y,
                               Eigen::MatrixXd etas,
                               Eigen::MatrixXd taus,
                               Eigen::MatrixXd cLambda,
                               Eigen::MatrixXd globalP,
                               Eigen::MatrixXd localP){
    Eigen::MatrixXd Lambda(Y.rows(),etas.rows()); // specify (g x m)
    Eigen::MatrixXd tLambda = cLambda.block(0,1,Lambda.rows(),Lambda.cols()-1); // subset Lambda by removing the first column
    Eigen::MatrixXd tEta    = etas.block(1,0,etas.rows()-1,etas.cols()); 

   
    Eigen::MatrixXd R_Y = Y; // R_Y is the residuals of Y - Lambda(-i)Eta(-i)
    for (int i = 0; i < etas.rows(); i++){
          R_Y = Y - tLambda*tEta; 
          Eigen::MatrixXd  sq_eta = taus.array()*Eigen::pow(etas.row(i).array(),2).sum(); 
          sq_eta = sq_eta.array() + globalP(i,0)*localP.col(i).array(); // We now have X'WX + pI
          sq_eta = 1/sq_eta.array().eval(); 
          Eigen::MatrixXd sqrt_eta = sq_eta.array().sqrt(); 
          Eigen::MatrixXd  tY = R_Y.array().colwise()*taus.col(0).array();  //note taus.col(0) turns taus into a vector
          Eigen::MatrixXd  tSum = (tY.array().rowwise()*etas.row(i).array()); //X'W Y
          tSum = tSum.array().rowwise().sum().eval(); 

          Eigen::MatrixXd ranNorm(Lambda.rows(),1); // =  rnorm(X.cols(), 0, 1); 
          for (int i = 0; i < ranNorm.rows(); i++){
              ranNorm(i,0) = R::rnorm(0,1); 
          }
          
          Lambda.col(i) = sqrt_eta.array()*ranNorm.col(0).array() + sq_eta.row(0).array()*tSum.row(0).array(); 
          
          // This if-then works because Lambda and Eta are one column and row (respectively)
          // bigger than tLambda and tEta.  So I am able to modify one row 
          if (i < etas.rows()-1){
            tLambda.col(i) = Lambda.col(i); 
            tEta.row(i) = etas.row(i); 
          }
    }

    return Lambda;
}

////////////////////////////////////////////////
/// @brief  Return the Eta random effects
/// @param Y - A (g x n) matrix of data 
/// @param L - A (g x m) matrix of a factor loading matrix
/// @param T - A (g x 1) matrix of precision values corresponding to 
///            the g genes. 
/// @return  - A (m x n) matrix of a latent loading matrix mean centered with 
///            variance 1.
// [[Rcpp::export]] 
Eigen::MatrixXd compute_eta(Eigen::MatrixXd Y,
                            Eigen::MatrixXd L,
                            Eigen::MatrixXd T){
      
      
      Eigen::MatrixXd    X     =  L.array().colwise()*T.col(0).array().sqrt();
      Eigen::MatrixXd  XtX     =  Eigen::MatrixXd(X.cols(), X.cols()).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint()); // should only need to do this once          
      Eigen::MatrixXd  chol_TV =  XtX_plusP_chol(XtX, Eigen::MatrixXd::Identity(Y.cols(),Y.cols()), 1.0);
      Eigen::MatrixXd       TY = Y.array().colwise()*T.col(0).array();
      Eigen::MatrixXd       TM =  chol_TV*chol_TV.transpose()*L.transpose()*TY;
      Eigen::MatrixXd ran_norm(L.rows(),Y.cols()); 
      for (int i = 0; i < L.rows(); i++){
        for (int j = 0; j < Y.cols(); j++){
            ran_norm(i,j) = R::rnorm(0,1);
        }
      }
      Eigen::MatrixXd Etas  = TM + chol_TV*ran_norm; 
      return Etas;
}


////////////////////////////////////////////////////////
/// @brief  Computes the test statistic between the mean
/// @param X - A (g x n) matrix of data 
/// @param B - A (g x m) matrix of a factor loading matrix
/// @return  - A (m x n) matrix of a latent loading matrix mean centered with 
///            variance 1.
// [[Rcpp::export]] 
Eigen::MatrixXd compute_test_stat(Eigen::MatrixXd X,
                                  Eigen::MatrixXd B
                                  ){
      Eigen::MatrixXd    rVal     = X*B; 
      rVal = rVal.transpose().eval();
      //Rcout << rVal.cols() << ":" << rVal.rows() << std::endl; 
      std::vector<int> ind{0,5};
      std::vector<int> bob{0,3};
      Rcout << X(1,1) << std:: endl; 
      return X; 
}
