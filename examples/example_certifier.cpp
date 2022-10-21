#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "generateCorrespondences.h"

#include "Essential.h"
#include "EssentialTypes.h"
#include "EssentialUtils.h"
#include "EssentialConstraints.h"

// For the certifier
#include "IterCertAlg/SymmCert.h"
#include "IterCertAlg/SymmCert.cpp"
#include "IterCertAlg/RankYCert.h"
#include "IterCertAlg/RankYCert.cpp"


#include <Eigen/Eigenvalues> 



using namespace std;
using namespace Eigen;
using namespace Essential;

typedef Eigen::Matrix<double, 7, 7> Matrix7;
typedef Eigen::Matrix<double, 8, 8> Matrix8;
typedef Eigen::Matrix<double, 10, 10> Matrix10;
typedef Eigen::Matrix<double, 11, 11> Matrix11;
typedef Eigen::Matrix<double, 15, 15> Matrix15;

typedef Eigen::Matrix<double, 13, 1> Vector13;
typedef Eigen::Matrix<double, 15, 1> Vector15;
typedef Eigen::Matrix<double, 28, 1> Vector28;


int main(int argc, char** argv)
{
    std::cout << "Essential Matrix Estimation with certification!\n";


  double total_time = 0;
  //set experiment parameters
  double noise = 0.05;
  size_t n_points = 100;
  double FoV = 100;  // in degrees
  double max_parallax = 2.;  // in meters
  double min_depth = 1.;     // in meters
  double max_depth = 8.;       // in meters
  double outlier_fraction = 0;


    std::srand(std::time(nullptr));

  


       Vector3 translation;
       Matrix3 rotation;
       bearing_vectors_t points_correspondences;
       Eigen::MatrixXd points_3D(3, n_points);
       std::vector<int> indices_outliers(1, n_points);
       // FoV in degrees
       createSyntheticExperiment(n_points, noise, outlier_fraction, FoV,
                        max_parallax, min_depth, max_depth, translation, rotation,
                        points_correspondences, points_3D, indices_outliers);

      // n_point
      std::cout << "Computing essential matrix GT\n";
      Matrix3 E = computeEfromRt(rotation, translation);



      

      std::cout << "Ground truth rotation:\n" << rotation << std::endl;
      std::cout << "Ground truth translation Tgt:\n" << translation << std::endl;
      std::cout << " \n-------------------------\n";

      EssentialEstimationOptions options;
      options.use_preconditioning = Preconditioner::Dominant_eigenvalues;
      options.verbose = 1;
      options.estimation_verbose=1;
              
      Matrix3 E_init = initialize8pts(points_correspondences);
        
      EssentialClass my_essential_estimation(points_correspondences, options, E_init);

      // run the estimation
      
      std::cout << "Running algotihm on-manifold estimation\n";
      EssentialEstimationResult my_result = my_essential_estimation.getResults();
         

      my_essential_estimation.printResult(my_result);

      std::cout << "Total time: " << my_result.elapsed_estimation_time << std::endl;
      std::cout << "Error rotation: " << distR(rotation, my_result.R_opt) << std::endl;
      std::cout << "Error translation: " << distT(translation, my_result.t_opt) << std::endl;
      
      // create C
      Matrix9 C9 = constructDataMatrix(points_correspondences);
      
      Matrix12 C12 = Matrix12::Zero(); 
      C12.block<9,9>(0,0) = C9; 
      
      Matrix15 C15 = Matrix15::Zero(); 
      C15.block<9,9>(0,0) = C9; 
      
      
      /* Checkin optimality */
      /** Left formulation **/
      std::cout << "########################\nLEFT  FORMULATION\n";
      Vector9 eopt = vec(my_result.E_opt);
      Vector3 t_opt = my_result.t_opt; 
      Matrix3 R_opt = my_result.R_opt; 
      Vector3 q_opt = R_opt.transpose() * t_opt;
      
      Vector12 x_left; 
      x_left << eopt, t_opt;       
      Vector7 mult_left_init = Vector7::Zero();
       
      std::vector<Matrix12> Aleft; 
      createLeftEConstraints(Aleft);     
      
      std::vector<Matrix12> Aright; 
      createRightEConstraints(Aright);                           
   
          
      SymmCert::SymmCertOptions options_symm;

      options_symm.verbose = 1; 
      options_symm.estimation_verbose = 1;
      SymmCert::SymmCertClass<Matrix12, Vector12, Vector7> cert_iter_symm(options_symm);
      
      SymmCert::SymmCertResult<Matrix12, Vector12, Vector7> res_left = cert_iter_symm.getResults(x_left, C12, Aleft, mult_left_init, C12);
      // cert_iter_symm.printResult(res_left); 
      
      
      Eigen::SelfAdjointEigenSolver<Matrix12> eig_left(res_left.opt_Hessian);
      std::cout << "Eigenvalues Hessian LEFT:\n" << eig_left.eigenvalues() << std::endl; 
      
      const int r_left = 8;
      Eigen::Matrix<double, 12, r_left> Yleft; 
      Matrix12 U = eig_left.eigenvectors();  
      Eigen::Matrix<double, 12, r_left> Ured = U.block<12,r_left>(0,4); 
      Matrix8 Dleft = (eig_left.eigenvalues().block<r_left,1>(4,0)).cwiseSqrt().asDiagonal();
      Yleft = Ured * Dleft;
        
      
      
      RankYCert::RankYCertOptions options_ranky;
      
      RankYCert::RankYCertClass<Matrix12, Eigen::Matrix<double, 12, r_left>, Vector12, Vector7> cert_iter_ranky(options_ranky);
     
      
      RankYCert::RankYCertResult<Matrix12, Eigen::Matrix<double, 12, r_left>, Vector7> res_left_ranky = cert_iter_ranky.getResults(x_left, 
                                                                                                           C12, Aleft, mult_left_init, Yleft);
      

      /** Right formulation **/
      std::cout << "########################\nRIGHT FORMULATION\n";
      Vector12 x_right; 
      x_right << eopt, q_opt;       
      Vector7 mult_right_init = Vector7::Zero();

      
      SymmCert::SymmCertResult<Matrix12, Vector12, Vector7> res_right = cert_iter_symm.getResults(x_right, C12, Aright, mult_right_init, C12);
      
            
      // cert_iter_symm.printResult(res_right);  
      
      Eigen::SelfAdjointEigenSolver<Matrix12> eig_right(res_right.opt_Hessian);
      std::cout << "Eigenvalues Hessian RIGHT:\n" << eig_right.eigenvalues() << std::endl; 
      
      const int r_right = 8;
      Eigen::Matrix<double, 12, r_right> Yright; 
      Matrix12 U_right = eig_right.eigenvectors();  
      Eigen::Matrix<double, 12, r_right> Ured_right = U_right.block<12,r_right>(0,4); 
      Matrix8 Dright = (eig_right.eigenvalues().block<r_right,1>(4,0)).cwiseSqrt().asDiagonal();
      Yright = Ured_right * Dright;
      
     
      RankYCert::RankYCertResult<Matrix12, Eigen::Matrix<double, 12, r_right>, Vector7> res_right_ranky = cert_iter_ranky.getResults(x_right, 
                                                                                                          C12, Aright, mult_right_init, Yright);
                                                                                                           
      
      /** Both formulation **/
      std::cout << "########################\nBOTH  FORMULATION\n";
      Vector15 x_both; 
      x_both << eopt, t_opt, q_opt;       
      Vector13 mult_both_init = Vector13::Zero();
       
      std::vector<Matrix15> Aboth; 
      createBothEConstraints(Aboth);   
      

      options_symm.verbose = 1; 
      options_symm.estimation_verbose = 1;
      SymmCert::SymmCertClass<Matrix15, Vector15, Vector13> cert_iter_symm_both(options_symm);

      SymmCert::SymmCertResult<Matrix15, Vector15, Vector13> res_both = cert_iter_symm_both.getResults(x_both, C15, Aboth, mult_both_init, C15);
      
      // cert_iter_symm_both.printResult(res_both); 
      
      Eigen::SelfAdjointEigenSolver<Matrix15> eig_both(res_both.opt_Hessian);
      std::cout << "Eigenvalues Hessian BOTH:\n" << eig_both.eigenvalues() << std::endl; 
      
      const int r_both = 8;
      Eigen::Matrix<double, 15, r_both> Yboth; 
      Matrix15 U_both = eig_both.eigenvectors();  
      Eigen::Matrix<double, 15, r_both> Ured_both = U_both.block<15,r_both>(0,7); 
      Matrix8 Dboth = (eig_both.eigenvalues().block<r_both,1>(7,0)).cwiseSqrt().asDiagonal();
      Yboth = Ured_both * Dboth;
      
      RankYCert::RankYCertClass<Matrix15, Eigen::Matrix<double, 15, r_both>, Vector15, Vector13> cert_iter_ranky_both(options_ranky);
     
      RankYCert::RankYCertResult<Matrix15, Eigen::Matrix<double, 15, r_both>, Vector13> res_both_ranky = cert_iter_ranky_both.getResults(x_both, 
                                                                                                          C15, Aboth, mult_both_init, Yboth);
                                                                                                            
                                                                                                            
      
      /** Adjugate formulation **/
      std::cout << "########################\nADJUGATE  FORMULATION\n";
      Vector15 x_adj; 
      x_adj = x_both;  
      Vector28 mult_adj_init = Vector28::Zero();
       
      std::vector<Matrix15> Aadj; 
      createAdjugateEConstraints(Aadj);   
     
      SymmCert::SymmCertClass<Matrix15, Vector15, Vector28> cert_iter_symm_adj(options_symm);
      SymmCert::SymmCertResult<Matrix15, Vector15, Vector28> res_adj = cert_iter_symm_adj.getResults(x_adj, C15, Aadj, mult_adj_init, C15);
      
      // cert_iter_symm_adj.printResult(res_adj); 
      
      Eigen::SelfAdjointEigenSolver<Matrix15> eig_adj(res_adj.opt_Hessian);
      std::cout << "Eigenvalues Hessian ADJ:\n" << eig_adj.eigenvalues() << std::endl;
      
      const int r_adj = 8;
      Eigen::Matrix<double, 15, r_adj> Yadj; 
      Matrix15 U_adj = eig_adj.eigenvectors();  
      Eigen::Matrix<double, 15, r_adj> Ured_adj = U_adj.block<15,r_adj>(0,7); 
      Matrix8 Dadj = (eig_adj.eigenvalues().block<r_adj,1>(7,0)).cwiseSqrt().asDiagonal();
      Yadj = Ured_adj * Dadj;
      
      RankYCert::RankYCertClass<Matrix15, Eigen::Matrix<double, 15, r_adj>, Vector15, Vector28> cert_iter_ranky_adj(options_ranky);
     
      RankYCert::RankYCertResult<Matrix15, Eigen::Matrix<double, 15, r_adj>, Vector28> res_adj_ranky = cert_iter_ranky_adj.getResults(x_adj, 
                                                                                                        C15, Aadj, mult_adj_init, Yadj);
                                                                                                          


  return 0;

}
