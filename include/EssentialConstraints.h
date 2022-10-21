#include "EssentialTypes.h"


namespace Essential
{

        typedef Eigen::Matrix<double, 15, 15> Matrix15;
        
        void createLeftEConstraints (std::vector<Matrix12> & As);
        
        void createRightEConstraints(std::vector<Matrix12> & As);

        void createBothEConstraints(std::vector<Matrix15> & As);
                
        void createAdjugateEConstraints(std::vector<Matrix15> & As);
        
}  // end of namespace Essential
