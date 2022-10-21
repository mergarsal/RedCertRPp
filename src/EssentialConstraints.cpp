#include "EssentialConstraints.h"




namespace Essential
{        
    void createLeftEConstraints(std::vector<Matrix12> & As)
    {
        
        size_t e1=0, e4=e1+1, e7=e4+1, e2=e7+1,e5=e2+1,e8=e5+1,e3=e8+1,e6=e3+1,e9=e6+1; 
        size_t t1=e9+1, t2=t1+1, t3=t2+1; 
        
        As.clear();
        /** Left formulation **/
        Matrix12 Ai; 
        Ai.setZero(); 
        Ai(t1,t1)=1;
        Ai(t2,t2)=1;
        Ai(t3,t3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();         
        Ai(e1,e1)=1;   Ai(e2,e2)=1;Ai(e3,e3)=1;
        Ai(t2, t2)=-1; Ai(t3, t3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  
        
        Ai.setZero();  Ai(e1,e4)=1; Ai(e2,e5)=1;Ai(e3,e6)=1;
        Ai(t1, t2)=1;           
        As.push_back(0.5 * (Ai + Ai.transpose()));      
                   
        Ai.setZero();  Ai(e1,e7)=1; Ai(e2,e8)=1;Ai(e3,e9)=1;
        Ai(t1, t3)=1;       
        As.push_back(0.5 * (Ai + Ai.transpose()));          
                   
        Ai.setZero();  Ai(e4,e4)=1; Ai(e5,e5)=1;Ai(e6,e6)=1;
        Ai(t1, t1)=-1;  Ai(t3, t3)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero();  Ai(e4,e7)=1; Ai(e5,e8)=1;Ai(e6,e9)=1;
        Ai(t2, t3)=1;     
        As.push_back(0.5 * (Ai + Ai.transpose()));            
                   
        Ai.setZero();  Ai(e7,e7)=1; Ai(e8,e8)=1;Ai(e9,e9)=1;
        Ai(t1, t1)=-1;  Ai(t2, t2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
    }
    
    
    void createRightEConstraints(std::vector<Matrix12> & As)
    {
        
        size_t e1=0, e4=e1+1, e7=e4+1, e2=e7+1,e5=e2+1,e8=e5+1,e3=e8+1,e6=e3+1,e9=e6+1; 
        size_t q1=e9+1, q2=q1+1, q3=q2+1; 
        
        /** Left formulation **/
        Matrix12 Ai; 
        Ai.setZero(); 
        Ai(q1,q1)=1;
        Ai(q2,q2)=1;
        Ai(q3,q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero(); Ai(e1,e1)=1; Ai(e4,e4)=1; Ai(e7,e7)=1; Ai(q2, q2)=-1; Ai(q3, q3)=-1;   
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero(); Ai(e2,e2)=1; Ai(e5,e5)=1; Ai(e8,e8)=1; Ai(q1, q1)=-1; Ai(q3, q3)=-1;   
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero(); Ai(e3,e3)=1; Ai(e6,e6)=1; Ai(e9,e9)=1; Ai(q1, q1)=-1; Ai(q2, q2)=-1;   
        As.push_back(0.5 * (Ai + Ai.transpose()));           
                   
        Ai.setZero(); Ai(e1,e2)=1; Ai(e4,e5)=1; Ai(e7,e8)=1; Ai(q1, q2)=1;                  
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero(); Ai(e1,e3)=1; Ai(e4,e6)=1; Ai(e7,e9)=1; Ai(q1, q3)=1;                
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero(); Ai(e2,e3)=1; Ai(e5,e6)=1; Ai(e8,e9)=1; Ai(q2, q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));


    }
    
    
    void createBothEConstraints(std::vector<Matrix15> & As)
    {
        
        size_t e1=0, e4=e1+1, e7=e4+1, e2=e7+1,e5=e2+1,e8=e5+1,e3=e8+1,e6=e3+1,e9=e6+1; 
        size_t t1=e9+1, t2=t1+1, t3=t2+1;
        size_t q1=t3+1, q2=q1+1, q3=q2+1; 
        
        /** Left formulation **/
        Matrix15 Ai; 
        Ai.setZero(); 
        Ai(t1,t1)=1;
        Ai(t2,t2)=1;
        Ai(t3,t3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero(); 
        Ai(q1,q1)=1;
        Ai(q2,q2)=1;
        Ai(q3,q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        /** Left formulation **/
        /*
        Ai.setZero();         
        Ai(e1,e1)=1;   Ai(e2,e2)=1;Ai(e3,e3)=1;
        Ai(t2, t2)=-1; Ai(t3, t3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  
        */
        
        Ai.setZero();  Ai(e1,e4)=1; Ai(e2,e5)=1;Ai(e3,e6)=1;
        Ai(t1, t2)=1;           
        As.push_back(0.5 * (Ai + Ai.transpose()));      
                   
        Ai.setZero();  Ai(e1,e7)=1; Ai(e2,e8)=1;Ai(e3,e9)=1;
        Ai(t1, t3)=1;       
        As.push_back(0.5 * (Ai + Ai.transpose()));          
                   
        Ai.setZero();  Ai(e4,e4)=1; Ai(e5,e5)=1;Ai(e6,e6)=1;
        Ai(t1, t1)=-1;  Ai(t3, t3)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero();  Ai(e4,e7)=1; Ai(e5,e8)=1;Ai(e6,e9)=1;
        Ai(t2, t3)=1;     
        As.push_back(0.5 * (Ai + Ai.transpose()));            
                   
        Ai.setZero();  Ai(e7,e7)=1; Ai(e8,e8)=1;Ai(e9,e9)=1;
        Ai(t1, t1)=-1;  Ai(t2, t2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        
        /** Right formulation **/
        Ai.setZero(); Ai(e1,e1)=1; Ai(e4,e4)=1; Ai(e7,e7)=1; Ai(q2, q2)=-1; Ai(q3, q3)=-1;   
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero(); Ai(e2,e2)=1; Ai(e5,e5)=1; Ai(e8,e8)=1; Ai(q1, q1)=-1; Ai(q3, q3)=-1;   
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero(); Ai(e3,e3)=1; Ai(e6,e6)=1; Ai(e9,e9)=1; Ai(q1, q1)=-1; Ai(q2, q2)=-1;   
        As.push_back(0.5 * (Ai + Ai.transpose()));           
                   
        Ai.setZero(); Ai(e1,e2)=1; Ai(e4,e5)=1; Ai(e7,e8)=1; Ai(q1, q2)=1;                  
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero(); Ai(e1,e3)=1; Ai(e4,e6)=1; Ai(e7,e9)=1; Ai(q1, q3)=1;                
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero(); Ai(e2,e3)=1; Ai(e5,e6)=1; Ai(e8,e9)=1; Ai(q2, q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

    }    

    
    void createAdjugateEConstraints(std::vector<Matrix15> & As)
    {
        
        size_t e1=0, e4=e1+1, e7=e4+1, e2=e7+1,e5=e2+1,e8=e5+1,e3=e8+1,e6=e3+1,e9=e6+1; 
        size_t t1=e9+1, t2=t1+1, t3=t2+1;
        size_t q1=t3+1, q2=q1+1, q3=q2+1; 
        
        /** Left formulation **/
        Matrix15 Ai; 
        Ai.setZero(); 
        Ai(t1,t1)=1;
        Ai(t2,t2)=1;
        Ai(t3,t3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero(); 
        Ai(q1,q1)=1;
        Ai(q2,q2)=1;
        Ai(q3,q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        /** Left formulation **/
        /*
        Ai.setZero();         
        Ai(e1,e1)=1;   Ai(e2,e2)=1;Ai(e3,e3)=1;
        Ai(t2, t2)=-1; Ai(t3, t3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  
        */
        
        Ai.setZero();  Ai(e1,e4)=1; Ai(e2,e5)=1;Ai(e3,e6)=1;
        Ai(t1, t2)=1;           
        As.push_back(0.5 * (Ai + Ai.transpose()));      
                   
        Ai.setZero();  Ai(e1,e7)=1; Ai(e2,e8)=1;Ai(e3,e9)=1;
        Ai(t1, t3)=1;       
        As.push_back(0.5 * (Ai + Ai.transpose()));          
                   
        Ai.setZero();  Ai(e4,e4)=1; Ai(e5,e5)=1;Ai(e6,e6)=1;
        Ai(t1, t1)=-1;  Ai(t3, t3)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero();  Ai(e4,e7)=1; Ai(e5,e8)=1;Ai(e6,e9)=1;
        Ai(t2, t3)=1;     
        As.push_back(0.5 * (Ai + Ai.transpose()));            
                   
        Ai.setZero();  Ai(e7,e7)=1; Ai(e8,e8)=1;Ai(e9,e9)=1;
        Ai(t1, t1)=-1;  Ai(t2, t2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        
        /** Right formulation **/
        /*
        Ai.setZero(); Ai(e1,e1)=1; Ai(e4,e4)=1; Ai(e7,e7)=1; Ai(q2, q2)=-1; Ai(q3, q3)=-1;   
        As.push_back(0.5 * (Ai + Ai.transpose()));
        */
        
        Ai.setZero(); Ai(e2,e2)=1; Ai(e5,e5)=1; Ai(e8,e8)=1; Ai(q1, q1)=-1; Ai(q3, q3)=-1;   
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero(); Ai(e3,e3)=1; Ai(e6,e6)=1; Ai(e9,e9)=1; Ai(q1, q1)=-1; Ai(q2, q2)=-1;   
        As.push_back(0.5 * (Ai + Ai.transpose()));           
                   
        Ai.setZero(); Ai(e1,e2)=1; Ai(e4,e5)=1; Ai(e7,e8)=1; Ai(q1, q2)=1;                  
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero(); Ai(e1,e3)=1; Ai(e4,e6)=1; Ai(e7,e9)=1; Ai(q1, q3)=1;                
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero(); Ai(e2,e3)=1; Ai(e5,e6)=1; Ai(e8,e9)=1; Ai(q2, q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        
        /** Adjugate constraints **/
        Ai.setZero(); 
        Ai(e1,e1)=1;Ai(e2,e2)=1;Ai(e3,e3)=1;Ai(e4,e4)=1;Ai(e5,e5)=1;
        Ai(e6,e6)=1;Ai(e7,e7)=1;Ai(e8,e8)=1;Ai(e9,e9)=1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));   
                   
        // leftnullspace: E^T t = 0
        Ai.setZero(); Ai(e1,t1)=1;Ai(e4,t2)=1;Ai(e7,t3)=1;   
        As.push_back(0.5 * (Ai + Ai.transpose()));
                 
        Ai.setZero(); Ai(e2,t1)=1;Ai(e5,t2)=1;Ai(e8,t3)=1;   
        As.push_back(0.5 * (Ai + Ai.transpose()));
                
        Ai.setZero(); Ai(e3,t1)=1;Ai(e6,t2)=1;Ai(e9,t3)=1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
                
        // rightnullspace: E q = 0
        Ai.setZero(); Ai(e1,q1)=1;Ai(e2,q2)=1;Ai(e3,q3)=1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
                     
        Ai.setZero(); Ai(e4,q1)=1;Ai(e5,q2)=1;Ai(e6,q3)=1;  
        As.push_back(0.5 * (Ai + Ai.transpose()));
                    
        Ai.setZero(); Ai(e7,q1)=1;Ai(e8,q2)=1;Ai(e9,q3)=1;     
        As.push_back(0.5 * (Ai + Ai.transpose()));
                 
        //trace(adj(E*E')) == 1
        Ai.setZero();  Ai(e5, e9)=1;Ai(e6, e8)=-1;Ai(t1,q1)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero();  Ai(e3, e8)=1;Ai(e2, e9)=-1;Ai(t2,q1)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
                   
        Ai.setZero();  Ai(e2, e6)=1;Ai(e3, e5)=-1;Ai(t3,q1)=-1;    
        As.push_back(0.5 * (Ai + Ai.transpose()));
                 
        Ai.setZero();  Ai(e6, e7)=1;Ai(e4, e9)=-1;Ai(t1,q2)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
                    
        Ai.setZero();  Ai(e1, e9)=1;Ai(e3, e7)=-1;Ai(t2,q2)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
                    
        Ai.setZero();  Ai(e3, e4)=1;Ai(e1, e6)=-1;Ai(t3,q2)=-1;    
        As.push_back(0.5 * (Ai + Ai.transpose()));
                 
        Ai.setZero();  Ai(e4, e8)=1;Ai(e5, e7)=-1;Ai(t1,q3)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
                    
        Ai.setZero();  Ai(e2, e7)=1;Ai(e1, e8)=-1;Ai(t2,q3)=-1;  
        As.push_back(0.5 * (Ai + Ai.transpose()));
                  
        Ai.setZero();  Ai(e1, e5)=1;Ai(e2, e4)=-1;Ai(t3,q3)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
                 
    }    

}  // end of namespace Essential
