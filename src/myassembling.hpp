#ifndef FILE_MYASSEMBLING_HPP
#define FILE_MYASSEMBLING_HPP


#include <comp.hpp>


/*
  Assembling matrix and vector
*/


namespace ngcomp
{
  
  shared_ptr<BaseSparseMatrix>
  MyAssembleMatrix(shared_ptr<FESpace> fes,
                   shared_ptr<BilinearFormIntegrator> bfi);
    
  shared_ptr<BaseVector>
  MyAssembleVector(shared_ptr<FESpace> fes,
                   shared_ptr<LinearFormIntegrator> lfi);

  
  
}

#endif
