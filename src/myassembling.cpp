
/*
  
  Assemble the system matrix

  The input is
  - a finite element space, which provides the basis functions
  - an integrator, which computes the element matrices
  
  The result is a sparse matrix
*/

#include "myassembling.hpp"

namespace ngcomp
{
  shared_ptr<BaseSparseMatrix> MyAssembleMatrix(shared_ptr<FESpace> fes,
                                                shared_ptr<BilinearFormIntegrator> bfi)
  {
    cout << "We assemble matrix" << endl;

    auto ma = fes->GetMeshAccess();
    
    int ndof = fes->GetNDof();
    int ne = ma->GetNE(VOL);

    // we build a sparse matrix
    // the non-zero pattern is given by the connectivity pattern provided by the FESpace
    // el2dof[i] stores dofs of element i
    Table<int> el2dof = fes->CreateDofTable(VOL);
    
    // generate sparse matrix of size ndof x ndof
    // from element-to-dof table for rows and columns
    auto mat = make_shared<SparseMatrix<double>> (ndof, ndof, el2dof, el2dof, false);
    mat -> SetZero();
    
    LocalHeap lh(1000*1000); // reserve 1MB 
    Array<int> dnums;
    
    // loop over all volume elements
    for (int i = 0; i < ne; i++)
      {
        HeapReset hr(lh);  // cleanup heap at end of scope
        ElementId ei(VOL, i);
        
        // let FESpace generate the finite element
        FiniteElement & fel = fes->GetFE (ei, lh);
        
        // the global dof numbers of the element
        fes->GetDofNrs (ei, dnums);
        
        // the mesh knows the geometry of the element
        const ElementTransformation & trafo = ma->GetTrafo (ei, lh);
        
        // compute the element matrix
        FlatMatrix<> elmat (fel.GetNDof(), lh);
        bfi->CalcElementMatrix (fel, trafo, elmat, lh);
        
        mat->AddElementMatrix (dnums, dnums, elmat);
      }

    return mat;
  }



  /*
    Exercise: implement a corresponding function for assembling the right hand side vector
   */
  shared_ptr<BaseVector> MyAssembleVector(shared_ptr<FESpace> fes,
                                          shared_ptr<LinearFormIntegrator> lfi)
  {
    // A VVector (virtual vector) is derived from BaseVector
    shared_ptr<BaseVector> vec = make_shared<VVector<double>> (fes->GetNDof());

    // adding an element vector to the global vector
    FlatVector<double> elvec;
    Array<int> dnums;
    vec->AddIndirect (dnums, elvec);
    
    return vec;
  }
}
