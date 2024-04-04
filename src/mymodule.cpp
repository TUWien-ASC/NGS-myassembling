#include <comp.hpp>
#include <python_comp.hpp>

#include "myintegrator.hpp"
#include "myassembling.hpp"


PYBIND11_MODULE(myassembling, m)
{
  cout << "Loading myassembling library" << endl;

  py::class_<ngfem::MyLaplaceIntegrator, shared_ptr<ngfem::MyLaplaceIntegrator>, ngfem::BilinearFormIntegrator>
    (m, "MyLaplace")
    .def(py::init<shared_ptr<ngfem::CoefficientFunction>>())
    ;
  
  py::class_<ngfem::MySourceIntegrator, shared_ptr<ngfem::MySourceIntegrator>, ngfem::LinearFormIntegrator>
    (m, "MySource")
    .def(py::init<shared_ptr<ngfem::CoefficientFunction>>())
    ;
  
  m.def ("MyAssembleMatrix",
         &ngcomp::MyAssembleMatrix,
         py::arg("fes"),py::arg("integrator"));
}    

