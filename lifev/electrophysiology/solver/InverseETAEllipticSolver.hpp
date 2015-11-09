#ifndef _INVERSEETAELLIPTICSOLVER_HPP_H_
#define _INVERSEETAELLIPTICSOLVER_HPP_H_


// Include general libraries
#include <iostream>
#include <fstream>
#include <string>

// Include Teuchos tools
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

// Include Epetra tools
#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// Include LifeV core and ETA
#include <lifev/core/LifeV.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/fem/BCManage.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>

#include <lifev/core/solver/ADRAssembler.hpp>

#include <lifev/eta/expression/BuildGraph.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>

namespace LifeV
{

//! Elliptic Control Problem solver - Class featuring the solver for the stochastic control problem

template <typename Mesh >
class InverseETAEllipticSolver
{
public :

    //! @name Type definitions
    //@{

    //! Mesh
    typedef Mesh                                    mesh_Type;
    typedef boost::shared_ptr<mesh_Type>            meshPtr_Type;

    //! Distributed vectors and matrices
    typedef VectorEpetra                            vector_Type;
    typedef boost::shared_ptr<VectorEpetra>         vectorPtr_Type;

    typedef MatrixEpetra<Real>                      matrix_Type;
    typedef boost::shared_ptr<matrix_Type>          matrixPtr_Type;

    //! 3x3 matrix
    typedef MatrixSmall<3,3>                        matrixSmall_Type;

    //! Epetra Communicator
    typedef Epetra_Comm                             comm_Type;
    typedef boost::shared_ptr<Epetra_Comm>          commPtr_Type;

    //! FE space
    typedef FESpace<mesh_Type, MapEpetra>           FESpace_Type;
    typedef boost::shared_ptr<FESpace_Type>         FESpacePtr_Type;

    //! Expression template FE space
    typedef ETFESpace<mesh_Type,MapEpetra,3,1>      ETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace_Type>       ETFESpacePtr_Type;

    //! Preconditioner
    typedef LifeV::Preconditioner                   basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>        basePrecPtr_Type;
    typedef LifeV::PreconditionerML                 prec_Type;
    typedef boost::shared_ptr<prec_Type>            precPtr_Type;

    //! Linear solver
    typedef LinearSolver                            linearSolver_Type;
    typedef boost::shared_ptr<linearSolver_Type>    linearSolverPtr_Type;

    //! Exporter
    typedef Exporter<mesh_Type>                     IOFile_Type;
    typedef boost::shared_ptr<IOFile_Type>          IOFilePtr_Type;
    typedef ExporterData<mesh_Type>                 IOData_Type;
    typedef boost::shared_ptr<IOData_Type>          IODataPtr_Type;
#ifdef HAVE_HDF5
    typedef ExporterHDF5<mesh_Type>                 hdf5IOFile_Type;
#endif

    //@}

    //! @name Constructors and Destructor
    //@{

    //! Empty Constructor
    /*!
     */
    InverseETAEllipticSolver();


    //! Constructor
    /*!
     * @param GetPot dataFile
     * @param std::shared_ptr<Mesh> Pointer to the partitioned mesh
     */
    InverseETAEllipticSolver(GetPot& dataFile, meshPtr_Type meshPtr);


    //! Constructor
    /*!
     * @param meshName file name of the mesh
     * @param meshPath path to the mesh
     * @param GetPot dataFile
     */
    InverseETAEllipticSolver(std::string meshName, std::string meshPath ,
                             GetPot& dataFile);


    //! Constructor
    /*!
     * @param meshName file name of the mesh
     * @param meshPath path to the mesh
     * @param GetPot dataFile
     * @param std::shared_ptr<Epetra_Comm> Epetra communicator
     */
    InverseETAEllipticSolver(std::string meshName, std::string meshPath ,
                             GetPot& dataFile , commPtr_Type comm);




    //! Destructor
    ~InverseETAEllipticSolver()
    {
    }

    //@}



    //! @name methods
    //@{

    //@}


    //! @name SET methods
    //@{

    inline void setVerbosity(bool verbose) { M_verbose = verbose ;}

    //@}

    //! @name GET methods
    //@{

    //@}



private:

    // GetPot data file
    GetPot M_dataFile;

    // partitioned mesh
    meshPtr_Type M_localMeshPtr;
    // full mesh
    meshPtr_Type M_fullMeshPtr;

    // FE space
    FESpacePtr_Type M_FESpacePtr;
    // ET FE space
    ETFESpacePtr_Type M_ETFESpacePtr;

    // stiffness matrix
    matrixPtr_Type M_stiffMatrixPtr;


    // verbosity
    bool M_verbose;

};
// class InverseETAEllipicSolver

//
// IMPLEMENTATION
//
// ===================================================
//! Constructors
// ===================================================

//! Empty constructor
template<typename Mesh>
InverseETAEllipticSolver<Mesh>::InverseETAEllipticSolver()
{
    std::cout << "Empty constructor is called" << std::endl;
}

//! Constructor with dataFile and partitioned mesh
template<typename Mesh>
InverseETAEllipticSolver<Mesh>::InverseETAEllipticSolver(GetPot& dataFile, meshPtr_Type meshPtr)
{

}


//! Constructor with datafile and name and path of mesh
template<typename Mesh>
InverseETAEllipticSolver<Mesh>::InverseETAEllipticSolver(std::string meshName, std::string meshPath ,
                         GetPot& dataFile)
{

}


//! Constructor with communicator
template<typename Mesh>
InverseETAEllipticSolver<Mesh>::InverseETAEllipticSolver(std::string meshName, std::string meshPath ,
                         GetPot& dataFile , commPtr_Type comm)
{

}


} // namespace LifeV



#endif // _INVERSEETAELLIPTICSOLVER_HPP_H_

















































