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


//! Scalar Functor Handler class
//! It would be useful to use boost functions instead of functions pointers
//! so that it becomes less prone to errors
template <typename Return_Type>
class InverseETAScalarFunctorHandler
{
public :

    typedef Return_Type return_Type;

    InverseETAScalarFunctorHandler( return_Type (*newFunction)(const Real & ,const Real &x , const Real &y ,
                                                               const Real &z , const ID & )) :
        myFunction(newFunction)
    {

    }

    virtual ~InverseETAScalarFunctorHandler()
    {

    }

    return_Type operator() (VectorSmall<3> spaceCoord)
    {
        Real x = spaceCoord[0];
        Real y = spaceCoord[1];
        Real z = spaceCoord[2];

        return myFunction(0,x,y,z,0);
    }

private :

    return_Type (*myFunction) (const Real &, const Real &x , const Real &y , const Real &z , const ID &) ;

};

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

    //! MultiLevel preconditioner
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

    //! Boost function
    typedef Real (*function_Type)(const Real &, const Real &x , const Real &y , const Real &z , const ID &);

    typedef boost::shared_ptr<function_Type>        functionPtr_Type;

    //! Functor handler

    typedef boost::shared_ptr<InverseETAScalarFunctorHandler<Real>>  functor_Type;


    //@}

    //! @name Constructors and Destructor
    //@{

    //! Constructor
    /*!
     * @param GetPot dataFile
     * @param std::shared_ptr<Mesh> Pointer to the partitioned mesh
     */
    InverseETAEllipticSolver(GetPot& dataFile,
                             meshPtr_Type meshPtr ,
                             const std::string & solverParamFile,
                             const BCHandler &bcHandler,
                             const function_Type torsoLocation);


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
    virtual ~InverseETAEllipticSolver() { }

    //@}



    //! @name methods
    //@{

    //! Setup forward global matrix
    void setupFwdMatrix();

    //! Setup forward rhs
    void setupFwdRhs();

    //! Setup forward boundary conditions
    void setupFwdBC(BCHandler& bcHandler );

    //! Solve forward problem
    void solveFwd();

#ifdef HAVE_HDF5
    //! Solve forward problem and export solution to given exporter
    void solveFwd( hdf5IOFile_Type& exporter );
#endif

    //@}


    //! @name SET methods
    //@{

    inline void setVerbosity(bool verbose) { M_verbose = verbose ;}

    //@}

    //! @name GET methods
    //@{

    //@}

    //! @name TEST methods
    //@{

    struct exampleRhs {

        typedef Real return_Type;

        Real operator() ( VectorSmall<3> spaceCoord)
        {
            // This solution to store functions is
            // hopefully only temporary

            Real x = spaceCoord[0] ;
            Real y = spaceCoord[1] ;
            Real z = spaceCoord[2] ;

            return 3 * M_PI * M_PI
                     * sin (M_PI * x ) * sin (M_PI * y ) * sin (M_PI * z ) ;
        }
    } ;


    //@}



private:

    // GetPot data file
    GetPot M_dataFile;

    // Communicator
    commPtr_Type M_commPtr;

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

    // rhs
    vectorPtr_Type M_fwdRhsPtr;

    // solution
    vectorPtr_Type M_fwdSolPtr;

    // preconditioner pointers
    prec_Type * M_precRawPtr;
    basePrecPtr_Type M_precPtr;

    // linear solver
    linearSolverPtr_Type M_linearSolverPtr;

    // identity matrix
    matrixSmall_Type M_identity;

    // diffusion coefficients
    VectorSmall<3> M_diffusionTensorTorso;
    VectorSmall<3> M_intDiffusionTensorHeartHealthy;
    VectorSmall<3> M_extDiffusionTensorHeartHealthy;
    VectorSmall<3> M_extDiffusionTensorHeartIsch;

    // function handlers for torso location
    functor_Type   M_torsoLocationFunctor;

    // bc handler (in theory not needed !! )
    BCHandler M_bcHandler;

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

//! Constructor with dataFile and partitioned mesh
template<typename Mesh>
InverseETAEllipticSolver<Mesh>::InverseETAEllipticSolver(GetPot& dataFile,
                                                         meshPtr_Type meshPtr ,
                                                         const std::string & solverParamFile ,
                                                         const BCHandler & bcHandler ,
                                                         const function_Type torsoLocation) :
    M_dataFile ( dataFile ) ,
    M_localMeshPtr ( meshPtr ) ,
    M_commPtr (meshPtr->comm() ) ,
    M_fullMeshPtr( new mesh_Type(M_commPtr) ) ,
    M_stiffMatrixPtr ( ) ,
    M_fwdRhsPtr ( ) ,
    M_fwdSolPtr ( ) ,
    M_precRawPtr( new prec_Type) ,
    M_linearSolverPtr( new linearSolver_Type(M_commPtr) ),
    M_bcHandler( bcHandler ) ,
    M_torsoLocationFunctor( new InverseETAScalarFunctorHandler<Real>( torsoLocation ) ) ,
    M_verbose ( true )
{
    if (M_verbose && M_commPtr->MyPID() == 0)
        std::cout << "[Solver class instantiation]" << std::endl;

    // Read physical parameters
    if (M_verbose && M_commPtr->MyPID() == 0)
        std::cout << "  [Importing physical parameters]" << std::endl;

    M_diffusionTensorTorso(0)           = dataFile("parameters/conductivity_torso", 0.0);
    M_intDiffusionTensorHeartHealthy(0) = dataFile("parameters/conductivity_heart_healthy_int", 0.0);
    M_extDiffusionTensorHeartHealthy(0) = dataFile("parameters/conductivity_heart_healthy_ext", 0.0);
    M_extDiffusionTensorHeartIsch(0)    = dataFile("parameters/conductivity_heart_ischemic_ext", 0.0);


    // Setup identity matrix
    M_identity (0, 0) = 1.0;
    M_identity (0, 1) = 0.0;
    M_identity (0, 2) = 0.0;
    M_identity (1, 0) = 0.0;
    M_identity (1, 1) = 1.0;
    M_identity (1, 2) = 0.0;
    M_identity (2, 0) = 0.0;
    M_identity (2, 1) = 0.0;
    M_identity (2, 2) = 1.0;

    // Define finite elements spaces
    if (M_verbose && M_commPtr->MyPID() == 0)
        std::cout << "  [Setting up FE spaces]" << std::endl;

    M_FESpacePtr.reset( new FESpace_Type( M_localMeshPtr ,
                                     M_dataFile("finite_element/degree","P1") ,
                                     1 ,
                                     M_commPtr) );

    M_ETFESpacePtr.reset( new ETFESpace_Type ( M_localMeshPtr ,
                                          & (M_FESpacePtr->refFE() ) ,
                                          & (M_FESpacePtr->fe().geoMap()) ,
                                          M_commPtr ) ) ;


    // Setting solver parameters
    if (M_verbose && M_commPtr->MyPID() == 0)
        std::cout << "  [Setting up solver parameters]" << std::endl;

    Teuchos::RCP < Teuchos::ParameterList > aztecList = Teuchos::rcp (
                new Teuchos::ParameterList ) ;
    aztecList = Teuchos::getParametersFromXmlFile(solverParamFile);

    M_linearSolverPtr->setParameters(*aztecList);

    // Setting preconditioner
    if (M_verbose && M_commPtr->MyPID() == 0)
        std::cout << "  [Setting up preconditioner]" << std::endl;

    M_precRawPtr->setDataFromGetPot(dataFile , "prec");

    M_precPtr.reset (M_precRawPtr);

    M_linearSolverPtr->setPreconditioner(M_precPtr);

}


//! Constructor with datafile and name and path of mesh
template<typename Mesh>
InverseETAEllipticSolver<Mesh>::InverseETAEllipticSolver(std::string meshName,
                                                         std::string meshPath ,
                                                         GetPot& dataFile) :
    M_dataFile ( dataFile ) ,
    M_localMeshPtr ( ) ,
    M_commPtr () ,
    M_fullMeshPtr(  ) ,
    M_stiffMatrixPtr ( ) ,
    M_fwdRhsPtr ( ) ,
    M_fwdSolPtr ( ) ,
    M_verbose ( true )
{





}


//! Constructor with communicator
template<typename Mesh>
InverseETAEllipticSolver<Mesh>::InverseETAEllipticSolver(std::string meshName,
                                                         std::string meshPath ,
                                                         GetPot& dataFile ,
                                                         commPtr_Type comm)
{

}


//
//
// ===================================================
//! General methods
// ===================================================
//
//

//! Setup global matrix for forward problem
template<typename Mesh>
void InverseETAEllipticSolver<Mesh>::setupFwdMatrix()
{
    if (M_verbose && M_commPtr->MyPID() == 0) {
        std::cout << "[Setting up forward matrix]" << std::endl;
        std::cout << "  Integration ..." << std::endl; }


    // Initialize

    M_stiffMatrixPtr.reset( new matrix_Type( M_ETFESpacePtr->map() ) );

    M_stiffMatrixPtr->zero();


    // Integrate

    {
        using namespace ExpressionAssembly;


        auto I = value(M_identity) ;

        Real torsoSigma = M_diffusionTensorTorso(0);                    // M_0

        Real healthyHeartSigma = M_extDiffusionTensorHeartHealthy(0)    // M_e + M_i
                               + M_intDiffusionTensorHeartHealthy(0);

        auto torsoD = value(torsoSigma) * I ;

        auto heartD = value(healthyHeartSigma) * I ;

        auto D = torsoD * eval( M_torsoLocationFunctor , X)
                + heartD * (1 - eval( M_torsoLocationFunctor ,X)) ;


        integrate( elements(M_localMeshPtr) , M_FESpacePtr->qr() ,
                   M_ETFESpacePtr , M_ETFESpacePtr ,
                   dot ( D * grad(phi_i) , grad(phi_j) ) )
                >> M_stiffMatrixPtr;
    }

    // Remark : matrix will be globally assembled after BCs are imposed
}

//! Setup global rhs for forward problem
template<typename Mesh>
void InverseETAEllipticSolver<Mesh>::setupFwdRhs()
{
    if (M_verbose && M_commPtr->MyPID() == 0)
        std::cout << "[Setting up forward right-hand-side]" << std::endl;

    M_fwdRhsPtr.reset (new vector_Type( M_ETFESpacePtr->map() , Unique) );

    M_fwdRhsPtr->zero();

    {
        using namespace ExpressionAssembly;

//        functionPtr_Type rhsFunctorPtr( new function_Type(exampleRhs() ) ) ;

        boost::shared_ptr<exampleRhs> rhsFunctorPtr(new exampleRhs() ) ;

        integrate( elements(M_localMeshPtr) ,
                   M_FESpacePtr->qr() ,
                   M_ETFESpacePtr ,
                   eval(rhsFunctorPtr , X) * phi_i )
                >> M_fwdRhsPtr ;
    }

    // Remark : vector will be globally assembled after BCs are imposed

}

//! Setup boundary conditions for the forward problem
template<typename Mesh>
void InverseETAEllipticSolver<Mesh>::setupFwdBC(BCHandler &bcHandler)
{
    if (M_verbose && M_commPtr->MyPID() == 0)
        std::cout << "[Setting up forward boundary conditions]" << std::endl;

    bcHandler.bcUpdate( *M_FESpacePtr->mesh() ,
                        M_FESpacePtr->feBd() ,
                        M_FESpacePtr->dof() );


    bcManage( *M_stiffMatrixPtr , *M_fwdRhsPtr ,
              *M_FESpacePtr->mesh() , M_FESpacePtr->dof() ,
              bcHandler , M_FESpacePtr->feBd() ,
              1.0 , 0.0 ) ;

    }

//! Solve forward problem
template<typename Mesh>
void InverseETAEllipticSolver<Mesh>::solveFwd()
{
    // Setup matrix
    setupFwdMatrix();

    // Setup rhs
    setupFwdRhs();

    // Impose boundary conditions
    setupFwdBC(M_bcHandler);

    // Assemble matrix and vector
    if (M_verbose && M_commPtr->MyPID() == 0)
        std::cout << "Global assembly ..." << std::endl;

    M_stiffMatrixPtr->globalAssemble();

//    M_stiffMatrixPtr->spy("stiffness_matrix");

    M_fwdRhsPtr->globalAssemble();

    // Setup linear solver
    M_linearSolverPtr->setOperator(M_stiffMatrixPtr);

    M_linearSolverPtr->setRightHandSide(M_fwdRhsPtr);

    // Initialize the solution vector
    M_fwdSolPtr.reset(new vector_Type(M_FESpacePtr->map() , Unique) );
    M_fwdSolPtr->zero();


    // SOLVE
    if (M_verbose && M_commPtr->MyPID() == 0)
        std::cout << "Solving the problem ..." << std::endl;

    M_linearSolverPtr->solve(M_fwdSolPtr) ;

}

#ifdef HAVE_HDF5
template<typename Mesh>
void InverseETAEllipticSolver<Mesh>::solveFwd( hdf5IOFile_Type& exporter )
{
    // Solve problem

    solveFwd();

    // Export solution

    if (M_verbose && M_commPtr->MyPID() == 0)
        std::cout << "Adding solution to given exporter ..." << std::endl;

    exporter.addVariable( ExporterData<mesh_Type>::ScalarField , "potential" ,
                          M_FESpacePtr , M_fwdSolPtr , UInt(0) ) ;
}
#endif


//
//
// ===================================================
//! Test methods
// ===================================================
//
//





} // namespace LifeV



#endif // _INVERSEETAELLIPTICSOLVER_HPP_H_

















































