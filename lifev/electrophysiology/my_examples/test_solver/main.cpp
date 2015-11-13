// Matteo Limberto
// test of mesh loading for the heart problem

// Include general libraries
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

#include <lifev/core/filter/ExporterHDF5.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>


#include <lifev/core/solver/ADRAssembler.hpp>

#include <lifev/eta/expression/BuildGraph.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>

#include "lifev/electrophysiology/solver/InverseETAEllipticSolver.hpp"

#include "myFunctor.hpp"

// Namespaces
using namespace LifeV;

// Type definitions
typedef RegionMesh<LinearTetra> mesh_Type;
typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

typedef FESpace<mesh_Type , MapEpetra >           uSpaceStd_Type;
typedef boost::shared_ptr<uSpaceStd_Type>         uSpaceStdPtr_Type;
typedef ETFESpace<mesh_Type,MapEpetra , 3 , 1>    uSpaceETA_Type;
typedef boost::shared_ptr<uSpaceETA_Type>         uSpaceETAPtr_Type;

typedef uSpaceStd_Type::function_Type             function_Type;

typedef Epetra_FECrsGraph                         graph_Type;
typedef boost::shared_ptr<graph_Type>             graphPtr_Type;
typedef MatrixEpetra<Real>                        matrix_Type;
typedef boost::shared_ptr<matrix_Type>            matrixPtr_Type;

typedef VectorEpetra                              vector_Type;
typedef boost::shared_ptr<vector_Type>            vectorPtr_Type;

typedef LinearSolver::SolverType                  solver_Type;
typedef LifeV::Preconditioner                     basePrec_Type;
typedef boost::shared_ptr<basePrec_Type>          basePrecPtr_Type;
typedef PreconditionerIfpack                      prec_Type;
typedef boost::shared_ptr<prec_Type>              precPtr_Type;



// sourceFunction definition
Real sourceFunction (const Real & t , const Real & x , const Real & y , 
                     const Real & z , const ID&  )
{
    return 3 * M_PI * M_PI
             * sin (M_PI * x ) * sin (M_PI * y ) * sin (M_PI * z ) ; 
}

Real zeroFunction (const Real & t , const Real & x , const Real & y , 
                   const Real & z , const ID&  )
{
    return 0;
}

// Main

int main(int argc, char *argv[])
{
  #ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD));
  #else
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm) ;
  #endif

    typedef RegionMesh<LinearTetra>             mesh_Type;

    if (Comm->MyPID()==0)
        std::cout << "Testing solver class..." << std::endl;

    bool verbose = true;

    // Load data
    GetPot command_line(argc,argv);
    const std::string dataFileName = command_line.follow("data",2,"-f","--file");
    GetPot dataFile(dataFileName);


    // Build mesh

    if(verbose && Comm->MyPID() == 0 )
        std::cout << "[Loading mesh]" << std::endl;

    meshPtr_Type fullMeshPtr(new mesh_Type(Comm));

    regularMesh3D( *fullMeshPtr , 0 ,
                   dataFile("mesh/nx",15) , dataFile("mesh/ny",15) , dataFile("mesh/nz",15) ,
                   dataFile("mesh/verbose",false) ,
                   2.0 , 2.0 , 2.0 , -1.0 , -1.0 , -1.0 ) ;

    const UInt overlap(dataFile("mesh/overlap",0));
  
    meshPtr_Type localMeshPtr;
    {
        MeshPartitioner<mesh_Type> meshPart;
        if (overlap)
            meshPart.setPartitionOverlap(overlap);
        meshPart.doPartition(fullMeshPtr , Comm);
        localMeshPtr = meshPart.meshPartition();
    }

    fullMeshPtr.reset();

    // manage bc
    const int BACK = 1;
    const int FRONT = 2;
    const int LEFT = 3;
    const int RIGHT = 4;
    const int TOP = 5;
    const int BOTTOM = 6;

    BCHandler bcHandler;

    BCFunctionBase ZeroBC (zeroFunction) ;

    bcHandler.addBC("Back", BACK , Essential , Scalar , ZeroBC , 1) ;
    bcHandler.addBC("Left" , LEFT , Essential , Scalar , ZeroBC , 1);
    bcHandler.addBC("Top" , TOP , Essential , Scalar , ZeroBC , 1);
    bcHandler.addBC("Front" , FRONT , Essential , Scalar , ZeroBC , 1);
    bcHandler.addBC("Right" , RIGHT , Essential , Scalar , ZeroBC , 1);
    bcHandler.addBC("Bottom" , BOTTOM , Essential , Scalar , ZeroBC , 1);

    // Instantiate class
    InverseETAEllipticSolver<mesh_Type> solver( dataFile ,
                                                localMeshPtr ,
                                                "SolverParamList2.xml" ,
                                                bcHandler
                                                );

    solver.solveFwd();


//  // Assemble matrix and rhs

//  if(verbose)
//    std::cout << "[Building graph and matrix]" << std::endl;
  
//  graphPtr_Type systemGraph;
//  matrixPtr_Type systemMatrix;

//  if (overlap)
//  {
//    systemGraph.reset( new graph_Type( Copy , * (uFESpace->map().map(Unique)),
//                                        50,true ) );
//  }
//  else
//  {
//    systemGraph.reset( new graph_Type( Copy , * (uFESpace->map().map(Unique)),
//                                        50 ) );
//  }

//  {
//    using namespace ExpressionAssembly ;

//    buildGraph (
//              elements(localMeshPtr) ,
//              uFESpace->qr() ,
//              ETuFESpace,
//              ETuFESpace,
//              dot(grad(phi_i) , grad(phi_j) )
//              )
//              >> systemGraph ;
//  }

//  systemGraph->GlobalAssemble();


//  if(overlap)
//  {
//    systemMatrix.reset( new matrix_Type ( ETuFESpace->map() , *systemGraph , true));
//  }
//  else
//  {
//    systemMatrix.reset( new matrix_Type ( ETuFESpace->map() , *systemGraph ));
//  }

//  systemMatrix->zero();

//  {
//    using namespace ExpressionAssembly;

//    integrate (
//      elements(localMeshPtr) ,
//      uFESpace->qr() ,
//      ETuFESpace,
//      ETuFESpace,
//      dot (grad ( phi_i) , grad(phi_j) )
//    )
//    >> systemMatrix ;
//  }


//  // Vectors
//  vectorPtr_Type fwdRhs;
//  vectorPtr_Type fwdSol;

//  if(overlap)
//  {
//    fwdRhs.reset(new vector_Type( uFESpace->map() , Unique, Zero));
//    fwdSol.reset(new vector_Type( uFESpace->map() , Unique, Zero));
//  }
//  else
//  {
//    fwdRhs.reset(new vector_Type( uFESpace->map() , Unique));
//    fwdSol.reset(new vector_Type( uFESpace->map() , Unique));
//  }

//  fwdRhs->zero();
//  fwdSol->zero();

//  boost::shared_ptr<myFunctor<Real> > mySourceFunctor( new myFunctor<Real>(sourceFunction));

//  {
//    using namespace ExpressionAssembly;

//    integrate(
//      elements(localMeshPtr) ,
//      uFESpace->qr() ,
//      ETuFESpace ,
//      eval(mySourceFunctor,X) * phi_i
//      )
//      >> fwdRhs ;

//  }

//  // Boundary conditions
//  const int BACK = 1;
//  const int FRONT = 2;
//  const int LEFT = 3;
//  const int RIGHT = 4;
//  const int TOP = 5;
//  const int BOTTOM = 6;

//  BCHandler bcHandler;

//  BCFunctionBase ZeroBC (zeroFunction) ;

//  bcHandler.addBC("Back", BACK , Essential , Scalar , ZeroBC , 1) ;
//  bcHandler.addBC("Left" , LEFT , Essential , Scalar , ZeroBC , 1);
//  bcHandler.addBC("Top" , TOP , Essential , Scalar , ZeroBC , 1);
//  bcHandler.addBC("Front" , FRONT , Essential , Scalar , ZeroBC , 1);
//  bcHandler.addBC("Right" , RIGHT , Essential , Scalar , ZeroBC , 1);
//  bcHandler.addBC("Bottom" , BOTTOM , Essential , Scalar , ZeroBC , 1);

//  bcHandler.bcUpdate( *uFESpace->mesh() , uFESpace->feBd() , uFESpace->dof() );
//  bcManage (*systemMatrix , *fwdRhs , *uFESpace->mesh() , uFESpace->dof()  ,
//            bcHandler , uFESpace->feBd() , 1.0 , 0.0 ) ;

//  systemMatrix->globalAssemble();
//  fwdRhs->globalAssemble();

//  // Setting solver
//  LinearSolver linearSolver( Comm ) ;
//  linearSolver.setOperator(systemMatrix) ;

//  Teuchos::RCP < Teuchos::ParameterList > aztecList = Teuchos::rcp (
//    new Teuchos::ParameterList) ;

//  aztecList = Teuchos::getParametersFromXmlFile ("SolverParamList2.xml");

//  linearSolver.setParameters(*aztecList) ;

//  // Setting preconditioner
//  prec_Type* precRawPtr( new prec_Type) ;
//  precRawPtr->setDataFromGetPot(dataFile , "prec");

//  basePrecPtr_Type precPtr;
//  precPtr.reset(precRawPtr) ;

//  linearSolver.setPreconditioner(precPtr) ;

//  // Set rhs
//  linearSolver.setRightHandSide (fwdRhs) ;

//  if (verbose)
//    std::cout << "Solving the problem... " << std::flush;

//  linearSolver.solve(fwdSol) ;

//  if (verbose)
//    std::cout << "done!" << std::endl << std::flush ;


//  // Export results
//  if (verbose)
//    std::cout << "[Export output]" << std::endl;

//  ExporterHDF5 <mesh_Type > exporter (dataFile , "fwd");
//  exporter.setMeshProcId (localMeshPtr , Comm->MyPID());
//  exporter.setPrefix("laplace");
//  exporter.setPostDir("./");
//  exporter.addVariable( ExporterData<mesh_Type>::ScalarField , "temperature" ,
//                        uFESpace , fwdSol , UInt(0));

//  exporter.postProcess(0);
//  exporter.closeFile();


  #ifdef HAVE_MPI
    MPI_Finalize();
  #endif
  return EXIT_SUCCESS;
}










































