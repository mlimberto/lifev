// Matteo Limberto
// test of mesh loading for the heart problem

// Include Epetra tools
#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// Include Teuchos tools
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"


// Include general libraries
#include <fstream>
#include <string>

// Include LifeV core and ETA
#include <lifev/core/LifeV.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/filter/ExporterHDF5.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/core/solver/ADRAssembler.hpp>

#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>


// Namespaces
using namespace LifeV;

// Function declarations
Real sourceFunction(const Real &t , const Real &x , const Real &y , 
										const Real &z , const ID &);


int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);
  {
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD));

    if (Comm->MyPID()==0)
      std::cout <<"Test mesh, using MPI" << std::endl;


    typedef RegionMesh<LinearTetra>            mesh_Type;
    typedef boost::shared_ptr<mesh_Type>       meshPtr_Type;
    typedef VectorEpetra                        vector_Type;
    typedef boost::shared_ptr<VectorEpetra>    vectorPtr_Type;
    typedef MatrixEpetra<Real>                  matrix_Type;
    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh<LifeV::LinearTetra> > >    filterPtr_Type;
    typedef LifeV::ExporterHDF5< RegionMesh<LinearTetra> >                                  hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                                              hdf5FilterPtr_Type;
    typedef ExporterData<mesh_Type>                           exporterData_Type;
    typedef Exporter< mesh_Type >                             IOFile_Type;
    typedef boost::shared_ptr< IOFile_Type >                  IOFilePtr_Type;
    typedef ExporterHDF5< mesh_Type >                         hdf5IOFile_Type;

    // Get the data file using GetPot
    GetPot command_line (argc,argv);
    const std::string data_file_name = command_line.follow("data",2,"-f","--file");
    GetPot dataFile (data_file_name);


    // Load the mesh
    if (Comm->MyPID()==0)
      std::cout <<"[Loading the mesh]" << std::endl ;

    MeshData meshData;
    meshData.setup(dataFile,"discretization/space");

    boost::shared_ptr<mesh_Type> fullMeshPtr (new mesh_Type(Comm) );
    readMesh(*fullMeshPtr , meshData);

    // Partition the mesh
    boost::shared_ptr<mesh_Type> localMeshPtr;
    {
      MeshPartitioner<mesh_Type> meshPart(fullMeshPtr , Comm) ;
      localMeshPtr = meshPart.meshPartition();
    }


    // FEspace
    typedef FESpace< mesh_Type , MapEpetra > uSpaceStd_Type;
    typedef boost::shared_ptr<uSpaceStd_Type> uSpaceStdPtr_Type;
    typedef ETFESpace < mesh_Type, MapEpetra , 3, 1 > uSpaceETA_Type;
    typedef boost::shared_ptr < uSpaceETA_Type > uSpaceETAPtr_Type;
    typedef FESpace<mesh_Type , MapEpetra>::function_Type function_Type;

    // Define standard FE and ETA spaces
    uSpaceStdPtr_Type uFESpace (new uSpaceStd_Type (localMeshPtr , 
    														dataFile("finite_element/degree","P1") , 1 , Comm));

    uSpaceETAPtr_Type ETuFESpace ( new uSpaceETA_Type ( localMeshPtr , 
    															& (uFESpace->refFE() ) , & (uFESpace->fe().geoMap() ) , Comm ));


    // Create a vector
    vectorPtr_Type out;

    out.reset(new vector_Type (uFESpace->map() , Unique)) ;
    out->zero();




    // Set generic exporter postprocessing
    boost::shared_ptr< ExporterHDF5 <mesh_Type > > exporter;
    
    exporter.reset (new ExporterHDF5<mesh_Type> (dataFile,"heart"));
    exporter->setPostDir("./");
    exporter->setMeshProcId( localMeshPtr , Comm->MyPID());
    exporter->addVariable( ExporterData<mesh_Type>::ScalarField , "example" ,
        											uFESpace , out , UInt(0) );
 
    exporter->postProcess(0);

    exporter->closeFile();

    MPI_Barrier(MPI_COMM_WORLD);



  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}


// source function definition
Real sourceFunction(const Real &t , const Real &x , const Real &y , 
										const Real &z , const ID&)
{
	return 1.0;
}






















