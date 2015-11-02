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

// Include general libraries
#include <fstream>
#include <string>

// Include LifeV core
#include <lifev/core/array/VectorSmall.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
// #include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

// #include <lifev/electrophysiology/solver/IonicModels/IonicAlievPanfilov.hpp>
// #include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

//#include <lifev/electrophysiology/examples/example_ECG/Norm.hpp>
#include <lifev/core/solver/ADRAssembler.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>


// Namespaces
using namespace LifeV;


int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);
	{
		boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD));

		if (Comm->MyPID()==0)
			std::cout <<"Test mesh, using MPI" << std::endl;


    typedef RegionMesh<LinearTetra>            mesh_Type;
    typedef boost::shared_ptr<mesh_Type>       meshPtr_Type;
    typedef boost::shared_ptr<VectorEpetra>    vectorPtr_Type;
    typedef FESpace< mesh_Type, MapEpetra >    feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>    feSpacePtr_Type;
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >   function_Type;
    typedef VectorEpetra                        vector_Type;
    typedef MatrixEpetra<Real>                  matrix_Type;
    typedef LifeV::Preconditioner               basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>    basePrecPtr_Type;
    typedef LifeV::PreconditionerML             prec_Type;
    typedef boost::shared_ptr<prec_Type>        precPtr_Type;
    // typedef ElectroETAMonodomainSolver< mesh_Type, IonicMinimalModel >        monodomainSolver_Type;
    // typedef boost::shared_ptr< monodomainSolver_Type >                        monodomainSolverPtr_Type;
    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh<LifeV::LinearTetra> > >    filterPtr_Type;
    typedef LifeV::ExporterHDF5< RegionMesh<LinearTetra> >                                  hdf5Filter_Type;
    typedef boost::shared_ptr<hdf5Filter_Type>                                              hdf5FilterPtr_Type;
    typedef ExporterData<mesh_Type>                           exporterData_Type;
    typedef Exporter< mesh_Type >                             IOFile_Type;
    typedef boost::shared_ptr< IOFile_Type >                  IOFilePtr_Type;
    typedef ExporterHDF5< mesh_Type >                         hdf5IOFile_Type;

		if (Comm->MyPID()==0)
			std::cout <<"Importing parameter list... " ;

		Teuchos::ParameterList monodomainList = *(Teuchos::getParametersFromXmlFile ("MonodomainSolverParamList.xml"));

		if (Comm->MyPID()==0)
			std::cout <<"done!" << std::endl ;

		// Load the mesh
		std::string meshName = monodomainList.get ("mesh_name" , "lid16.mesh");
		std::string meshPath = monodomainList.get ("mesh_path" , "./");

		if (Comm->MyPID()==0)
			std::cout <<"[Loading the mesh]" << std::endl ;

		meshPtr_Type fullMeshPtr (new mesh_Type (Comm));

		std::vector<Real> meshDim(3,0);
		meshDim[0] = monodomainList.get("meshDim_X", 10);
		meshDim[1] = monodomainList.get("meshDim_Y", 10);
		meshDim[2] = monodomainList.get("meshDim_Z", 10);
		std::vector<Real> domain(3,0);
		domain[0] = monodomainList.get("domain_X", 1.);
		domain[1] = monodomainList.get("domain_Y", 1.);
		domain[2] = monodomainList.get("domain_Z", 1.);

		regularMesh3D ( *fullMeshPtr,
										1,
										meshDim[0],meshDim[1],meshDim[2],
										false,
										domain[0],domain[1],domain[2],
										0.0 , 0.0 , 0.0 );

    if ( Comm->MyPID() == 0 )
    	std::cout << "Mesh size: " << MeshUtility::MeshStatistics::computeSize(* fullMeshPtr).maxH << std::endl;

		if ( Comm->MyPID() == 0 )
    	std::cout << "Partitioning the mesh ... " << std::endl;
    
    meshPtr_Type meshPtr;
    {
    	MeshPartitioner<mesh_Type> meshPart(fullMeshPtr,Comm);
    	meshPtr = meshPart.meshPartition();
    }
    fullMeshPtr.reset(); // freeing global mesh





		// Setup the preconditioner
		GetPot command_line (argc,argv);
		const std::string data_file_name = command_line.follow("data",2,"-f","--file");
		GetPot dataFile (data_file_name);


	}	

	MPI_Finalize();
  return EXIT_SUCCESS;
}
























