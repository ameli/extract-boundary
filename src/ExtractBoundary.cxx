/*
 * =====================================================================================
 *
 *       Filename:  ExtractBoundary.cxx
 *
 *    Description:  Extract Boundary curve of 2D data
 *
 *        Version:  1.0
 *        Created:  12/09/2012 10:07:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

// =======
// Headers
// =======

#include "ExtractBoundary.h"

// Pipeline
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>

// Computation
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolygon.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkUnsignedIntArray.h>
#include <vtkGenericCell.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>

// ======
// Macros
// ======

#define EPSILON 1e-10
vtkStandardNewMacro(ExtractBoundary);
// vtkCxxRevisionMacro(ExtractBoundary,"$Revision 1.0$");

// ===========
// Constructor
// ===========

ExtractBoundary::ExtractBoundary()
{
    // Default Member Data
    this->SearchBoundaryMode = SEARCH_OUTER_BOUNDARY;
    this->ExtractVelocityBoundary = true;
    this->CalculateNormals = true;
    this->NormalsOutward = false;
    this->BoundaryPointIdsInInputGrid = NULL;

    // Pipeline
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

// ==========
// Destructor
// ==========

ExtractBoundary::~ExtractBoundary()
{
}

// ==========
// Print Self
// ==========

void ExtractBoundary::PrintSelf(ostream &os, vtkIndent indent)
{
    // Superclass member data
    this->Superclass::PrintSelf(os,indent);

    // Class member data
    os << indent << "SearchBoundaryMode: " << this->GetSearchBoundaryModeAsString() << std::endl;
    os << indent << "ExtractVelocityBoundary: " << this->GetExtractVelocityBoundary() << std::endl;
    os << indent << "CalculateNormals: " << this->GetCalculateNormals() << std::endl;
    os << indent << "NormalsOutward: " << this->GetNormalsOutward() << std::endl; 
}

// ===================
// Accessors, Mutators
// ===================

// Cast from int to enum for Paraview Plugin parameters
void ExtractBoundary::SetSearchBoundaryMode(int IntSearchBoundaryMode)
{
    // this->SearchBoundaryMode = static_cast<SearchBoundaryModeType>(IntSearchBoundaryMode);
    switch(IntSearchBoundaryMode)
    {
        case 0:
            {
                this->SearchBoundaryMode = SEARCH_OUTER_BOUNDARY;
                break;
            }

        case 1:
            {
                this->SearchBoundaryMode = SEARCH_OUTER_AND_INNER_BOUNDARIES;
                std::cout << "Using Outer and Inner method." << std::endl;
                break;
            }

        default:
            {
                vtkErrorMacro("Search mode is not supported.");
            }
    }
}

// Search Outer Boundary
void ExtractBoundary::SetSearchBoundaryModeToOuterBoundary()
{
    this->SearchBoundaryMode = SEARCH_OUTER_BOUNDARY;
}

// Search Outer and Inner Boundary
void ExtractBoundary::SetSearchBoundaryModeToOuterAndInnerBoundaries()
{
    this->SearchBoundaryMode = SEARCH_OUTER_AND_INNER_BOUNDARIES;
}

// ===========================
// Fill Input Port Information
// ===========================

int ExtractBoundary::FillInputPortInformation(int port, vtkInformation *info)
{
    if(port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),"vtkDataSet");
        return 1;
    }
    return 0;
}

// ============================
// Fill Output Port Information
// ============================

int ExtractBoundary::FillOutputPortInformation(int port, vtkInformation *info)
{
    if(port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkPolyData");
        return 1;
    }
    return 0;
}

// ============
// Request Data
// ============

int ExtractBoundary::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input
    vtkInformation *inputInfo = inputVector[0]->GetInformationObject(0);
    vtkDataSet *input = vtkDataSet::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Output
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);
    vtkPolyData *output = vtkPolyData::SafeDownCast(outputInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Find Whole Boundary Points
    if(this->SearchBoundaryMode == SEARCH_OUTER_BOUNDARY)
    {
        this->FindWholeOuterBoundary(input,output);
    }
    else if(this->SearchBoundaryMode == SEARCH_OUTER_AND_INNER_BOUNDARIES)
    {
        this->FindWholeOuterAndInnerBoundaries(input,output);
    }
    
    // Extract Velocity Boundary
    if(this->ExtractVelocityBoundary == true)
    {
        this->FindVelocityBoundary(input,output);
    }

    // Caclulate Normals
    if(this->CalculateNormals == true)
    {
        this->FindNormals(input,output);
    }

    return 1;
}

// ===================
// Find A Corner Point
// ===================

vtkIdType ExtractBoundary::FindACornerPoint(vtkDataSet *input)
{
    // Direction of search fo corner point
    unsigned int CheckDimension = 0; // 0=x, 1=y
    if(CheckDimension >1)
    {
        std::cerr << "Check min dimension should be only in x or y direction";
    }

    // Initial guess for corner point
    unsigned int MinXPointId = 0;
    double MinXPoint[3];
    input->GetPoint(MinXPointId,MinXPoint);
    double MinXValue = MinXPoint[CheckDimension];

    // Loop over all points for binary comparison
    double *point;
    unsigned int NumberOfPoints = input->GetNumberOfPoints();
    for(unsigned int index=1; index<NumberOfPoints; index++)
    {
        point = input->GetPoint(index);
        if(point[CheckDimension] < MinXValue)
        {
            // update with more minimum point
            MinXValue = point[CheckDimension];
            MinXPointId = index;
        }
    }

    // MinPointId has the least MinXValue
    return MinXPointId;
}

// ============================
// Find Adjacent Boundary Point
// ============================

vtkIdType ExtractBoundary::FindAdjacentBoundaryPoint(
        vtkIdType CurrentPointId,
        vtkIdType PreviousAdjacentPointId,
        vtkDataSet *input)
{
    // Cell neighbours of current point
    vtkSmartPointer<vtkIdList> CellListIds = vtkSmartPointer<vtkIdList>::New();
    input->GetPointCells(CurrentPointId,CellListIds);
   
    // Points neighbor of current point
    vtkSmartPointer<vtkIdList> Connectivities = vtkSmartPointer<vtkIdList>::New();
    this->GetPointConnectivity(CurrentPointId,input,Connectivities);

    for(vtkIdType ConnectivityIndex=0; ConnectivityIndex<Connectivities->GetNumberOfIds(); ConnectivityIndex++)
    {
        unsigned int CommonCells = 0;
        vtkSmartPointer<vtkIdList> NeighborCellListIds = vtkSmartPointer<vtkIdList>::New();
        input->GetPointCells(Connectivities->GetId(ConnectivityIndex),NeighborCellListIds);
        for(vtkIdType NeighborCellIndex=0; NeighborCellIndex<NeighborCellListIds->GetNumberOfIds(); NeighborCellIndex++)
        {
            for(vtkIdType CellListIndex = 0; CellListIndex<CellListIds->GetNumberOfIds(); CellListIndex++)
            {
                vtkIdType CellId = CellListIds->GetId(CellListIndex);
                vtkIdType NeighborCellId = NeighborCellListIds->GetId(NeighborCellIndex);
                if(CellId == NeighborCellId)
                {
                    CommonCells++;
                    break;
                }
            }
            if(CommonCells > 1)
            {
                break;
            }
        }
        if(CommonCells == 1)
        {
            vtkIdType NextAdjacentPointId = Connectivities->GetId(ConnectivityIndex);
            if(NextAdjacentPointId != PreviousAdjacentPointId)
            {
                return NextAdjacentPointId;
            }
        }
    }

    // Else, return an error with negative Id
    return -1;
}

// ======================
// Get Point Connectivity
// ======================

void ExtractBoundary::GetPointConnectivity(
        vtkIdType CurrentPointId,
        vtkDataSet *input,
        vtkIdList *Connectivities)
{
    // Neighbor cells of current point
    vtkSmartPointer<vtkIdList> CellListIds = vtkSmartPointer<vtkIdList>::New();
    input->GetPointCells(CurrentPointId,CellListIds);
    vtkIdType NumberOfNeighborCells = CellListIds->GetNumberOfIds();

    // Initialize Number of Connectivities
    vtkIdType NumberOfConnectivities = 0;

    // Loop Over neighbor cells
    for(vtkIdType CellListIndex=0; CellListIndex<NumberOfNeighborCells; CellListIndex++)
    {
        // Select a cell from neighbor cells
        vtkIdType CellId = CellListIds->GetId(CellListIndex);

        // Get points of selected cell
        vtkSmartPointer<vtkIdList> CellPointsList = vtkSmartPointer<vtkIdList>::New();
        input->GetCellPoints(CellId,CellPointsList);
        vtkIdType NumberOfCellPoints = CellPointsList->GetNumberOfIds();

        // Loop over points in selected neighbor cell
        for(vtkIdType PointIndex=0; PointIndex<NumberOfCellPoints; PointIndex++)
        {
            // Select a point from points of selected neighbor cell
            vtkIdType CheckPointId = CellPointsList->GetId(PointIndex);

            // Check if selected point is not the original inquiry point
            if(CheckPointId != CurrentPointId)
            {
                bool RepeatedBefore = false;
                if(NumberOfConnectivities != 0)
                {
                    // Loop over previous connected points
                    for(vtkIdType ConnectivityIndex=0;
                            ConnectivityIndex<NumberOfConnectivities;
                            ConnectivityIndex++)
                    {
                        vtkIdType ConnectedPointId = Connectivities->GetId(ConnectivityIndex);
                        if(CheckPointId == ConnectedPointId)
                        {
                            RepeatedBefore = true;
                        }
                    }
                }
                if(RepeatedBefore == false)
                {
                    Connectivities->InsertNextId(CheckPointId);
                    NumberOfConnectivities++;
                }
            }
        }
    }
}

// =========================
// Find Whole Outer Boundary
// =========================

void ExtractBoundary::FindWholeOuterBoundary(vtkDataSet *input,vtkPolyData *output)
{
    // Initialize Boundary Points, Vertices
    vtkSmartPointer<vtkPoints> BoundaryPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> BoundaryVertices = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolygon> BoundaryPolygon = vtkSmartPointer<vtkPolygon>::New();

    // Create id list of Boundary Points
    if(this->ExtractVelocityBoundary == true)
    {
        this->BoundaryPointIdsInInputGrid = vtkSmartPointer<vtkIdList>::New();
    }

    // Find a corner point
    vtkIdType CornerPointId = this->FindACornerPoint(input);

    // Add corner point to BoundaryPoints
    double *CornerPoint = input->GetPoint(CornerPointId);
    vtkIdType BoundaryPointId[1];
    BoundaryPointId[0] = BoundaryPoints->InsertNextPoint(CornerPoint[0],CornerPoint[1],0);
    BoundaryVertices->InsertNextCell(1,BoundaryPointId);
    BoundaryPolygon->GetPointIds()->InsertNextId(BoundaryPointId[0]);

    // Update BoundaryPointIdsInInputGrid
    if(this->ExtractVelocityBoundary == true)
    {
        this->BoundaryPointIdsInInputGrid->InsertNextId(CornerPointId);
    }

    // Initialize Previous, Current and Next Adjacent points
    vtkIdType PreviousAdjacentPointId = CornerPointId;
    vtkIdType CurrentBoundaryPointId = CornerPointId;
    vtkIdType NextAdjacentPointId;

    // Find Next Adjacent Point
    bool BoundaryRepeated = false;
    unsigned int NumberOfProcessedPoints = 0;
    while(!BoundaryRepeated)
    {
        if(NumberOfProcessedPoints >= static_cast<unsigned int>(input->GetNumberOfPoints()))
        {
            vtkErrorMacro("Can not find a closed boundary curve.");
        }

        NextAdjacentPointId = this->FindAdjacentBoundaryPoint(CurrentBoundaryPointId,PreviousAdjacentPointId,input);

        // Check not to repead boundary
        if(NextAdjacentPointId != CornerPointId)
        {
            double *BoundaryPoint = input->GetPoint(NextAdjacentPointId);
            vtkIdType BoundaryPointId[1];
            BoundaryPointId[0] = BoundaryPoints->InsertNextPoint(BoundaryPoint[0],BoundaryPoint[1],0);
            BoundaryVertices->InsertNextCell(1,BoundaryPointId);
            BoundaryPolygon->GetPointIds()->InsertNextId(BoundaryPointId[0]);

            // Update BoundaryPointIdsInInputGrid
            if(this->ExtractVelocityBoundary == true)
            {
                this->BoundaryPointIdsInInputGrid->InsertNextId(NextAdjacentPointId);
            }

            // Update for next iteration
            PreviousAdjacentPointId = CurrentBoundaryPointId;
            CurrentBoundaryPointId = NextAdjacentPointId;
        }
        else
        {
            BoundaryRepeated = true;
        }
        NumberOfProcessedPoints++;
    }

    // Create List of Polygons
    vtkSmartPointer<vtkCellArray> BoundaryPolygons = vtkSmartPointer<vtkCellArray>::New();
    BoundaryPolygons->InsertNextCell(BoundaryPolygon);

    // Set boundary points,vertices to output polydata
    output->SetPoints(BoundaryPoints);
    output->SetVerts(BoundaryVertices);
    output->SetPolys(BoundaryPolygons);
}

// =====================================
// Find Whole Outer and Inner Boundaries
// =====================================

void ExtractBoundary::FindWholeOuterAndInnerBoundaries(vtkDataSet *input, vtkPolyData *output)
{
    // Initialize Boundary Points, Vertices
    vtkSmartPointer<vtkPoints> BoundaryPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> BoundaryVertices = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolygon> BoundaryPolygon = NULL;
    vtkSmartPointer<vtkCellArray> BoundaryPolygons = vtkSmartPointer<vtkCellArray>::New();

    // Creat Id list of boundary points
    if(this->ExtractVelocityBoundary == true)
    {
        this->BoundaryPointIdsInInputGrid = vtkSmartPointer<vtkIdList>::New();
    }

    // Create PointProcessed array
    vtkIdType NumberOfPoints = input->GetNumberOfPoints();
    bool *PointsProcessed = new bool[NumberOfPoints];
    for(unsigned int i=0; i < static_cast<unsigned int>(NumberOfPoints); i++)
    {
        PointsProcessed[i] = false;
    }

    // Search in All Points for a boundary Point
    for(vtkIdType CurrentPointId=0; CurrentPointId<NumberOfPoints; CurrentPointId++)
    {
        if(PointsProcessed[CurrentPointId] == true)
        {
            continue;
        }

        vtkIdType AdjacentPointId[1];
        bool IsBoundaryPointStatus = this->IsBoundaryPoint(CurrentPointId,AdjacentPointId,input);
        PointsProcessed[CurrentPointId] = true;
        // std::cout << IsBoundaryPointStatus << std::endl;

        if(IsBoundaryPointStatus == true)
        {
            // Define two known points on boundary
            vtkIdType CurrentBoundaryPointId = CurrentPointId;
            vtkIdType PreviousAdjacentPointId = AdjacentPointId[0];
            vtkIdType FirstPointInCurveId = AdjacentPointId[0];
            vtkIdType NextAdjacentPointId;

            // Update BoundaryPointIdsInInputGrid
            if(this->ExtractVelocityBoundary == true)
            {
                this->BoundaryPointIdsInInputGrid->InsertNextId(PreviousAdjacentPointId);
                this->BoundaryPointIdsInInputGrid->InsertNextId(CurrentBoundaryPointId);
            }

            PointsProcessed[PreviousAdjacentPointId] = true;

            // Add Previous point to BoundaryPoints and Vertices
            vtkIdType PreviousAdjacentPointBoundaryId[1];
            double *PreviousAdjacentPoint = input->GetPoint(PreviousAdjacentPointId);
            PreviousAdjacentPointBoundaryId[0] = BoundaryPoints->InsertNextPoint(PreviousAdjacentPoint[0],PreviousAdjacentPoint[1],0);
            BoundaryVertices->InsertNextCell(1,PreviousAdjacentPointBoundaryId);

            // Add Current Point to BoundaryPoint and Vertices
            vtkIdType CurrentBoundaryPointBoundaryId[1];
            double *CurrentBoundaryPoint = input->GetPoint(CurrentBoundaryPointId);
            CurrentBoundaryPointBoundaryId[0] = BoundaryPoints->InsertNextPoint(CurrentBoundaryPoint[0],CurrentBoundaryPoint[1],0);
            BoundaryVertices->InsertNextCell(1,CurrentBoundaryPointBoundaryId);

            // Initialize new Polygon for each curve
            BoundaryPolygon = vtkSmartPointer<vtkPolygon>::New();

            // Add Previous Piont to Polygon
            BoundaryPolygon->GetPointIds()->InsertNextId(PreviousAdjacentPointBoundaryId[0]);

            // Add Current Point to Polygon
            BoundaryPolygon->GetPointIds()->InsertNextId(CurrentBoundaryPointBoundaryId[0]);

            // Find whole boundary curve
            bool BoundaryRepeated = false;
            unsigned int NumberOfProcessedPoints = 0;
            while(!BoundaryRepeated)
            {
                if(NumberOfProcessedPoints >= static_cast<unsigned int>(NumberOfPoints))
                {
                    vtkErrorMacro("Can not find a closed boundary curve.");
                }

                NextAdjacentPointId = this->FindAdjacentBoundaryPoint(CurrentBoundaryPointId,PreviousAdjacentPointId,input);

                // Check not to repeat boundary
                if(NextAdjacentPointId != FirstPointInCurveId)
                {
                    PointsProcessed[NextAdjacentPointId] = true;

                    double *BoundaryPoint = input->GetPoint(NextAdjacentPointId);
                    vtkIdType BoundaryPointId[1];
                    BoundaryPointId[0] = BoundaryPoints->InsertNextPoint(BoundaryPoint[0],BoundaryPoint[1],0);
                    BoundaryVertices->InsertNextCell(1,BoundaryPointId);
                    BoundaryPolygon->GetPointIds()->InsertNextId(BoundaryPointId[0]);

                    // Update BoundaryPointIdsInInputGrid
                    if(this->ExtractVelocityBoundary == true)
                    {
                        this->BoundaryPointIdsInInputGrid->InsertNextId(NextAdjacentPointId);
                    }

                    // Update for next  iteration
                    PreviousAdjacentPointId = CurrentBoundaryPointId;
                    CurrentBoundaryPointId = NextAdjacentPointId;
                }
                else
                {
                    BoundaryRepeated = true;
                }
                NumberOfProcessedPoints++;
            }

            // Add Polygon to list of Polygons
            BoundaryPolygons->InsertNextCell(BoundaryPolygon);
        }
    }

    // Set Boundary Points and Vertices to output polydata
    output->SetPoints(BoundaryPoints);
    output->SetVerts(BoundaryVertices);
    output->SetPolys(BoundaryPolygons);

    // Freed Memory
    if(PointsProcessed != NULL)
    {
        delete [] PointsProcessed;
        PointsProcessed = NULL;
    }
}

// =================
// Is Boundary Point
// =================

bool ExtractBoundary::IsBoundaryPoint(
        vtkIdType CurrentPointId,
        vtkIdType AdjacentPointId[1],
        vtkDataSet *input)
{
    // Points neighbor of current point
    vtkSmartPointer<vtkIdList> Connectivities = vtkSmartPointer<vtkIdList>::New();
    this->GetPointConnectivity(CurrentPointId,input,Connectivities);

    for(vtkIdType ConnectivityIndex=0; ConnectivityIndex < Connectivities->GetNumberOfIds(); ConnectivityIndex++)
    {
        vtkIdType AdjacentTestPointId = Connectivities->GetId(ConnectivityIndex);
        int SharedCellId = this->CheckSharedCells(CurrentPointId,AdjacentTestPointId,input);
        if(SharedCellId >= 0)
        {
            // unique common cell if found
            AdjacentPointId[0] = AdjacentTestPointId;
            return true;
        }
    }
    
    // If no unique common cell is found
    return false;
}

// ==================
// Check Shared Cells
// ==================

// Description:
// Finds shared cells between two points.
// Return value >  0 : if shared cell is unique then it returns that cell id
// Return value = -1 : if shared cells are not unique
// Return value = -2 : if no shared cell is found

int ExtractBoundary::CheckSharedCells(
        vtkIdType PointId1,
        vtkIdType PointId2,
        vtkDataSet *input)
{
    unsigned int SharedCells = 0;

    // Cells of Point 1
    vtkSmartPointer<vtkIdList> CellIdsOfPoint1 = vtkSmartPointer<vtkIdList>::New();
    input->GetPointCells(PointId1,CellIdsOfPoint1);
    vtkIdType NumberOfCellsOfPoint1 = CellIdsOfPoint1->GetNumberOfIds();

    // Cells of Point 2
    vtkSmartPointer<vtkIdList> CellIdsOfPoint2 = vtkSmartPointer<vtkIdList>::New();
    input->GetPointCells(PointId2,CellIdsOfPoint2);
    vtkIdType NumberOfCellsOfPoint2 = CellIdsOfPoint2->GetNumberOfIds();

    // Check Shared Cells
    int SharedCellId = -2;
    for(unsigned int CellIdIndexOfPoint1=0; CellIdIndexOfPoint1 < static_cast<unsigned int>(NumberOfCellsOfPoint1); CellIdIndexOfPoint1++)
    {
        vtkIdType CellIdOfPoint1 = CellIdsOfPoint1->GetId(CellIdIndexOfPoint1);
        for(unsigned int CellIdIndexOfPoint2=0; CellIdIndexOfPoint2 < static_cast<unsigned int>(NumberOfCellsOfPoint2); CellIdIndexOfPoint2++)
        {
            vtkIdType CellIdOfPoint2 = CellIdsOfPoint2->GetId(CellIdIndexOfPoint2);

            // Comparision of Ids
            if(CellIdOfPoint1 == CellIdOfPoint2)
            {
                SharedCellId = CellIdOfPoint1;
                SharedCells++;
            }

            // Check if Shared Cells are not unique
            if(SharedCells > 1)
            {
                return -1;
            }
        }
    }

    // else, sharedCells are unique
    return SharedCellId;
}

// ======================
// Find Velocity Boundary
// ======================

void ExtractBoundary::FindVelocityBoundary(vtkDataSet *input, vtkPolyData *output)
{
    // Velocity Boundary Array
    vtkSmartPointer<vtkUnsignedIntArray> VelocityBoundary = vtkSmartPointer<vtkUnsignedIntArray>::New();
    VelocityBoundary->SetNumberOfComponents(1);
    VelocityBoundary->SetNumberOfTuples(output->GetNumberOfPoints());
    VelocityBoundary->SetName("VelocityBoundary");
    unsigned int *pVelocityBoundary = VelocityBoundary->GetPointer(0);

    // Check BoundaryPointIdsInInputGrid
    if(this->BoundaryPointIdsInInputGrid == NULL)
    {
        vtkErrorMacro("BoundaryPointIdsInInputGrid is NULL");
    }
    else if(this->BoundaryPointIdsInInputGrid->GetNumberOfIds() != output->GetNumberOfPoints())
    {
        vtkErrorMacro("BoundaryPointIdsInInputGrid does not match output point numbers.");
    }

    // Abstract Velocity Array
    vtkSmartPointer<vtkDataArray> VelocityArray = input->GetPointData()->GetVectors();
    if(VelocityArray == NULL)
    {
        vtkErrorMacro("input data does not have vector array.");
    }
    else if(VelocityArray->GetNumberOfTuples() == 0)
    {
        vtkErrorMacro("VelocityArray is empty.");
    }

    // Cast Abstract array to double type
    // double *pVelocities = static_cast<double*>(VelocityArray->GetVoidPointer(0));

    // vtkSmartPointer<vtkFloatArray> fVelocityArray = vtkFloatArray::SafeDownCast(VelocityArray);
    // double *pVelocities = Velocities->GetPointer(0);
    // float *pVelocities = fVelocityArray->GetPointer(0);

    // vtkSmartPointer<vtkDataArray> Velocities = input->GetPointData()->GetVectors();
    // float *pVelocities = static_cast<float*>(Velocities->GetVoidPointer(0));

    // Look for near zero velocity
    for(vtkIdType i=0; i<output->GetNumberOfPoints(); i++)
    {
        vtkIdType IdInInputGrid = this->BoundaryPointIdsInInputGrid->GetId(i);
        double *pVelocities = VelocityArray->GetTuple3(IdInInputGrid);
        // PROBLEM: how data arranged in array?
        // unsigned int VelocityMagnitude =
        //     fabs(pVelocities[3*IdInInputGrid]) +
        //     fabs(pVelocities[3*IdInInputGrid+1]) +
            fabs(pVelocities[3*IdInInputGrid+2]);
        unsigned int VelocityMagnitude =
            fabs(pVelocities[0]) +
            fabs(pVelocities[1]) +
            fabs(pVelocities[2]);

        if(VelocityMagnitude < EPSILON)
        {
            pVelocityBoundary[i] = 1;
        }
        else
        {
            pVelocityBoundary[i] = 0;
        }
    }

    // Write to output
    output->GetPointData()->AddArray(VelocityBoundary);
}

// ============
// Find Normals
// ============

void ExtractBoundary::FindNormals(vtkDataSet *input,vtkPolyData *output)
{
    // Polygons and Polygon
    vtkSmartPointer<vtkCellArray> Polygons = output->GetPolys();
    unsigned int NumberOfPolygons = Polygons->GetNumberOfCells();
    vtkSmartPointer<vtkIdList> Polygon[NumberOfPolygons]; // MAC problem, use line below instead
    // vtkIdList *Polygon[NumberOfPolygons];

    // initialize Normals Array
    vtkSmartPointer<vtkDoubleArray> Normals = vtkSmartPointer<vtkDoubleArray>::New();
    Normals->SetNumberOfComponents(3);
    Normals->SetName("NormalVectors");

    // for each Polygon
    for(unsigned int i=0; i<NumberOfPolygons; i++)
    {
        Polygon[i] = vtkSmartPointer<vtkIdList>::New(); // Mac problem, use line below instead, them do deletion
        // Polygon[i] = vtkIdList::New();
        int PolygonStatus = Polygons->GetNextCell(Polygon[i]);
        if(PolygonStatus == 0)
        {
            vtkErrorMacro("Polygon is empty");
        }

        // for each point in Polygon
        vtkIdType NumberOfPoints = Polygon[i]->GetNumberOfIds();
        for(vtkIdType j=0; j<NumberOfPoints; j++)
        {
            // Point indices in Loop
            vtkIdType PreviousPointIdInLoop = ((j-1)>=0 ? (j-1) : (NumberOfPoints-j-1));
            vtkIdType CurrentPointIdInLoop = j;
            vtkIdType NextPointIdInLoop = ((j+1)<NumberOfPoints ? (j+1) : (j+1-NumberOfPoints));

            // Point indices in Polygon
            vtkIdType PreviousPointIdInPolygon = Polygon[i]->GetId(PreviousPointIdInLoop);
            vtkIdType CurrentPointIdInPolygon = Polygon[i]->GetId(CurrentPointIdInLoop);
            vtkIdType NextPointIdInPolygon = Polygon[i]->GetId(NextPointIdInLoop);

            // Point Indices in Input Grid
            vtkIdType PreviousPointIdInInputGrid = this->BoundaryPointIdsInInputGrid->GetId(PreviousPointIdInPolygon);
            vtkIdType CurrentPointIdInInputGrid = this->BoundaryPointIdsInInputGrid->GetId(CurrentPointIdInPolygon);
            vtkIdType NextPointIdInInputGrid = this->BoundaryPointIdsInInputGrid->GetId(NextPointIdInPolygon);

            // Get Boundary Cell Ids
            vtkIdType PreviousBoundaryCellId = this->CheckSharedCells(CurrentPointIdInInputGrid,PreviousPointIdInInputGrid,input);
            if(PreviousBoundaryCellId < 0)
            {
                vtkErrorMacro("These points are not on boundary or not adjacent.");
            }
            vtkIdType NextBoundaryCellId = this->CheckSharedCells(NextPointIdInInputGrid,CurrentPointIdInInputGrid,input);
            if(NextBoundaryCellId < 0)
            {
                vtkErrorMacro("These points are not on boundary or not adjacent.");
            }

            // Get Boundary Cells
            vtkSmartPointer<vtkGenericCell> PreviousBoundaryCell = vtkSmartPointer<vtkGenericCell>::New();
            vtkSmartPointer<vtkGenericCell> NextBoundaryCell = vtkSmartPointer<vtkGenericCell>::New();

            input->GetCell(PreviousBoundaryCellId,PreviousBoundaryCell);
            input->GetCell(NextBoundaryCellId,NextBoundaryCell);

            PreviousBoundaryCell->SetCellTypeToTriangle();
            NextBoundaryCell->SetCellTypeToTriangle();

            if(PreviousBoundaryCell == NULL)
            {
                vtkErrorMacro("Previous cell is NULL.");
            }
            if(NextBoundaryCell == NULL)
            {
                vtkErrorMacro("Next cell is NULL.");
            }

            // Get Boundary Cell Centers in Parametric Coordinates
            double PreviousCellCenterPCoords[3];
            double NextCellCenterPCoords[3];
            PreviousBoundaryCell->GetParametricCenter(PreviousCellCenterPCoords);
            NextBoundaryCell->GetParametricCenter(NextCellCenterPCoords);

            // Get Boundary Cell Centers in Global Coordinates
            int SubId = 0;
            double PreviousCellCenterGCoords[3];
            double NextCellCenterGCoords[3];
            double PreviousWeights[3];
            double NextWeights[3];
            PreviousBoundaryCell->EvaluateLocation(SubId,PreviousCellCenterPCoords,PreviousCellCenterGCoords,PreviousWeights);
            NextBoundaryCell->EvaluateLocation(SubId,NextCellCenterPCoords,NextCellCenterGCoords,NextWeights);

            // Get cell point global coordinates
            double PreviousPointGCoords[3];
            double CurrentPointGCoords[3];
            double NextPointGCoords[3];
            input->GetPoint(PreviousPointIdInInputGrid,PreviousPointGCoords);
            input->GetPoint(CurrentPointIdInInputGrid,CurrentPointGCoords);
            input->GetPoint(NextPointIdInInputGrid,NextPointGCoords);

            // Get Altitude of each cell
            double PreviousCellAltitude[3];
            double NextCellAltitude[3];
            this->GetTriangleAltitude(CurrentPointGCoords,PreviousPointGCoords,PreviousCellCenterGCoords,PreviousCellAltitude);
            this->GetTriangleAltitude(CurrentPointGCoords,NextPointGCoords,NextCellCenterGCoords,NextCellAltitude);

            // Normalalize Altitudes
            vtkMath::Normalize(PreviousCellAltitude);
            vtkMath::Normalize(NextCellAltitude);

            // Normal to polygon
            double NormalVectorInCurrentPoint[3];
            vtkMath::Add(PreviousCellAltitude,NextCellAltitude,NormalVectorInCurrentPoint);

            // Normalize Normal vector
            vtkMath::Normalize(NormalVectorInCurrentPoint);

            // Outward normals
            if(this->NormalsOutward == true)
            {
                vtkMath::MultiplyScalar(NormalVectorInCurrentPoint,-1);
            }

            // put result in Normals Array
            Normals->InsertNextTuple3(
                    NormalVectorInCurrentPoint[0],
                    NormalVectorInCurrentPoint[1],
                    NormalVectorInCurrentPoint[2]);
        }
    }

    // Append Normals to output
    output->GetPointData()->SetNormals(Normals);
}

// =====================
// Get Triangle Altitude
// =====================

void ExtractBoundary::GetTriangleAltitude(
        double BasePoint1[3],
        double BasePoint2[3],
        double AltitudePoint[3],
        double AltitudeVector[3])
{
    // Base of triangle
    double BaseVector[3];
    vtkMath::Subtract(BasePoint2,BasePoint1,BaseVector);

    // one side of triangle
    double SideVector[3];
    vtkMath::Subtract(AltitudePoint,BasePoint1,SideVector);
    AltitudePoint[2] = 0;

    // projection of side on base
    double SideOnBaseProjectionVector[3];
    vtkMath::ProjectVector(SideVector,BaseVector,SideOnBaseProjectionVector);
    SideOnBaseProjectionVector[2] = 0;
    
    // Altitude
    vtkMath::Subtract(SideVector,SideOnBaseProjectionVector,AltitudeVector);
    AltitudeVector[2] = 0;

    // Considering 2D
    AltitudeVector[2] = 0.0;
}

// ==================================
// Get Search Boundary Mode As String
// ==================================

const char * ExtractBoundary::GetSearchBoundaryModeAsString()
{
    switch(this->SearchBoundaryMode)
    {
        // Search Outer Boundary
        case SEARCH_OUTER_BOUNDARY:
        {
            return "SEARCH_OUTER_BOUNDARY";
        }

        // Search Outer and Inner Boundaries
        case SEARCH_OUTER_AND_INNER_BOUNDARIES:
        {
            return "SEARCH_OUTER_AND_INNER_BOUNDARIES";
        }

        // Not supported
        default:
        {
            vtkErrorMacro("No such mode exists.");
            return NULL;
        }
    }
}
