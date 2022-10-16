/*
 * =====================================================================================
 *
 *       Filename:  ExtractBoundary.h
 *
 *    Description:  Extract Boundary curve of 2D data
 *
 *        Version:  1.0
 *        Created:  12/09/2012 10:08:15 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __ExtractBoundary_h
#define __ExtractBoundary_h

// ======
// Macros
// ======

// Debugging
#define HERE std::cout << "DEBUG: " << __FILE__ << " at " << __LINE__ << std::endl;

// ===========
// Enumerators
// ===========

enum SearchBoundaryModeType
{
    SEARCH_OUTER_BOUNDARY = 0,
    SEARCH_OUTER_AND_INNER_BOUNDARIES,
    SEARCH_NUMBERS
};

// ====================
// Forward Declarations
// ====================

class vtkIdList;
#include <vtkSmartPointer.h>
#include <vtkPolyDataAlgorithm.h>

// ======================
// Class Extract Boundary
// ======================

class ExtractBoundary : public vtkPolyDataAlgorithm
{
    public:
        static ExtractBoundary * New();
        // vtkTypeRevisionMacro(ExtractBoundary,vtkPolyDataAlgorithm);
        vtkTypeMacro(ExtractBoundary,vtkPolyDataAlgorithm);
        virtual void PrintSelf(ostream & os, vtkIndent indent);

        // Accessors, Mutators
        vtkGetMacro(SearchBoundaryMode,SearchBoundaryModeType);
        vtkSetMacro(SearchBoundaryMode,SearchBoundaryModeType);
        void SetSearchBoundaryMode(int IntSearchBoundaryMode);

        void SetSearchBoundaryModeToOuterBoundary();
        void SetSearchBoundaryModeToOuterAndInnerBoundaries();

        vtkGetMacro(ExtractVelocityBoundary,bool);
        vtkSetMacro(ExtractVelocityBoundary,bool);

        vtkGetMacro(CalculateNormals,bool);
        vtkSetMacro(CalculateNormals,bool);

        vtkGetMacro(NormalsOutward,bool);
        vtkSetMacro(NormalsOutward,bool);

        void SetNormalsInward() { this->NormalsOutward = false;} 
        void SetNormalsOutward() { this->NormalsOutward = true; }

    protected:
        ExtractBoundary();
        ~ExtractBoundary();

        // Pipeline Executives
        virtual int FillInputPortInformation(int port, vtkInformation *info);

        virtual int FillOutputPortInformation(int port, vtkInformation *info);

        virtual int RequestData(
                vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        // Member Methods
        vtkIdType FindACornerPoint(vtkDataSet *input);

        bool IsBoundaryPoint(
                vtkIdType CurrentPointId,
                vtkIdType AdjacentPointId[1],
                vtkDataSet *input);

        vtkIdType FindAdjacentBoundaryPoint(
                vtkIdType CurrentPointId,
                vtkIdType PreviousAdjacentPointId,
                vtkDataSet *input);

        void GetPointConnectivity(
                vtkIdType CurrentPointId,
                vtkDataSet *input,
                vtkIdList *Connectivities);

        int CheckSharedCells(
                vtkIdType PointId1,
                vtkIdType PointId2,
                vtkDataSet *input);

        void FindWholeOuterBoundary(
                vtkDataSet *input,
                vtkPolyData *output);

        void FindWholeOuterAndInnerBoundaries(
                vtkDataSet *input,
                vtkPolyData *output);

        void FindVelocityBoundary(
                vtkDataSet *input,
                vtkPolyData *output);

        void FindNormals(
                vtkDataSet *input,
                vtkPolyData *output);

        void GetTriangleAltitude(
                double BasePoint1[3],
                double BasePoint2[3],
                double AltitudePoint[3],
                double AltitudeVector[3]);

        const char * GetSearchBoundaryModeAsString();

        // Member Data
        SearchBoundaryModeType SearchBoundaryMode;
        bool ExtractVelocityBoundary;
        bool CalculateNormals;
        bool NormalsOutward;

        // Internal Member Data
        vtkSmartPointer<vtkIdList> BoundaryPointIdsInInputGrid;

    private:
        ExtractBoundary(const ExtractBoundary & rhs);   // Not implemented.
        void operator=(const ExtractBoundary & rhs);   // Not implemented.
};

#endif
