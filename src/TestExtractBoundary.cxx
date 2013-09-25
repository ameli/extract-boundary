/*
 * =====================================================================================
 *
 *       Filename:  test.cxx
 *
 *    Description:  Test for Extract Boundary
 *
 *        Version:  1.0
 *        Created:  12/09/2012 10:06:32 PM
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

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include "ExtractBoundary.h"
#include <vtkPolyDataWriter.h>

// ====
// Main
// ====

int main(int argc, char *argv[])
{
    // Reader
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    // reader->SetFileName("/home/sia/code/projvtk/Filter/ExtractBoundary/data/InputData1.vtk");
    reader->SetFileName("/home/sia/code/projvtk/Filter/ExtractBoundary/data/InputData2.vtk");
    reader->Update();

    // Extract Boundary
    vtkSmartPointer<ExtractBoundary> extractboundary = vtkSmartPointer<ExtractBoundary>::New();
    extractboundary->SetInputConnection(reader->GetOutputPort());
    extractboundary->SetSearchBoundaryModeToOuterAndInnerBoundaries();
    extractboundary->Update();

    // Writer
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputConnection(extractboundary->GetOutputPort());
    writer->SetFileName("/home/sia/code/projvtk/Filter/ExtractBoundary/bin/boundary.vtk");
    writer->Write();

    return EXIT_SUCCESS;
}
