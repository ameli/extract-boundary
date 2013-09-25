/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description:  Extract boundary of 2D data
 *
 *        Version:  1.0
 *        Created:  01/14/2013 02:13:08 PM
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

#include "CommandLineParser.h"
#include "ExtractBoundary.h"
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkPolyDataWriter.h>

// ====
// Main
// ====

int main(int argc, char *argv[])
{
    // Parse command line
    CommandLineParser UserInput(argc,argv);

    // Reader
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(UserInput.GetInputFileName());
    reader->Update();

    // Extract Boundary
    vtkSmartPointer<ExtractBoundary> extractboundary = vtkSmartPointer<ExtractBoundary>::New();
    extractboundary->SetSearchBoundaryModeToOuterAndInnerBoundaries();
    extractboundary->SetExtractVelocityBoundary(true);
    extractboundary->SetCalculateNormals(true);
    extractboundary->SetNormalsOutward(false);
    extractboundary->SetInputConnection(reader->GetOutputPort());
    extractboundary->Update();
 
    // Writer
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputConnection(extractboundary->GetOutputPort());
    writer->SetFileName(UserInput.GetOutputFileName());
    writer->Write();

    return EXIT_SUCCESS;
}
