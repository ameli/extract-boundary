/*
 * =====================================================================================
 *
 *       Filename:  CommandLineParser.cxx
 *
 *    Description:  Parser class for command line inputs
 *
 *        Version:  1.0
 *        Created:  01/15/2013 12:46:28 AM
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
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <cstdlib>

// ======
// Macros
// ======

#ifndef MAX_CHAR_LENGTH
#define MAX_CHAR_LENGTH 512
#endif

// ===========
// Constructor
// ===========

CommandLineParser::CommandLineParser(int ArgC, char *ArgV[])
{
    // Member data
    this->argc = ArgC;
    this->argv = ArgV;

    this->InputFileName = NULL;
    this->InputFileBaseName = NULL;
    this->InputFileExtension = NULL;

    // Check
    this->CheckArguments();
    this->CheckOptions();
    this->ProcessOptions();
    this->FindInputFileArgument();
    this->FindOutputFileArgument();
}

// ==========
// Destructor
// ==========

CommandLineParser::~CommandLineParser()
{
    // Input file base name
    if(this->InputFileBaseName != NULL)
    {
        delete [] this->InputFileBaseName;
        this->InputFileBaseName = NULL;
    }

    // Input file extension
    if(this->InputFileExtension != NULL)
    {
        delete [] this->InputFileExtension;
        this->InputFileExtension = NULL;
    }
}

// ===============
// Check Arguments
// ===============

void CommandLineParser::CheckArguments()
{
    if(argc < 2)
    {
        this->PrintHelp();
        exit(EXIT_FAILURE);
    }
}

// =============
// Check Options
// =============

void CommandLineParser::CheckOptions()
{
    bool UnknownOptionEntered = false;

    // Check each input argument
    for(unsigned int i=1; i< static_cast<unsigned int>(this->argc); i++)
    {
        // Check if input is option
        if(strncmp(argv[i],"-",1))
        {
            // not an option
            continue;
        }

        // Find map of option to OptionListType
        std::string InputOptionAsString = std::string(argv[i]);
        std::map<std::string,OptionsListEnumType>::iterator OptionListIterator;
        OptionListIterator = CommandLineParser::OptionsList.find(InputOptionAsString);
        OptionsListEnumType InputOptionAsIndex = OptionListIterator->second;

        // Unknown option
        if(InputOptionAsIndex < 0 || InputOptionAsIndex > NUMBER_OF_OPTIONS)
        {
            std::cerr << "Unknown option: " << argv[i] << std::endl;
            UnknownOptionEntered = true;
        }
    }

    // Unknown option entered
    if(UnknownOptionEntered == true)
    {
        std::cerr << std::endl;
        this->PrintOptions();
        exit(EXIT_FAILURE);
    }
}

// ===============
// Process Options
// ===============

void CommandLineParser::ProcessOptions()
{
    bool UnknownOptionEntered = false;

    // Check each input argument
    for(unsigned int i=1; i<static_cast<unsigned int>(this->argc); i++)
    {
        // Check if input is option
        if(strncmp(argv[i],"-",1))
        {
            // not an option
            continue;
        }

        // Find map of option to OptionListType
        std::string InputOptionAsString = std::string(argv[i]);
        std::map<std::string,OptionsListEnumType>::iterator OptionListIterator;
        OptionListIterator = CommandLineParser::OptionsList.find(InputOptionAsString);
        OptionsListEnumType InputOptionAsIndex = OptionListIterator->second;

        switch(InputOptionAsIndex)
        {
            // Help -h
            case OPTION_HELP:
            {
                this->PrintHelp();
                this->PrintOptions();
                exit(EXIT_SUCCESS);
            }

            // License -l
            case OPTION_LICENSE:
            {
                this->PrintLicense();
                exit(EXIT_SUCCESS);
            }

            // About -a
            case OPTION_ABOUT:
            {
                this->PrintAbout();
                exit(EXIT_SUCCESS);
            }

            // Version -v
            case OPTION_VERSION:
            {
                this->PrintVersion();
                exit(EXIT_SUCCESS);
            }

            // Unknown option
            default:
            {
                std::cerr << "Unknown option: " << argv[i] << std::endl;
                UnknownOptionEntered = true;
                break;
            }
        }
    }

    // Unknown option entered
    if(UnknownOptionEntered == true)
    {
        std::cerr << std::endl;
        this->PrintOptions();
        exit(EXIT_FAILURE);
    }
}

// ========================
// Find Input File Argument
// ========================

void CommandLineParser::FindInputFileArgument()
{
    if(!strncmp(this->argv[1],"-",1))
    {
        // first arg is an option
        std::cout << std::endl;
        std::cerr << "Input file name is not specified in first argument." << std::endl;
        this->PrintHelp();
        exit(EXIT_FAILURE);
    }

    // first arg appears to be input file
    this->InputFileName = this->argv[1];

    // Allocate mempry for member data
    this->InputFileBaseName = new char[MAX_CHAR_LENGTH];
    this->InputFileExtension = new char[MAX_CHAR_LENGTH];

    // Extract file extension
    CommandLineParser::DecomposeFileName(
            this->InputFileName,
            this->InputFileBaseName,
            this->InputFileExtension);
}

// =========================
// Find Output File Argument
// =========================

void CommandLineParser::FindOutputFileArgument()
{
    // TODO
}

// ==========
// Print Help
// ==========

void CommandLineParser::PrintHelp()
{
    std::cout << "Usages: " << std:: endl;
    std::cout << argv[0] << " </FullPath/InputFileName> </FullPath/OutputFileName>"  << std::endl;
    std::cout << argv[0] << " </FullPath/InputFileName> <OutputFileName>"  << std::endl;
    std::cout << argv[0] << " <InputFileName> <OutputFileName>"  << std::endl;
    std::cout << argv[0] << " </FullPath/InputFileName>" << std::endl;
    std::cout << argv[0] << " <InputFileName>" << std::endl;
    std::cout << std::endl;
}

// =============
// Print Options
// =============

void CommandLineParser::PrintOptions()
{
    std::cout << "Options: " << std::endl;
    std::cout << "-a, --about   :  show information about code and exit" << std::endl;
    std::cout << "-h, --help    :  show help and exit" << std::endl;
    std::cout << "-l, --license :  show license and exit" << std::endl;
    std::cout << "-v, --version :  show version and exit" << std::endl;
}

// ===========
// Print About
// ===========

void CommandLineParser::PrintAbout()
{
    std::cout << "Author: Siavash Ameli <sameli at berkeley dot edu>" << std::endl;
    std::cout << "Shadden Research Group" << std::endl;
    std::cout << "University of California, Berkeley" << std::endl;
}

// =============
// Print License
// =============

void CommandLineParser::PrintLicense()
{
    std::ifstream LicenseFile;
    LicenseFile.open("../License.txt",std::ios::in);
    if(!LicenseFile.is_open())
    {
        std::cerr << "Can not find License file." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Read file
    std::string Line = "";
    while(LicenseFile.good())
    {
        std::getline(LicenseFile,Line);
        std::cout << Line << std::endl;
    }
    
    LicenseFile.close();
}

// =============
// Print Version
// =============

void CommandLineParser::PrintVersion()
{
    std::cout << "extractboundary " << VERSION_MAJOR << "." << VERSION_MINOR << std::endl;
}

// ===================
// Create Options List
// ===================

CommandLineParser::OptionsListMapType CommandLineParser::CreateOptionsList()
{
    OptionsListMapType tempOptionsList;

    tempOptionsList["-h"]        = OPTION_HELP;
    tempOptionsList["--help"]    = OPTION_HELP;
    tempOptionsList["-a"]        = OPTION_ABOUT;
    tempOptionsList["--about"]   = OPTION_ABOUT;
    tempOptionsList["-l"]        = OPTION_LICENSE;
    tempOptionsList["--license"] = OPTION_LICENSE;
    tempOptionsList["-v"]        = OPTION_VERSION;
    tempOptionsList["--version"] = OPTION_VERSION;

    return tempOptionsList;
}

// ==================
// Static Member Data
// ==================

CommandLineParser::OptionsListMapType CommandLineParser::OptionsList(CommandLineParser::CreateOptionsList());

// ===================
// Get Input File Name
// ===================

const char * CommandLineParser::GetInputFileName() const
{
    return this->argv[1];
}

// ====================
// Get Output File Name
// ====================

const char * CommandLineParser::GetOutputFileName() const
{
    char *OutputFileName = NULL;
   
    // No output filename specified
    if(this->argc<3 || !strncmp(argv[2],"-",1))
    {
        std::string InputFileName(this->GetInputFileName());
        size_t DotPosition = InputFileName.find(".");
        std::string InputFileNameWithoutExtention = (std::string::npos == DotPosition) ? InputFileName : InputFileName.substr(0,DotPosition);
        std::stringstream OutputFileNameString;
        OutputFileNameString << InputFileNameWithoutExtention << "-boundary" << ".vtk";
        OutputFileName = new char[OutputFileNameString.str().size()+1];
        strcpy(OutputFileName,OutputFileNameString.str().c_str());
    }

    // output filename specified
    else
    {
        OutputFileName = new char[strlen(this->GetInputFileName())+1];
        strcpy(OutputFileName,this->argv[2]);
    }

    return OutputFileName;
}

// ===================
// Decompose File Name
// ===================

void CommandLineParser::DecomposeFileName(
        char * FileName,
        char * FileBaseName,
        char * FileExtension)
{
    // File name
    std::string FileNameAsString(FileName);

    // Find dot
    size_t DotPosition = FileNameAsString.find_last_of(".");

    // Base name of file
    std::string FileBaseNameAsString = (std::string::npos == DotPosition) ? FileNameAsString : FileNameAsString.substr(0,DotPosition);
    strcpy(FileBaseName,FileBaseNameAsString.c_str());

    // Extension of file
    std::string FileExtensionAsString = (std::string::npos == DotPosition) ? "" : FileNameAsString.substr(DotPosition+1);
    strcpy(FileExtension,FileExtensionAsString.c_str());
}
