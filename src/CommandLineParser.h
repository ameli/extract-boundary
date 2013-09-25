/*
 * =====================================================================================
 *
 *       Filename:  CommandLineParser.h
 *
 *    Description:  Parser class for command line inputs
 *
 *        Version:  1.0
 *        Created:  01/15/2013 12:21:17 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __CommandLineParser_h
#define __CommandLineParser_h

// ======
// Macros
// ======

// Debug
#ifndef HERE
#define HERE std::cout << "DEBUG: " << __FILE__ << " at " << __LINE__ << std::endl;
#endif

// Version Major
#ifndef VERSION_MAJOR
#define VERSION_MAJOR 0
#endif

// Version Minor
#ifndef VERSION_MINOR
#define VERSION_MINOR 1
#endif

// Max Filename Size
#ifndef MAX_FILENAME_SIZE
#define MAX_FILENAME_SIZE
#endif

// =======
// Headers
// =======

#include <map>
#include <string>

// =====
// Types
// =====

enum OptionsListEnumType
{
    OPTION_HELP = 0,     // mapped to; -h
    OPTION_ABOUT,        // mapped to: -a
    OPTION_LICENSE,      // mapped to: -l
    OPTION_VERSION,      // mapped to: -v
    NUMBER_OF_OPTIONS
};

// =========================
// Command Line Parser Class
// =========================

class CommandLineParser
{
    public:
        CommandLineParser(int argc, char *argv[]);
        ~CommandLineParser();

        // Accessors, Mutators
        const char * GetInputFileName() const;
        const char * GetOutputFileName() const;

    protected:
        // member methods
        void CheckArguments();
        void CheckOptions();
        void ProcessOptions();
        void FindInputFileArgument();
        void FindOutputFileArgument();
        void PrintHelp();
        void PrintOptions();
        void PrintLicense();
        void PrintAbout();
        void PrintVersion();
        static void DecomposeFileName(
                char * FileName,
                char * FileBaseName,
                char * FileExtansion);

        // Member data
        int argc;
        char **argv;

        // Internal Member Data
        char * InputFileName;
        char * InputFileBaseName;
        char * InputFileExtension;
        typedef std::map<std::string,OptionsListEnumType> OptionsListMapType;
        static OptionsListMapType OptionsList;
        static OptionsListMapType CreateOptionsList();
};

#endif
