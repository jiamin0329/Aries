/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for storing processor data
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    01-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */


#include "ProcData.hpp"

#include <fstream>
#include <string>
#include <iostream>

using namespace std;

namespace ARIES
{
    ProcData::ProcData()
    {
    }

    ProcData::~ProcData()
    {
    }

    unsigned short ProcData::GetNumDim()
    {
        string text_line, Marker_Tag;
        ifstream meshFile;
        short nDim = 3;
        unsigned short iLine, nLine = 10;
        string::size_type position;
        
        //Open grid file
        meshFile.open(d_meshFileName, ios::in);
        if (meshFile.fail())
        {
            cout << "Mesh File Name: " << d_meshFileName << endl;
            cout << "There is no geometry file (GetnZone))!" << endl;
#ifndef ARIES_HAVE_MPI
            exit(EXIT_FAILURE);
#else
            MPI_Abort(MPI_COMM_WORLD,1);
            MPI_Finalize();
#endif
        }
  
        switch (d_meshFileType)
        {
            case MeshFile_SU2:
                /*--- Read SU2 mesh file ---*/
                for (iLine = 0; iLine < nLine ; iLine++)
                {
                    getline (meshFile, text_line);
                    
                    /*--- Search for the "NDIM" keyword to see if there are multiple Zones ---*/
                    position = text_line.find ("NDIME=",0);
                    if (position != string::npos)
                    {
                        text_line.erase (0,6);
                        nDim = atoi(text_line.c_str());
                    }
                }
                break;
      
            case MeshFile_CGNS:
#ifdef ARIES_HAVE_CGNS
                /*--- Local variables which are needed when calling the CGNS mid-level API. ---*/
                int fn, nbases = 0, nzones = 0, file_type;
                int cell_dim = 0, phys_dim = 0;
                char basename[CGNS_STRING_SIZE];
      
                /*--- Check whether the supplied file is truly a CGNS file. ---*/
                if ( cg_is_cgns(val_mesh_filename.c_str(), &file_type) != CG_OK )
                {
                    printf( "\n\n   !!! Error !!!\n" );
                    printf( " %s is not a CGNS file.\n", val_mesh_filename.c_str());
                    printf( " Now exiting...\n\n");
                    exit(EXIT_FAILURE);
                }
      
                /*--- Open the CGNS file for reading. The value of fn returned
                  is the specific index number for this file and will be
                  repeatedly used in the function calls. ---*/
                if (cg_open(val_mesh_filename.c_str(), CG_MODE_READ, &fn))
                    cg_error_exit();
      
                /*--- Get the number of databases. This is the highest node
                  in the CGNS heirarchy. ---*/
                if (cg_nbases(fn, &nbases))
                    cg_error_exit();
      
                /*--- Check if there is more than one database. Throw an
                  error if there is because this reader can currently
                  only handle one database. ---*/
                if ( nbases > 1 )
                {
                    printf("\n\n   !!! Error !!!\n" );
                    printf("CGNS reader currently incapable of handling more than 1 database.");
                    printf("Now exiting...\n\n");
                    exit(EXIT_FAILURE);
                }
      
                /*--- Read the databases. Note that the indexing starts at 1. ---*/
                for ( int i = 1; i <= nbases; i++ )
                {
                    if (cg_base_read(fn, i, basename, &cell_dim, &phys_dim))
                        cg_error_exit();
        
                    /*--- Get the number of zones for this base. ---*/
                    if (cg_nzones(fn, i, &nzones))
                        cg_error_exit();
        
                }
      
                /*
                 *  Set the problem dimension as read from the CGNS file
                 */
                nDim = cell_dim;
#endif
                break;
            default:
                break;
        }
        
        meshFile.close();
        return (unsigned short) nDim;
    }
    
    unsigned short ProcData::GetNumZone()
    {
        string text_line, Marker_Tag;
        ifstream meshFile;
        short nZone = 1; // Default value
        unsigned short iLine, nLine = 10;
        string::size_type position;

        //Open grid file
        meshFile.open(d_meshFileName, ios::in);
        if (meshFile.fail())
        {
            cout << "Mesh File Name: " << d_meshFileName << endl;
            cout << "There is no geometry file (GetnZone))!" << endl;
#ifndef ARIES_HAVE_MPI
            exit(EXIT_FAILURE);
#else
            MPI_Abort(MPI_COMM_WORLD,1);
            MPI_Finalize();
#endif
        }
  
        /*--- Search the mesh file for the 'NZONE' keyword. ---*/
        switch (d_meshFileType)
        {
            case MeshFile_SU2:
                /*--- Read the SU2 mesh file ---*/
                for (iLine = 0; iLine < nLine ; iLine++)
                {
                    getline (meshFile, text_line);
                    /*--- Search for the "NZONE" keyword to see if there are multiple Zones ---*/
                    position = text_line.find ("NZONE=",0);
                    if (position != string::npos)
                    {
                        text_line.erase (0,6);
                        nZone = atoi(text_line.c_str());
                    }
                }
      
                break;
            default:
                break;
        }
  
        /*--- For time spectral integration, nZones = nTimeInstances. ---*/
        if (d_unsteadyType == TIME_SPECTRAL)
        {
            nZone = d_numTimeInstance;
        }
  
        return (unsigned short) nZone;
    }
}



