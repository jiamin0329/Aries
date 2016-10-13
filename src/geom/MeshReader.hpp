/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for mesh reader
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */


#ifndef ARIES_MESHREADER_HPP
#define ARIES_MESHREADER_HPP


#include "MeshData.hpp"

#include <string>

using namespace std;

namespace ARIES
{
    class MeshReader
    {
    public:
        MeshReader();
        MeshReader(MeshData *meshData);
        MeshReader();

        virtual ~MeshReader();


        virtual bool ReadMesh (IProcData* procData ) = 0;
    private:
        string d_fileName;
        

        long *Global_to_Local_Point;				/*!< \brief Global-local indexation for the points. */
        long *Local_to_Global_Point;				/*!< \brief Local-global indexation for the points. */
        unsigned short *Local_to_Global_Marker;	/*!< \brief Local to Global marker. */
        unsigned short *Global_to_Local_Marker;	/*!< \brief Global to Local marker. */
        unsigned long *adj_counter; /*!< \brief Adjacency counter. */
        unsigned long **adjacent_elem; /*!< \brief Adjacency element list. */



        

        
        


    };
}

#endif
 








