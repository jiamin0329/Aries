/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for cgns mesh manager.
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */


#ifndef ARIES_MESHREADERCGNS_HPP
#define ARIES_MESHREADERCGNS_HPP

#include <vector>
#include <string>

using namespace std;

namespace ARIES
{
    class MeshReaderCGNS: public MeshReader
    {
    public:
        MeshReaderCGNS();

        ~MeshReaderCGNS();


        bool ReadMesh(IProcData* procData);




    private:


        



    };

   
}

#endif
 








