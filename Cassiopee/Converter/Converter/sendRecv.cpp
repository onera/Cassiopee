/*    
    Copyright 2013-2025 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
// send/recv C

# include "converter.h"

#ifdef _MPI
# include "mpi.h"
# include "mpi4py/mpi4py.h"
#else
#define MPI_Comm void
#endif

using namespace std;

//==================================================
//           GETNAME : PyObject -> string
void getName(PyObject* obj, char* name)
{
    if (PyString_Check(obj)){
        //printf("[i] str python2\n");
        strcpy(name, PyString_AsString(obj));
    }
#if PY_VERSION_HEX >= 0x03000000
    if (PyUnicode_Check(obj)){
        //printf("[i] str python3\n");
        strcpy(name, PyBytes_AsString(PyUnicode_AsUTF8String(obj))); 
    }
#endif
}
//==================================================
//==================================================

//==================================================
//                C LAYER iSend
PyObject* K_CONVERTER::iSend(PyObject* self, PyObject* args)
{    
    // Recuperation des donnees
    PyObject* datas;
    PyObject* mpi4pyCom;
    E_Int oppNode;
    E_Int rank;
    if (!PYPARSETUPLE_(args, O_ II_ O_, &datas, &oppNode, &rank, &mpi4pyCom)) return NULL;


    // Recuperation du communicateur
#ifdef _MPI
    void* pt_comm = (void*)&(((PyMPICommObject*)mpi4pyCom)->ob_mpi);
    MPI_Comm comm = *((MPI_Comm*) pt_comm);
#endif
    
    // GESTION DU CAS ISEND(None)
    Py_INCREF(Py_None);
    if (datas==Py_None)
    {
        char* bufToSend = new char[1];
        *bufToSend = 'n';

        // Envoi du buf presque vide
#ifdef _MPI
        MPI_Request request;
        MPI_Isend(bufToSend, 1, MPI_CHAR, oppNode, 0, comm, &request);
        MPI_Request* crossRequest = new MPI_Request[1];
        *crossRequest = request;
#else
        void* crossRequest = NULL;
#endif
        
        PyObject* hook;
        void** packet = new void* [2]; // 2 pour request, bufToSend 
        packet[0] = crossRequest; // requete MPI
        packet[1] = bufToSend;
        #if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
        hook = PyCObject_FromVoidPtr(packet, NULL);
        #else
        hook = PyCapsule_New(packet, NULL, NULL);
        #endif

        return hook;
    }


    // Type de la donnee : premier ou deuxieme sendRecv
    E_Int dataType = 0;

    E_Int sizeDatas = PyList_Size(datas);

    // Pointeur de pointeurs
    char** bigBuf = new char* [sizeDatas];
    
    // Tableaux du nombre d'octets 
    std::vector<E_Int> tabOctets(sizeDatas);
    
    // Nombre d octets du type E_Int
    constexpr E_Int intSize = sizeof(E_Int);
    
    // Var des noms de zones
    char zoneName[256];
    char zoneDName[256];

    for (E_Int nData=0; nData<sizeDatas; nData++) // parcours des listes dans datas
    {
        // Recuperation de la liste des donnees
        PyObject* listData = PyList_GetItem(datas, nData);
        E_Int sizeList = PyList_Size(listData);
        if (sizeList == 6) // premier sendRecv de transfers
        {
            // TYPE 1 de donnees : [[zname, znamed, indicesI, XI, YI, ZI]]
            dataType = 1;

            // Nom de la zone a interpoler
            PyObject* zone = PyList_GetItem(listData, 0);
            getName(zone, zoneName);
            E_Int size_zoneName = strlen(zoneName);

            // Nom de la zone donneuse
            PyObject* zoneD = PyList_GetItem(listData, 1);
            getName(zoneD, zoneDName);
            E_Int size_zoneDName = strlen(zoneDName);

            // Indices des points à interpoler
            PyObject* PyIndices = PyList_GetItem(listData, 2);
            E_Int* indices; E_Int nIndices; E_Int nfld;
            K_NUMPY::getFromNumpyArray(PyIndices, indices, nIndices, nfld, true);

            // XCoordinates des points a interpoler
            PyObject* PyXCoord = PyList_GetItem(listData, 3);
            E_Float* xCoords; E_Int nXCoords; E_Int nfldx;
            K_NUMPY::getFromNumpyArray(PyXCoord, xCoords, nXCoords, nfldx, true);
            // YCoordinates des points a interpoler
            PyObject* PyYCoord = PyList_GetItem(listData, 4);
            E_Float* yCoords; E_Int nYCoords; E_Int nfldy;
            K_NUMPY::getFromNumpyArray(PyYCoord, yCoords, nYCoords, nfldy, true);
            // ZCoordinates des points a interpoler
            PyObject* PyZCoord = PyList_GetItem(listData, 5);
            E_Float* zCoords; E_Int nZCoords; E_Int nfldz;
            K_NUMPY::getFromNumpyArray(PyZCoord, zCoords, nZCoords, nfldz, true);
            
            // Calcul du nombre d'octets nécessaires :
            //  - nom de la zone : type, taille, nom
            //  - nom de la zone donneuse : type, taille, nom
            //  - indices des pts a interpo : type, taille, indices
            //  - coordonnées X des points a interpo : type, taille, X
            //  - coordonnées Y des points a interpo : type, taille, Y
            //  - coordonnées Z des points a interpo : type, taille, Z
            E_Int nOctets = 6*1 + 6*intSize + size_zoneName + size_zoneDName + nIndices*intSize + nXCoords*8*3;
            tabOctets[nData] = nOctets;
            // 6 char pour types + 6 entiers pour 6 tailles + n char nom x2 + nIndices entiers + nCoords floats *3 (3 coords) 

            // Initialisation du buffer
            char* buf = new char[nOctets];  

            // Placement dans le bigBuf
            bigBuf[nData] = buf; 

            // Placement nom de la zone
            (*buf) = 'c'                    ; buf+=1;
            (*(E_Int*)buf) = size_zoneName  ; buf+=intSize;
            for (E_Int i = 0; i < size_zoneName; i++) buf[i] = zoneName[i];
            buf += size_zoneName;

            // Placement nom de la zoneD
            (*buf) = 'c'                    ; buf+=1;
            (*(E_Int*)buf) = size_zoneDName ; buf+=intSize;
            for (E_Int i = 0; i < size_zoneDName; i++) buf[i] = zoneDName[i];
            buf += size_zoneDName;

            // Placement des indices
            (*buf) = 'i'                    ; buf+=1;
            (*(E_Int*)buf) = nIndices       ; buf+=intSize;
            E_Int* buf_indices = (E_Int*) buf;
            buf += intSize*nIndices;

            // Placement des X
            (*buf) = 'f'                    ; buf+=1;
            (*(E_Int*)buf) = nXCoords       ; buf+=intSize;
            E_Float* buf_X = (E_Float*) buf;
            buf += 8*nXCoords;
            
            // Placement des Y
            (*buf) = 'f'                    ; buf+=1;
            (*(E_Int*)buf) = nYCoords       ; buf+=intSize;
            E_Float* buf_Y = (E_Float*) buf;
            buf += 8*nYCoords;

            // Placement des Z
            (*buf) = 'f'                    ; buf+=1;
            (*(E_Int*)buf) = nZCoords       ; buf+=intSize;
            E_Float* buf_Z = (E_Float*) buf;
            buf += 8*nZCoords;

            // Placement parallele
            #pragma omp parallel
            {
                #pragma omp for
                for (E_Int i=0; i<nIndices; i++) { buf_indices[i] = indices[i];}

                #pragma omp for
                for (E_Int i=0; i<nXCoords; i++) { buf_X[i] = xCoords[i];}

                #pragma omp for
                for (E_Int i=0; i<nYCoords; i++) { buf_Y[i] = yCoords[i];}

                #pragma omp for
                for (E_Int i=0; i<nZCoords; i++) { buf_Z[i] = zCoords[i];}
            }

            // DECREF
            Py_DECREF(PyIndices);
            Py_DECREF(PyXCoord);
            Py_DECREF(PyYCoord);
            Py_DECREF(PyZCoord);
            
        }
        else if (sizeList == 3) // deuxieme sendRecv de _transfer
        {
            // Verification si 2e sendRecv de transfer ou si sendRecv de setInterpData
            if (PyString_Check(PyList_GetItem(listData, 0)))
            {
                dataType=2;
            }
#if PY_VERSION_HEX >= 0x03000000
            if (PyUnicode_Check(PyList_GetItem(listData, 0)))
            {
                dataType=2;
            }
#endif
            else { dataType=3; }
      
            if (dataType == 2)
            {
                // TYPE 2 de donnees : [[zrcvname,indicesR,fields]]
                
                // Nom de la zone recv
                PyObject* zone = PyList_GetItem(listData, 0);
                getName(zone, zoneName);
                E_Int size_zoneName = strlen(zoneName);
                // Indices des points a interpoler
                PyObject* PyIndices = PyList_GetItem(listData, 1);
                E_Int* indices; E_Int nIndices; E_Int nfld;
                K_NUMPY::getFromNumpyArray(PyIndices, indices, nIndices, nfld, true);

                // Fields a interpoler
                PyObject* PyFields = PyList_GetItem(listData, 2);
                char* fieldNames; FldArrayF* interpFields; E_Int nPts; E_Int nj, nk;
                FldArrayI* cn; char* eltType;
                E_Int oka = K_ARRAY::getFromArray(PyFields, fieldNames, 
                                                 interpFields, nPts, nj, nk, cn, eltType, true);
                E_Int nFlds = interpFields->getNfld();
                E_Float* fields = interpFields->begin();
                E_Int size_fieldNames = strlen(fieldNames);
                // Calcul du nombre d'octets necessaires :
                //  - nom de la zone 
                //  - indices des pts a interpo : type, taille, indices
                //  - valeurs des fields : type, taille, 6 fields
                E_Int nOctets = intSize*1 + 5*intSize + size_zoneName + size_fieldNames + nIndices*intSize + nPts*nFlds*8;
                tabOctets[nData] = nOctets;
                // 4 char pour types + 4 entiers pour 4 tailles + n char nom + nIndices entiers + nCoords floats *3 (3 coords) 

                // Initialisation du buffer
                char* buf = new char[nOctets]; 

                // Placement dans le bigBuf
                bigBuf[nData] = buf;

                // Placement nom de la zone
                (*buf) = 'c'                    ; buf += 1;
                (*(E_Int*)buf) = size_zoneName  ; buf += intSize;
                for (E_Int i = 0; i < size_zoneName; i++) buf[i] = zoneName[i];
                buf += size_zoneName;

                // Placement des indices
                (*buf) = 'i'                    ; buf += 1;
                (*(E_Int*)buf) = nIndices       ; buf += intSize;
                E_Int* buf_indices = (E_Int*) buf;
                buf += intSize*nIndices;
                
                // Placement nom des fields
                (*buf) = 'c'                     ; buf += 1;
                (*(E_Int*)buf) = size_fieldNames ; buf += intSize;
                for (E_Int i = 0; i < size_fieldNames; i++) buf[i] = fieldNames[i];
                buf += size_fieldNames;

                // Placement des fields
                (*buf) = 'f'                    ; buf += 1;
                (*(E_Int*)buf) = nPts           ; buf += intSize;
                (*(E_Int*)buf) = nFlds          ; buf += intSize;
                E_Float* buf_flds = (E_Float*) buf;
                buf += 8*nPts*nFlds;
                
                // Placement parallele
                #pragma omp parallel
                {
                    #pragma omp for
                    for (E_Int i=0; i<nIndices;   i++) { buf_indices[i] = indices[i];}

                    #pragma omp for
                    for (E_Int i=0; i<nPts*nFlds; i++) { buf_flds[i] = fields[i];}
                }

                // releaseshared
                Py_DECREF(PyIndices);
                RELEASESHAREDB(oka, PyFields, interpFields, cn);
            }
            else if (dataType == 3)
            {
                // TYPE 3 de donnees : [[int, fields, indices]]
                
                // Nom de la zone recv
                PyObject* PyVar = PyList_GetItem(listData, 0);
                E_Int var = PyLong_AsLong(PyVar);

                // Fields a interpoler
                PyObject* PyFields = PyList_GetItem(listData, 1);
                char* fieldNames; FldArrayF* interpFields; E_Int nPts; E_Int nj, nk;
                FldArrayI* cn; char* eltType;
                E_Int oka = K_ARRAY::getFromArray(PyFields, fieldNames, 
                                                 interpFields, nPts, nj, nk, cn, eltType, true);
                E_Int nFlds = interpFields->getNfld();
                E_Float* fields = interpFields->begin();
                E_Int size_fieldNames = strlen(fieldNames);

                // Indices des points a interpoler
                PyObject* PyIndices = PyList_GetItem(listData, 2);
                E_Int* indices; E_Int nIndices; E_Int nfld;
                K_NUMPY::getFromNumpyArray(PyIndices, indices, nIndices, nfld, true);

                // Calcul du nombre d'octets necessaires :
                //  - entier
                //  - indices des pts a interpo : type, taille, indices
                //  - valeurs des fields : type, taille, 6 fields
                E_Int nOctets = intSize*1 + 5*intSize + size_fieldNames + nIndices*intSize + nPts*nFlds*8;
                tabOctets[nData] = nOctets;
                // 4 char pour types + 4 entiers pour 4 tailles + n char nom + nIndices entiers + nCoords floats *3 (3 coords) 

                // Initialisation du buffer
                char* buf = new char[nOctets]; 

                // Placement dans le bigBuf
                bigBuf[nData] = buf;

                // Placement de l'entier
                (*buf) = 'i'          ; buf+=1;
                (*(E_Int*)buf) = var  ; buf+=intSize;

                // Placement nom des fields
                (*buf) = 'c'                     ; buf+=1;
                (*(E_Int*)buf) = size_fieldNames ; buf+=intSize;
                for (E_Int i = 0; i < size_fieldNames; i++) buf[i] = fieldNames[i];
                buf += size_fieldNames;

                // Placement des fields
                (*buf) = 'f'                    ; buf+=1;
                (*(E_Int*)buf) = nPts           ; buf+=intSize;
                (*(E_Int*)buf) = nFlds          ; buf+=intSize;
                E_Float* buf_flds = (E_Float*) buf;
                buf += 8*nPts*nFlds;

                // Placement des indices
                (*buf) = 'i'                    ; buf+=1;
                (*(E_Int*)buf) = nIndices       ; buf+=intSize;
                E_Int* buf_indices = (E_Int*) buf;
                buf += intSize*nIndices;

                // Placement parallele
                #pragma omp parallel
                {
                    #pragma omp for
                    for (E_Int i=0; i<nPts*nFlds; i++) { buf_flds[i] = fields[i];}
                    
                    #pragma omp for
                    for (E_Int i=0; i<nIndices;   i++) { buf_indices[i] = indices[i];}
                }
                
                // releaseshared
                Py_DECREF(PyIndices);
                RELEASESHAREDB(oka, PyFields, interpFields, cn);
            }
            else
            {
                printf("[" SF_D_ "][ERROR] Bad dataType (=" SF_D_ ")(= 2 or 3 normally)\n", rank, dataType); fflush(stdout);
                Py_INCREF(Py_None);
                return Py_None;
            }
        }
        else 
        {
            printf("[" SF_D_ "][ERROR] size of list = " SF_D_ " (= 3 or 6 normally)\n", rank, sizeList); fflush(stdout);
        }
    }

    // Calcul de la taille du buffer final
    E_Int nOctetsTot = 0;
    for (E_Int i=0; i < sizeDatas; i++)
    {
        nOctetsTot += tabOctets[i];
    }
    // sizeData entiers en plus pour connaitre la taille des mini buffers
    nOctetsTot += intSize*sizeDatas;
    // un entier en plus pour le nombre de data
    nOctetsTot += intSize;
    // un entier en plus pour le dataType
    nOctetsTot += intSize;

    // def du tableau a envoyer
    char* bufToSend = new char[nOctetsTot];
    char* initBufToSend = bufToSend;

    (*(E_Int*)bufToSend) = sizeDatas; bufToSend += intSize; // nombre de donnees
    (*(E_Int*)bufToSend) = dataType; bufToSend += intSize; // type de la donnees (1 ou 2)

    // Version by Clement
    #pragma omp parallel
    {
        for (E_Int i=0; i < sizeDatas; i++)
        {
            #pragma omp single
            { 
                // Placement du nb d'octets avant le buf
                (*(E_Int*)bufToSend) = tabOctets[i]; bufToSend += intSize;
            }

            char* buf = bigBuf[i];

            // Placement du buf
            #pragma omp for  
            for (E_Int j=0; j<tabOctets[i]; j++)
            {
                bufToSend[j] = buf[j];
            } 

            #pragma omp single
            {
                bufToSend += tabOctets[i];
            } 
        }
    }

    // Rewritten by CB : split explicite de la boucle sizeDatas
    /*
    E_Int nthreads = __NUMTHREADS__;
    std::vector<E_Int> ptrOctets(sizeDatas);
    ptrOctets[0] = 8;
    for (E_Int i = 1; i < sizeDatas; i++) ptrOctets[i] = ptrOctets[i-1]+intSize+tabOctets[i];

    #pragma omp parallel
    {
        E_Int ithread = __CURRENT_THREAD__;
        E_Int sizeDataLoc = sizeDatas / nthreads;
        for (E_Int i = sizeDataLoc*ithread; i < std::min(sizeDataLoc*(ithread+1), sizeDatas); i++)
        {
            char* bf = bufToSend+ptrOctets[i];
            char* buf = bigBuf[i];

            (*(E_Int*)bf) = tabOctets[i];
            for (E_Int j=0; j < tabOctets[i]; j++)
            {
                bf[j+intSize] = buf[j];
            }
        }
    }
    */
   
    // Envoi des donnees
#ifdef _MPI
    MPI_Request request;
    MPI_Isend(initBufToSend, nOctetsTot, MPI_CHAR, oppNode, 0, comm, &request);
#endif

    // Suppression des pointeurs (sauf celui de l'envoi)
    for (E_Int i=0; i < sizeDatas; i++)
    {
        delete [] bigBuf[i];
    }
    delete [] bigBuf;

#ifdef _MPI
    MPI_Request* crossRequest = new MPI_Request[1];
    *crossRequest = request;
#else
    void* crossRequest = NULL;
#endif

    PyObject* hook;
    void** packet = new void* [3]; // 4 pour request, bigBuf, nombre ptr dans bigBuf, initBufToSend 
    packet[0] = crossRequest; // requete MPI
    packet[1] = initBufToSend;
    #if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
    hook = PyCObject_FromVoidPtr(packet, NULL);
    #else
    hook = PyCapsule_New(packet, NULL, NULL);
    #endif

    return hook;
}
//==================================================
//==================================================



//==================================================
//                 C LAYER waitAll
PyObject* K_CONVERTER::waitAll(PyObject* self, PyObject* args)
{
    // Recuperation des donnees
    PyObject* reqs;
    if (!PYPARSETUPLE_(args, O_, &reqs)) return NULL;

    E_Int nRequests = PyList_Size(reqs);
    
    for (E_Int n=0; n<nRequests; n++)
    {
        PyObject* hook = PyList_GetItem(reqs, n);

        // recupere le hook
        void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
        packet = (void**) PyCObject_AsVoidPtr(hook);
#else
        packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
        
#ifdef _MPI
        MPI_Request* request = (MPI_Request*)packet[0];
        MPI_Wait(request, MPI_STATUS_IGNORE);
        // Suppression du ptr de la requete
        delete [] request;
#endif

        char* bufToSend = (char*)packet[1]; 

        // Suppression du buffer d'envoi
        delete [] bufToSend;
        // suppresion du paquet
        delete [] packet;
    }

    Py_INCREF(Py_None);
    return Py_None;
}
//==================================================
//==================================================




//==================================================
//                 C LAYER iRecv
PyObject* K_CONVERTER::recv(PyObject* self, PyObject* args)
{
    // Recuperation des donnees
    PyObject* mpi4pyCom;
    E_Int node;
    E_Int rank;
    if (!PYPARSETUPLE_(args, II_ O_, &node, &rank, &mpi4pyCom)) return NULL;

    // Recuperation du communicateur
#ifdef _MPI
    void* pt_comm = (void*)&(((PyMPICommObject*)mpi4pyCom)->ob_mpi);
    MPI_Comm comm = *((MPI_Comm*) pt_comm);
#endif

    // get MPI status
#ifdef _MPI
    MPI_Status status;
    MPI_Probe(node, 0, comm, &status);
#endif

    // get number of elements
    int nOctetsTot=0; // flaw!!
#ifdef _MPI
    MPI_Get_count(&status, MPI_CHAR, &nOctetsTot);
#endif

    // Definition des pointeurs
    char* recvBuf = new char[nOctetsTot];
    char* initRecvBuf = recvBuf;
    
    // Nombre d octets du type E_Int
    constexpr E_Int intSize = sizeof(E_Int);

    // reception du buffer
#ifdef _MPI
    MPI_Recv(recvBuf, nOctetsTot, MPI_CHAR, node, 0, comm, &status);
#endif

    // GESTION DU CAS ISEND(None)
    if (nOctetsTot == 1)
    {
        if (*recvBuf == 'n')
        {
            Py_INCREF(Py_None);
            return Py_None;  
        }
    }

    E_Int nOctets = 0;
    
    // recuperation du nombre de data
    E_Int* intRecvBuf = (E_Int*) initRecvBuf;
    
    E_Int sizeData = intRecvBuf[0]; recvBuf += intSize;
    // recuperation du type de donnees
    E_Int dataType = intRecvBuf[1]; recvBuf += intSize;

    // def de la liste des datas
    PyObject* datas = PyList_New(sizeData);

    if ((dataType != 1) && (dataType != 2) && (dataType != 3))
    {
        printf("Error: recv: unknown data type (return None).\n"); fflush(stdout);
        delete [] initRecvBuf;
        Py_INCREF(Py_None);
        return Py_None;
    } 

    for (E_Int nData=0; nData<sizeData; nData++)
    {
        intRecvBuf = (E_Int*) recvBuf;

        nOctets = intRecvBuf[0]; recvBuf+=intSize;

        // def buffer pour une data
        char* buf = new char[nOctets];
        char* initBuf = buf;

        // remplissage du buffer pour la data
        #pragma omp parallel for
        for (E_Int j=0; j<nOctets; j++)
        {
            buf[j] = recvBuf[j];
        }

        // Remplissage de la case nData de la liste selon dataType
        if (dataType == 1)
        {
            char typeData[1]; // type de ce qui vient dans buf
            E_Int size = 0; // taille de ce qui vient dans buf

            // Nom de la zone
            (*typeData)      = *buf      ; buf+=1;
            E_Int* intBuf    = (E_Int*) buf;
            size             = intBuf[0] ; buf+=intSize;
            char* zoneName   = new char[size+1];
            for (E_Int k=0; k<size; k++) { zoneName[k]   = buf[k]; }
            zoneName[size]='\0'; buf += size;

            // Nom de la zone donneuse
            (*typeData)       = *buf      ; buf+=1;
            intBuf            = (E_Int*) buf;
            size              = intBuf[0] ; buf+=intSize;
            char* zoneDName   = new char[size+1];
            for (E_Int k=0; k<size; k++) { zoneDName[k]   = buf[k]; }
            zoneDName[size]='\0'; buf += size;

            // Tableau des indices
            (*typeData)         = *buf    ; buf+=1; if ((*typeData)!='i'){printf("[" SF_D_ "][RECV] Probleme de type pour indices (!=integer)\n", rank); fflush(stdout);} ;
            intBuf              = (E_Int*) buf;
            size                = intBuf[0]; buf+=intSize;
            E_Int* indices      = new E_Int[size];
            E_Int* indicesBuf   = (E_Int*) buf;
            E_Int npts = size;
            buf += size*intSize;

            // Tableau des X
            (*typeData)         = *buf;      buf+=1; if ((*typeData)!='f'){printf("[" SF_D_ "][RECV] Probleme de type pour X (!=float)\n", rank); fflush(stdout);};
            intBuf              = (E_Int*) buf;
            size                = intBuf[0]; buf+=intSize;
            E_Float* xCoords    = new E_Float[size];
            E_Float* xBuf       = (E_Float*) buf;
            buf+=size*8;

            // Tableau des Y
            (*typeData)         = *buf;      buf+=1; if ((*typeData)!='f'){printf("[" SF_D_ "][RECV] Probleme de type pour Y (!=float)\n", rank); fflush(stdout);};
            intBuf              = (E_Int*) buf;
            size                = intBuf[0]; buf+=intSize;
            E_Float* yCoords    = new E_Float[size];
            E_Float* yBuf       = (E_Float*) buf;
            buf += size*8;

            // Tableau des Z
            (*typeData)         = *buf;      buf+=1; if ((*typeData)!='f'){printf("[" SF_D_ "][RECV] Probleme de type pour Z (!=float)\n", rank); fflush(stdout);};
            intBuf              = (E_Int*) buf;
            size                = intBuf[0]; buf+=intSize;
            E_Float* zCoords    = new E_Float[size];
            E_Float* zBuf       = (E_Float*) buf;
            buf += size*8;
            
            // Remplissage parallele
            #pragma omp parallel
            {
                #pragma omp for
                for (E_Int i=0; i<size; i++) { indices[i] = indicesBuf[i];}

                #pragma omp for
                for (E_Int i=0; i<size; i++) { xCoords[i] = xBuf[i];}

                #pragma omp for
                for (E_Int i=0; i<size; i++) { yCoords[i] = yBuf[i];}
            
                #pragma omp for
                for (E_Int i=0; i<size; i++) { zCoords[i] = zBuf[i];}
            }

            // Transformation des donnees C en donnees Python
            PyObject* PyZoneName = Py_BuildValue("s", zoneName);
            PyObject* PyZoneDName = Py_BuildValue("s", zoneDName);
            PyObject* PyIndices =  K_NUMPY::buildNumpyArray(indices, npts, 1, 1);
            PyObject* PyXCoords =  K_NUMPY::buildNumpyArray(xCoords, npts, 1, 1);
            PyObject* PyYCoords =  K_NUMPY::buildNumpyArray(yCoords, npts, 1, 1);
            PyObject* PyZCoords =  K_NUMPY::buildNumpyArray(zCoords, npts, 1, 1);
            
            // Liste finale
            PyObject* dataToFill = PyList_New(6);
            // passage des numpy dans la liste
            PyList_SET_ITEM(dataToFill, 0, PyZoneName);
            PyList_SET_ITEM(dataToFill, 1, PyZoneDName);
            PyList_SET_ITEM(dataToFill, 2, PyIndices);
            PyList_SET_ITEM(dataToFill, 3, PyXCoords);
            PyList_SET_ITEM(dataToFill, 4, PyYCoords);
            PyList_SET_ITEM(dataToFill, 5, PyZCoords);

            // liste dans datas
            PyList_SET_ITEM(datas, nData, dataToFill);

            delete [] zoneName; delete [] zoneDName;
            delete [] indices; delete [] xCoords; delete [] yCoords; delete [] zCoords;
        }
        else if (dataType == 2)
        {
            char typeData[1]; // type de ce qui vient dans buf
            E_Int size = 0; // taille de ce qui vient dans buf

            // Nom de la zone
            (*typeData)      = *buf      ; buf+=1; if ((*typeData)!='c'){printf("[" SF_D_ "][RECV] Probleme de type pour zoneName (!=char)\n", rank); fflush(stdout);};
            E_Int* intBuf    = (E_Int*) buf;
            size             = intBuf[0] ; buf+=intSize;
            char* zoneName   = new char[size+1];
            for (E_Int k=0; k<size; k++) { zoneName[k]   = buf[k]; }
            zoneName[size]='\0'; buf+=size;

            // Tableau des indices
            (*typeData)         = *buf    ; buf+=1; if ((*typeData)!='i'){printf("[" SF_D_ "][RECV] Probleme de type pour indices (!=integer)\n", rank); fflush(stdout);} ;
            intBuf              = (E_Int*) buf;
            size                = intBuf[0]; buf+=intSize;
            E_Int* indices      = new E_Int[size];
            E_Int* indicesBuf   = (E_Int*) buf;
            E_Int npts = size;
            PyObject* PyNpts = PyLong_FromLong(size);
            buf += size*intSize;

            // Nom des fields
            (*typeData)          = *buf      ; buf+=1; if ((*typeData)!='c'){printf("[" SF_D_ "][RECV] Probleme de type pour nameFields (!=char)\n", rank); fflush(stdout);};
            intBuf               = (E_Int*) buf;
            size                 = intBuf[0] ; buf+=intSize;
            char* fieldNames     = new char[size+1];
            for (E_Int k=0; k<size; k++) { fieldNames[k]   = buf[k]; }
            fieldNames[size]='\0'; buf+=size;

            // Tableau des fields
            (*typeData)         = *buf;      buf+=1; if ((*typeData)!='f'){printf("[" SF_D_ "][RECV] Probleme de type pour X (!=float)\n", rank); fflush(stdout);};
            intBuf              = (E_Int*) buf;
            size                = intBuf[0]; buf+=intSize;
            E_Int nFlds         = intBuf[1]; buf+=intSize;
            E_Float* fields     = new E_Float[size*nFlds];
            E_Float* floatBuf   = (E_Float*) buf;
            buf += size*nFlds*8;

            // Remplissage parallele
            #pragma omp parallel
            {
                #pragma omp for
                for (E_Int i=0; i<npts; i++)        { indices[i] = indicesBuf[i];}
                
                #pragma omp for
                for (E_Int i=0; i<size*nFlds; i++)  { fields[i] = floatBuf[i];}
            }

            // Transformation des données C en données Python
            PyObject* PyZoneName = Py_BuildValue("s", zoneName);
            PyObject* PyIndices =  K_NUMPY::buildNumpyArray(indices, npts, 1, 1);
            PyObject* PyFieldNames = Py_BuildValue("s", fieldNames);
            PyObject* PyFields =  K_NUMPY::buildNumpyArray(fields, npts, nFlds, 0);

            // Liste finale
            PyObject* dataToFill = PyList_New(3);
            PyObject* dataToFill2 = PyList_New(5);

            // passage des champs dans la liste
            PyList_SET_ITEM(dataToFill2, 0, PyFieldNames);
            PyList_SET_ITEM(dataToFill2, 1, PyFields);
            PyList_SET_ITEM(dataToFill2, 2, PyNpts);
            PyList_SET_ITEM(dataToFill2, 3, PyLong_FromLong(1));
            PyList_SET_ITEM(dataToFill2, 4, PyLong_FromLong(1));

            // passage des arguments dans la liste
            PyList_SET_ITEM(dataToFill, 0, PyZoneName);
            PyList_SET_ITEM(dataToFill, 1, PyIndices);
            PyList_SET_ITEM(dataToFill, 2, dataToFill2);

            // liste dans datas
            PyList_SET_ITEM(datas, nData, dataToFill);

            delete [] zoneName; delete [] fieldNames;
            delete [] indices; delete [] fields;

        }
        else if (dataType == 3)
        {
            char typeData[1]; // type de ce qui vient dans buf
            E_Int size = 0; // taille de ce qui vient dans buf

            // Entier
            (*typeData)      = *buf      ; buf+=1; if ((*typeData)!='i'){printf("[" SF_D_ "][RECV] Probleme de type pour l'entier (!=int)\n", rank); fflush(stdout);};
            E_Int* intBuf    = (E_Int*) buf;
            E_Int var        = intBuf[0] ; buf+=intSize;

            // Nom des fields
            (*typeData)          = *buf      ; buf+=1; if ((*typeData)!='c'){printf("[" SF_D_ "][RECV] Probleme de type pour nameFields (!=char)\n", rank); fflush(stdout);};
            intBuf               = (E_Int*) buf;
            size                 = intBuf[0] ; buf+=intSize;
            char* fieldNames     = new char[size+1];
            for (E_Int k=0; k<size; k++) { fieldNames[k] = buf[k]; }
            fieldNames[size] = '\0'; buf += size;

            // Tableau des fields
            (*typeData)         = *buf;      buf+=1; if ((*typeData)!='f'){printf("[" SF_D_ "][RECV] Probleme de type pour X (!=float)\n", rank); fflush(stdout);};
            intBuf              = (E_Int*) buf;
            size                = intBuf[0]; buf+=intSize;
            E_Int nFlds         = intBuf[1]; buf+=intSize;
            E_Float* fields     = new E_Float[size*nFlds];
            E_Float* floatBuf   = (E_Float*) buf;
            buf += size*nFlds*8;

            // Tableau des indices
            (*typeData)         = *buf    ; buf+=1; if ((*typeData)!='i'){printf("[" SF_D_ "][RECV] Probleme de type pour indices (!=integer)\n", rank); fflush(stdout);} ;
            intBuf              = (E_Int*) buf;
            E_Int npts          = intBuf[0]; buf+=intSize;
            E_Int* indices      = new E_Int[npts];
            E_Int* indicesBuf   = (E_Int*) buf;
            PyObject* PyNpts = PyLong_FromLong(npts);
            buf += size*intSize;
            
            // Remplissage parallele
            #pragma omp parallel
            {
                #pragma omp for
                for (E_Int i=0; i<size*nFlds; i++)  { fields[i] = floatBuf[i];}
                
                #pragma omp for
                for (E_Int i=0; i<npts; i++)        { indices[i] = indicesBuf[i];}
            }

            // Transformation des données C en données Python
            PyObject* PyVar = PyLong_FromLong(var);
            PyObject* PyIndices =  K_NUMPY::buildNumpyArray(indices, npts, 1, 1);
            PyObject* PyFieldNames = Py_BuildValue("s", fieldNames);
            if (size != npts) { nFlds=size; size=npts; } // forced a cause des fields a 1 point
            PyObject* PyFields =  K_NUMPY::buildNumpyArray(fields, size, nFlds, 0);

            // Liste finale
            PyObject* dataToFill = PyList_New(3);
            PyObject* dataToFill2 = PyList_New(5);

            // passage des champs dans la liste
            PyList_SET_ITEM(dataToFill2, 0, PyFieldNames);
            PyList_SET_ITEM(dataToFill2, 1, PyFields);
            PyList_SET_ITEM(dataToFill2, 2, PyNpts);
            PyList_SET_ITEM(dataToFill2, 3, PyLong_FromLong(1));
            PyList_SET_ITEM(dataToFill2, 4, PyLong_FromLong(1));

            // passage des arguments dans la liste
            PyList_SET_ITEM(dataToFill, 0, PyVar);
            PyList_SET_ITEM(dataToFill, 1, dataToFill2);
            PyList_SET_ITEM(dataToFill, 2, PyIndices);
            
            // liste dans datas
            PyList_SET_ITEM(datas, nData, dataToFill);

            delete [] indices; delete [] fields; delete [] fieldNames;
        }

        delete [] initBuf;

        recvBuf += nOctets; 
    }

    delete [] initRecvBuf;

    
    // Retour des datas
    return datas;
}
//==================================================
//==================================================
